library(Matrix)      # For handling sparse matrices
library(igraph)      # For co-authorship network analysis
library(tm)

#The MADstat dataset is publicly available at https://github.com/ZhengTracyKe/MADStat

# author_list <- readLines("author_name.txt")
# load("AuthorPaperInfo.RData")
# load("BibtexInfo.RData")



####################### Raw data filtering and splitting #######################
# Time window selection
AuPapMat_filtered <- subset(AuPapMat, year >= 2003 & year <= 2012)

# Journal classification
# 1. AuPapMat_4: contains only papers published in Bka, JRSSB, JASA, or AoS
AuPapMat_4 <- subset(AuPapMat_filtered, journal %in% c("Bka", "JRSSB", "JASA", "AoS"))
# 2. AuPapMat_other: contains papers from other journals
AuPapMat_other <- subset(AuPapMat_filtered, !(journal %in% c("Bka", "JRSSB", "JASA", "AoS")))

# Identifier isolation
idxAu_4 <- unique(AuPapMat_4$idxAu)
idxPap_4 <- unique(AuPapMat_4$idxPap)
idxAu_other <- unique(AuPapMat_other$idxAu)
idxPap_other <- unique(AuPapMat_other$idxPap)

# Data splitting
author_list_4 <- author_list[idxAu_4]
author_list_other <- author_list[idxAu_other]
paper_4 <- paper[idxPap_4, ]
paper_other <- paper[idxPap_other, ]
paper_author_4 <- paper_author[idxPap_4]
paper_author_other <- paper_author[idxPap_other]

####################### Further data filtering (for target journals) #######################

author_mapping_4 <- setNames(1:length(idxAu_4), idxAu_4)
n_authors_4 <- length(idxAu_4)
n_papers_4 <- length(idxPap_4)
rows_4 <- integer()
cols_4 <- integer()

# Binary association matrix construction
for (k in seq_along(paper_author_4)) {
  authors <- paper_author_4[[k]]$id
  if (length(authors) == 0) next  
  row_indices <- author_mapping_4[as.character(authors)]  # Convert author IDs to row indices
  rows_4 <- c(rows_4, row_indices)
  cols_4 <- c(cols_4, rep(k, length(row_indices)))
}
adj_matrix_4 <- sparseMatrix(
  i = rows_4,
  j = cols_4,
  x = 1,
  dims = c(n_authors_4, n_papers_4),
  dimnames = list(
    author = author_list[idxAu_4],  # Row names: author names
    paper = idxPap_4                # Column names: paper IDs
  )
)

####################### Further data filtering (for source journals) #######################

author_mapping_other <- setNames(1:length(idxAu_other), idxAu_other)
n_authors_other <- length(idxAu_other)
n_papers_other <- length(idxPap_other)
rows_other <- integer()
cols_other <- integer()

for (k in seq_along(paper_author_other)) {
  authors <- paper_author_other[[k]]$id
  if (length(authors) == 0) next
  row_indices <- author_mapping_other[as.character(authors)]
  rows_other <- c(rows_other, row_indices)
  cols_other <- c(cols_other, rep(k, length(row_indices)))
}
adj_matrix_other <- sparseMatrix(
  i = rows_other,
  j = cols_other,
  x = 1,
  dims = c(n_authors_other, n_papers_other),
  dimnames = list(
    author = author_list[idxAu_other],
    paper = idxPap_other
  )
)

####################### Text processing ####################### 

# Define processing function

create_author_term_avg <- function(adj_matrix, paper_subset, author_list_subset) {
 # 1: Filter authors with at least 2 papers
  author_paper_counts <- rowSums(adj_matrix)
  valid_authors <- which(author_paper_counts >= 2)
  adj_matrix_filtered <- adj_matrix[valid_authors, ]
  author_list_filtered <- author_list_subset[valid_authors]

 # 2: Process paper titles
  titles <- tolower(paper_subset$title)
  titles <- removePunctuation(titles)
  titles <- removeWords(titles, stopwords("en"))
  titles <- strsplit(titles, "\\s+")
  titles <- lapply(titles, function(x) x[x != ""])
  
 # 3: Build term frequency matrix
  term_freq <- table(unlist(titles))
  term_freq <- term_freq[names(term_freq) != ""]
  term_freq <- term_freq[term_freq >= 10]  # Keep terms appearing at least 10 times
  
 # 4: Build author-term matrix
  top_terms <- names(sort(term_freq, decreasing = TRUE))
  author_term_matrix <- Matrix(0, 
                               nrow = nrow(adj_matrix_filtered),
                               ncol = length(top_terms),
                               dimnames = list(author_list_filtered, top_terms))
  
  for (j in 1:ncol(adj_matrix_filtered)) {
    if (j %% 1000 == 0) message("Processing paper ", j)
    paper_terms <- intersect(titles[[j]], top_terms)
    if (length(paper_terms) > 0) {
      authors <- which(adj_matrix_filtered[, j] > 0)
      author_term_matrix[authors, paper_terms] <- author_term_matrix[authors, paper_terms] + 1
    }
  }
  
 # 5: Compute average term frequency
  author_total_papers <- rowSums(adj_matrix_filtered)
  author_term_avg <- author_term_matrix / author_total_papers
  
 # 6: Build coauthorship network and extract largest connected component
  coauthor_network <- tcrossprod(adj_matrix_filtered)
  diag(coauthor_network) <- 0
  g <- graph_from_adjacency_matrix(coauthor_network, mode = "undirected", weighted = TRUE)
  comp <- components(g)
  largest_component <- which(comp$membership == which.max(comp$csize))
  
  # Return results
  list(
    author_term_avg = author_term_avg[largest_component, ],
    coauthor_network = coauthor_network[largest_component, largest_component]
  )
}

# Run the function on both subsets
result_4 <- create_author_term_avg(
  adj_matrix = adj_matrix_4,
  paper_subset = paper_4,
  author_list_subset = author_list_4
)

result_other <- create_author_term_avg(
  adj_matrix = adj_matrix_other,
  paper_subset = paper_other,
  author_list_subset = author_list_other
)

####################### Final results #######################
                   
author_term_matrix_4 <- result_4$author_term_avg
author_term_matrix_other <- result_other$author_term_avg

target_authors_4 <- rownames(author_term_matrix_4)
target_indices_4 <- match(target_authors_4, author_list_4)
coauthor_4 <- tcrossprod(adj_matrix_4)
sub_adj_4 <- coauthor_4[target_indices_4, target_indices_4]
coauthor_network_4 <- (sub_adj_4 > 0) * 1  # Binarize matrix while preserving sparsity
diag(coauthor_network_4) <- 0  # Remove diagonal

target_authors_other <- rownames(author_term_matrix_other)
target_indices_other <- match(target_authors_other, author_list_other)
coauthor_other <- tcrossprod(adj_matrix_other)
sub_adj_other <- coauthor_other[target_indices_other, target_indices_other]
coauthor_network_other <- (sub_adj_other > 0) * 1  # Binarize matrix while preserving sparsity
diag(coauthor_network_other) <- 0  # Remove diagonal


# The file 'DataForGNC-Plot-Combined.Rda' contains the top 300 terms selected by tf-idf scores, as described in the work by Li et al. (2020) on high-dimensional Gaussian graphical models for network-linked data. This dataset is publicly available in the GNC repository at https://github.com/tianxili/GNC/blob/master/GNC-lasso.R.
# load("DataForGNC-Plot-Combined.Rda")
# x <- X[, (1:300)]
# author_term_matrix_4<- author_term_matrix_4[, colnames(x)]
# author_term_matrix_other <- author_term_matrix_other[, colnames(x)]
