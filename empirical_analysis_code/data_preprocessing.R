############################################################
# MADStat Dataset Preprocessing
#
# In this script, we describe the preprocessing procedure of
# the MADStat dataset used for empirical analysis.
#
# To ensure direct comparability with the results reported in
# Li et al. (2020), we strictly follow their preprocessing
# protocol, including temporal window selection and term
# extraction criteria.
#
# The preprocessing procedure uses four main data objects:
#
#   1. AuPapMat
#      Author-paper-journal-year association information
#
#   2. paper_author
#      Author information for each publication
#
#   3. paper
#      Publication metadata including titles
#
#   4. author_list
#      Author identification list
#
# The complete preprocessing procedure consists of:
#
# Step 1:
#   Raw data filtering and splitting
#   - Select publications within predefined time windows
#   - Separate target journals and source journals
#
# Step 2:
#   Author-paper network construction
#   - Construct binary author-paper adjacency matrices
#   - Remove authors with insufficient publications
#   - Extract the largest connected component
#
# Step 3:
#   Text preprocessing and feature construction
#   - Process publication titles
#   - Remove stopwords and rare terms
#   - Construct author-term frequency matrices
#   - Apply nonparanormal transformation
#
# The resulting objects are:
#
#   author_term_avg_4_list
#       Target-domain author-term feature matrices
#
#   author_term_avg_other_list
#       Source-domain author-term feature matrices
#
#   sub_adj_binary_4_list
#       Target-domain coauthor networks
#
#   sub_adj_binary_other_list
#       Source-domain coauthor networks
#
############################################################



############################################################
# Load required packages
############################################################

library(Matrix)
library(igraph)
library(tm)
library(huge)


#The MADstat dataset is publicly available at https://github.com/ZhengTracyKe/MADStat

# author_list <- readLines("author_name.txt")
# load("AuthorPaperInfo.RData")
# load("BibtexInfo.RData")





############################################################
# Function: Construct author-level textual features and
# coauthor networks
############################################################

create_author_term_avg <- function(
    adj_matrix,
    paper_subset,
    author_list_subset
) {
  
  library(Matrix)
  library(tm)
  library(igraph)
  
  
  ##########################################################
  # Step 1: Remove authors with only one publication
  ##########################################################
  
  author_paper_counts <- rowSums(adj_matrix)
  
  valid_authors <- which(author_paper_counts >= 2)
  
  adj_matrix_filtered <- adj_matrix[valid_authors, ]
  
  author_list_filtered <- author_list_subset[valid_authors]
  
  
  
  ##########################################################
  # Step 2: Text preprocessing
  # (No stemming is applied)
  ##########################################################
  
  titles <- tolower(paper_subset$title)
  
  titles <- removePunctuation(titles)
  
  titles <- removeNumbers(titles)
  
  titles <- removeWords(titles, stopwords("en"))
  
  
  # Tokenization
  
  titles <- strsplit(titles, "\\s+")
  
  
  # Remove empty tokens
  
  titles <- lapply(
    titles,
    function(x) x[x != ""]
  )
  
  
  
  ##########################################################
  # Step 3: Calculate document frequency
  ##########################################################
  
  term_doc_count <- table(
    unlist(lapply(titles, unique))
  )
  
  
  # Keep terms appearing in at least 10 papers
  
  top_terms <- names(
    term_doc_count[term_doc_count >= 10]
  )
  
  
  
  ##########################################################
  # Step 4: Construct author-term frequency matrix
  ##########################################################
  
  author_term_matrix <- Matrix(
    0,
    nrow = nrow(adj_matrix_filtered),
    ncol = length(top_terms),
    dimnames = list(
      author_list_filtered,
      top_terms
    ),
    sparse = TRUE
  )
  
  
  for (j in 1:ncol(adj_matrix_filtered)) {
    
    if (j %% 1000 == 0) {
      message("Processing paper ", j)
    }
    
    
    # Term frequency for the current paper
    
    term_counts <- table(titles[[j]])
    
    
    # Keep only selected terms
    
    valid_terms <- intersect(
      names(term_counts),
      top_terms
    )
    
    
    if (length(valid_terms) > 0) {
      
      authors <- which(
        adj_matrix_filtered[, j] > 0
      )
      
      
      for (term in valid_terms) {
        
        author_term_matrix[authors, term] <-
          author_term_matrix[authors, term] +
          term_counts[term]
        
      }
    }
  }
  
  
  
  ##########################################################
  # Step 5: Normalize by publication counts
  # nonparanormal transformation
  ##########################################################
  
  author_total_papers <- rowSums(
    adj_matrix_filtered
  )
  
  
  # Average term frequency per author
  
  author_term_avg <-
    author_term_matrix / author_total_papers
  
  
  
  # Apply nonparanormal transformation
  
  author_term_avg <- huge.npn(
    as.matrix(author_term_avg),
    npn.func = "truncation"
  )
  
  
  # Convert back to sparse matrix format
  
  author_term_avg <- Matrix(
    author_term_avg,
    sparse = TRUE
  )
  
  
  
  ##########################################################
  # Step 6: Construct coauthor network
  ##########################################################
  
  coauthor_network <- tcrossprod(
    adj_matrix_filtered
  )
  
  
  diag(coauthor_network) <- 0
  
  
  g <- graph_from_adjacency_matrix(
    coauthor_network,
    mode = "undirected",
    weighted = TRUE
  )
  
  
  
  ##########################################################
  # Step 7: Extract largest connected component
  ##########################################################
  
  comp <- components(g)
  
  largest_component_id <- which.max(
    comp$csize
  )
  
  largest_component <- which(
    comp$membership == largest_component_id
  )
  
  
  
  ##########################################################
  # Return processed data
  ##########################################################
  
  list(
    author_term_avg = 
      author_term_avg[largest_component, ],
    
    coauthor_network =
      coauthor_network[
        largest_component,
        largest_component
      ],
    
    authors =
      author_list_filtered[largest_component]
  )
}
############################################################
# Time periods
############################################################

time_periods <- list(
  "1981_2000" = 1981:2000,
  "2001_2010" = 2001:2010,
  "2011_2015" = 2011:2015
)



############################################################
# Initialize storage objects
############################################################

author_term_avg_4_list <- list()

author_term_avg_other_list <- list()

sub_adj_binary_4_list <- list()

sub_adj_binary_other_list <- list()



############################################################
# Process each time period
############################################################

for (period_name in names(time_periods)) {
  
  cat(
    "Processing time period:",
    period_name,
    "\n"
  )
  
  
  years <- time_periods[[period_name]]
  
  
  ##########################################################
  # Select publications within the current period
  ##########################################################
  
  AuPapMat_filtered <- subset(
    AuPapMat,
    year %in% years
  )
  
  
  
  ##########################################################
  # Split target journals and auxiliary journals
  ##########################################################
  
  AuPapMat_4 <- subset(
    AuPapMat_filtered,
    journal %in% c(
      "Bka",
      "JRSSB",
      "JASA",
      "AoS"
    )
  )
  
  
  AuPapMat_other <- subset(
    AuPapMat_filtered,
    !(journal %in% c(
      "Bka",
      "JRSSB",
      "JASA",
      "AoS"
    ))
  )
  
  
  
  ##########################################################
  # Extract author and paper indices
  ##########################################################
  
  idxAu_4 <- unique(
    AuPapMat_4$idxAu
  )
  
  idxPap_4 <- unique(
    AuPapMat_4$idxPap
  )
  
  
  idxAu_other <- unique(
    AuPapMat_other$idxAu
  )
  
  idxPap_other <- unique(
    AuPapMat_other$idxPap
  )
  
  
  
  ##########################################################
  # Extract author and paper information
  ##########################################################
  
  author_list_4 <- author_list[idxAu_4]
  
  author_list_other <- author_list[idxAu_other]
  
  
  paper_4 <- paper[idxPap_4, ]
  
  paper_other <- paper[idxPap_other, ]
  
  
  paper_author_4 <- paper_author[idxPap_4]
  
  paper_author_other <- paper_author[idxPap_other]
  
  
  
  ##########################################################
  # Construct author-paper adjacency matrix for target data
  ##########################################################
  
  author_mapping_4 <- setNames(
    1:length(idxAu_4),
    idxAu_4
  )
  
  
  rows_4 <- integer()
  
  cols_4 <- integer()
  
  
  for (k in seq_along(paper_author_4)) {
    
    authors <- paper_author_4[[k]]$id
    
    
    if (length(authors) == 0) {
      next
    }
    
    
    row_indices <- 
      author_mapping_4[
        as.character(authors)
      ]
    
    
    rows_4 <- c(
      rows_4,
      row_indices
    )
    
    
    cols_4 <- c(
      cols_4,
      rep(k, length(row_indices))
    )
  }
  
  
  adj_matrix_4 <- sparseMatrix(
    i = rows_4,
    j = cols_4,
    x = 1,
    dims = c(
      length(idxAu_4),
      length(idxPap_4)
    ),
    dimnames = list(
      author = author_list_4,
      paper = idxPap_4
    )
  )
  
  
  
  ##########################################################
  # Construct author-paper adjacency matrix for auxiliary data
  ##########################################################
  
  author_mapping_other <- setNames(
    1:length(idxAu_other),
    idxAu_other
  )
  
  
  rows_other <- integer()
  
  cols_other <- integer()
  
  
  for (k in seq_along(paper_author_other)) {
    
    authors <- paper_author_other[[k]]$id
    
    
    if (length(authors) == 0) {
      next
    }
    
    
    row_indices <-
      author_mapping_other[
        as.character(authors)
      ]
    
    
    rows_other <- c(
      rows_other,
      row_indices
    )
    
    
    cols_other <- c(
      cols_other,
      rep(k, length(row_indices))
    )
  }
  
  
  adj_matrix_other <- sparseMatrix(
    i = rows_other,
    j = cols_other,
    x = 1,
    dims = c(
      length(idxAu_other),
      length(idxPap_other)
    ),
    dimnames = list(
      author = author_list_other,
      paper = idxPap_other
    )
  )
  
  
  
  ##########################################################
  # Generate author textual features and networks
  ##########################################################
  
  result_4 <- create_author_term_avg(
    adj_matrix_4,
    paper_4,
    author_list_4
  )
  
  
  result_other <- create_author_term_avg(
    adj_matrix_other,
    paper_other,
    author_list_other
  )
  
  
  
  ##########################################################
  # Store author-term matrices
  ##########################################################
  
  author_term_avg_4_list[[period_name]] <-
    result_4$author_term_avg
  
  author_term_avg_other_list[[period_name]] <-
    result_other$author_term_avg
  
  
  
  ##########################################################
  # Extract binary coauthor networks
  ##########################################################
  
  target_authors_4 <-
    rownames(result_4$author_term_avg)
  
  
  target_indices_4 <-
    match(
      target_authors_4,
      author_list_4
    )
  
  
  coauthor_network_4 <-
    tcrossprod(adj_matrix_4)
  
  
  sub_adj_4 <-
    coauthor_network_4[
      target_indices_4,
      target_indices_4
    ]
  
  
  sub_adj_binary_4 <-
    (sub_adj_4 > 0) * 1
  
  
  diag(sub_adj_binary_4) <- 0
  
  
  sub_adj_binary_4_list[[period_name]] <-
    sub_adj_binary_4
  
  
  
  target_authors_other <-
    rownames(result_other$author_term_avg)
  
  
  target_indices_other <-
    match(
      target_authors_other,
      author_list_other
    )
  
  
  coauthor_network_other <-
    tcrossprod(adj_matrix_other)
  
  
  sub_adj_other <-
    coauthor_network_other[
      target_indices_other,
      target_indices_other
    ]
  
  
  sub_adj_binary_other <-
    (sub_adj_other > 0) * 1
  
  
  diag(sub_adj_binary_other) <- 0
  
  
  sub_adj_binary_other_list[[period_name]] <-
    sub_adj_binary_other
  
}



############################################################
# Align vocabulary across different time periods
############################################################

common_terms <- Reduce(
  intersect,
  lapply(
    author_term_avg_4_list,
    colnames
  )
)


cat(
  "Number of common terms:",
  length(common_terms),
  "\n"
)



############################################################
# Keep only common terms
############################################################

author_term_avg_4_list <-
  lapply(
    author_term_avg_4_list,
    function(mat) {
      
      mat[
        ,
        common_terms,
        drop = FALSE
      ]
    }
  )


author_term_avg_other_list <-
  lapply(
    author_term_avg_other_list,
    function(mat) {
      
      mat[
        ,
        common_terms,
        drop = FALSE
      ]
    }
  )



############################################################
# Save processed data
############################################################

save(
  author_term_avg_4_list,
  author_term_avg_other_list,
  sub_adj_binary_4_list,
  sub_adj_binary_other_list,
  common_terms,
  file = "processed_author_network_data.RData"
)