
#' FUNCTION TO COUNT TRIPLETS AGGREGATED ACROSS ALL SEQUENCES
#'
#' This function takes triplet counts for each locus (as a list) and aggregates them across all loci,
#' returning a global count for all triplets.
#'
#' @param triplet_counts A list of triplet counts, where each element corresponds to a locus.
#'   Each element is a named vector or table containing triplet counts for that locus.
#' @return A named vector (as a table) representing the aggregated counts of all triplets across all loci.
#'
#' @examples
#' loci <- c("chr1:1000000-1000100", "chr2:2000000-2000100")
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' sequence_df <- extract_sequences(loci, BSgenome.Hsapiens.UCSC.hg38)
#' print(sequence_df)
#' triplet_counts <- count_triplets_per_locus(sequence_df)
#' print(triplet_counts)
#'
#' global_triplet_counts <- count_triplets_aggregated(triplet_counts)
#' print(global_triplet_counts)
#' @export
count_triplets_aggregated <- function(triplet_counts) {
  # Find all possible triplets
  all_triplets <- unique(unlist(
    lapply(triplet_counts, function(x) {
      names(x)  # Extract triplet names
    })
  ))

  # Initialize a global count vector with zeros
  global_triplet_counts <- rep(0, length(all_triplets))
  names(global_triplet_counts) <- all_triplets  # Assign triplet names

  # Loop through each locus and sum counts
  for (locus in names(triplet_counts)) {
    # Get the current locus counts
    current_counts <- triplet_counts[[locus]]

    # Standardize the counts to match all possible triplets
    current_counts_vector <- rep(0, length(all_triplets))
    names(current_counts_vector) <- all_triplets
    current_counts_vector[names(current_counts)] <- current_counts  # Update counts

    # Add to global counts
    global_triplet_counts <- global_triplet_counts + current_counts_vector
  }

  # Convert to a table and return
  #global_triplet_counts_table <- as.table(global_triplet_counts)

  return(global_triplet_counts)
}




