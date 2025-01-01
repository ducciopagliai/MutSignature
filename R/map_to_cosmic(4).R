
#' FUNCTION TO MAP TRIPLET COUNTS TO THE COSMIC FORMAT
#'
#' This function maps triplet counts, either per locus or aggregated, to the COSMIC mutational signatures format.
#'
#' @param triplet_counts A list of triplet counts per locus (from count_triplets_per_locus),
#'   or a named numeric vector of aggregated triplet counts (from count_triplets_aggregated).
#' @param cosmic_path The path to the COSMIC signatures file in Alexandrov format
#' (e.g. COSMIC v3.4 SBS GRCh38).
#' @return A list of COSMIC counts for each locus if the input is a list of triplet counts,
#' or a single vector of aggregated COSMIC counts if the input is a numeric vector.
#'
#' @examples
#' loci <- c("chr1:1000000-1000100", "chr2:2000000-2000100")
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' sequence_df <- extract_sequences(loci, BSgenome.Hsapiens.UCSC.hg38)
#' print(sequence_df)
#' triplet_counts <- count_triplets_per_locus(sequence_df)
#' print(triplet_counts)
#' global_triplet_counts <- count_triplets_aggregated(triplet_counts)
#' print(global_triplet_counts)
#'
#' cosmic_path <- system.file("extdata", "COSMIC_v3.4_SBS_GRCh38 2.txt", package = "MutSignature")
#'
#' cosmic_counts_per_locus <- map_to_cosmic(triplet_counts, cosmic_path)
#' print(cosmic_counts_per_locus)
#' cosmic_counts_aggregated <- map_to_cosmic(global_triplet_counts, cosmic_path)
#' print(cosmic_counts_aggregated)
#' @export
map_to_cosmic <- function(triplet_counts, cosmic_path) {
  # Check if the file exists
  if (!file.exists(cosmic_path)) {
    stop("The specified cosmic_path does not exist.")
  }

  # Validate file extension
  if (tools::file_ext(cosmic_path) != "txt" && tools::file_ext(cosmic_path) != "tsv") {
    stop("The specified file must be a tab-delimited text file (.txt or .tsv).")
  }

  # Load the COSMIC matrix
  cosmic_matrix <- read.table(cosmic_path, header = TRUE, sep = "\t", row.names = 1)

  # Check COSMIC matrix
  if (ncol(cosmic_matrix) == 0) {
    stop("The COSMIC file appears to have no data or the separator is incorrect.")
  }

  # Check input type
  if (is.list(triplet_counts)) {
    # Case 1: Input is a list of triplet counts (per locus)
    cosmic_counts_list <- lapply(triplet_counts, function(single_triplet_counts) {
      map_single_triplet_set(single_triplet_counts, cosmic_matrix)
    })
    return(cosmic_counts_list)
  } else if (is.vector(triplet_counts) && is.numeric(triplet_counts)) {
    # Case 2: Input is a named vector (aggregated counts)
    return(map_single_triplet_set(triplet_counts, cosmic_matrix))
  } else {
    stop("Input must be either a list of triplet counts or a single named vector.")
  }
}

#' MAP A SINGLE SET OF TRIPLET COUNTS TO THE COSMIC FORMAT
#'
#' This function maps a single set of triplet counts to the COSMIC mutational signatures format.
#'
#' @param single_triplet_counts A numeric vector of triplet counts.
#' @param cosmic_matrix A COSMIC mutational signatures matrix loaded as a data frame.
#'
#' @return A numeric vector of COSMIC counts corresponding to the mutations.
#'
#' @keywords internal
map_single_triplet_set <- function(single_triplet_counts, cosmic_matrix) {
  # Vector count of zeros
  cosmic_counts <- setNames(rep(0, nrow(cosmic_matrix)), rownames(cosmic_matrix))

  # Transform each triplet to its mutations and map to the COSMIC matrix
  for (triplet in names(single_triplet_counts)) {
    # Generate all possible mutations for the triplet
    mutations <- generate_mutations(triplet)

    # Map the mutations to the COSMIC matrix
    for (mutation in mutations) {
      if (mutation %in% rownames(cosmic_matrix)) {
        cosmic_counts[mutation] <- cosmic_counts[mutation] + single_triplet_counts[triplet]
      }
    }
  }

  return(cosmic_counts)
}



