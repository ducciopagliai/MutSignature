
#' FUNCTION TO COUNT TRIPLETS FOR EACH LOCUS
#'
#' This function takes a data frame of DNA sequences from genomic loci and counts the triplets (3-mers) present in each sequence.
#' In addition, if the central base of a triplet is either A or G, its reverse complement is considered.
#'
#' @param sequence_df A data frame containing two columns:
#'   \itemize{
#'     \item locus: The genomic loci as input (e.g. "chr1:1000-2000").
#'     \item sequence: The DNA sequence corresponding to each locus.
#'   }
#' @return A list where each element corresponds to a locus and contains a table of triplet counts for that locus.
#'
#' @examples
#' loci <- c("chr1:1000000-1000100", "chr2:2000000-2000100")
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' sequence_df <- extract_sequences(loci, BSgenome.Hsapiens.UCSC.hg38)
#' print(sequence_df)
#'
#' triplet_counts <- count_triplets_per_locus(sequence_df)
#' print(triplet_counts)
#' @export
count_triplets_per_locus <- function(sequence_df) {

  # Check if the dataframe is not empty
  if (nrow(sequence_df) == 0) {
    warning("The input data frame is empty.")
    return(list())
  }

  triplet_counts <- list()

  # FOR loop for each sequence
  for (i in seq_len(nrow(sequence_df))) {
    sequence <- sequence_df$sequence[i]

    # Check the sequence length
    if (!is.character(sequence) || nchar(sequence) < 3) {
      warning(sprintf("Sequence at row %d is too short or invalid and will be skipped.", i))
      next
    }

    # Triplets filtering
    triplets <- sapply(1:(nchar(sequence) - 2), function(i) {
      substr(sequence, i, i + 2)
    })

    # Retrieve reverse complements of genomic loci if the central base is an A or G
    initial_triplets <- sapply(triplets, function(triplet) {
      central_base <- substr(triplet, 2, 2)
      if (central_base %in% c("A", "G")) {
        reverse_complement <- chartr("ACGT", "TGCA",
                                     rev(strsplit(triplet, NULL)[[1]]) |>
                                       paste(collapse = ""))
        return(reverse_complement)
      } else {
        return(triplet)
      }
    })

    sequence_triplet_counts <- table(initial_triplets)
    triplet_counts[[sequence_df$locus[i]]] <- sequence_triplet_counts
  }

  return(triplet_counts)
}


