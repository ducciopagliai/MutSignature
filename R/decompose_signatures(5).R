
#' FUNCTION TO DECOMPOSE MUTATIONAL SIGNATURES
#'
#' This function performs decomposition of mutational signatures based on COSMIC counts.
#' It takes as input either a list of COSMIC counts per locus or aggregated COSMIC counts,
#' and determines the contribution of each mutational signature.
#'
#' @param cosmic_counts Either a list of named numeric vectors (COSMIC counts per locus)
#' or a single named numeric vector (aggregated COSMIC counts).
#' @param cosmic_path The path to the COSMIC signatures file in Alexandrov format
#' (e.g. COSMIC v3.4 SBS GRCh38).
#'
#' @return A list of decomposed mutational signatures: if the input is a list of COSMIC counts,
#' a list of results per locus, if input is aggregated COSMIC counts, a single result.
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
#' cosmic_path <- system.file("extdata", "COSMIC_v3.4_SBS_GRCh38 2.txt", package = "MutSignature")
#' cosmic_counts_per_locus <- map_to_cosmic(triplet_counts, cosmic_path)
#' print(cosmic_counts_per_locus)
#' cosmic_counts_aggregated <- map_to_cosmic(global_triplet_counts, cosmic_path)
#' print(cosmic_counts_aggregated)
#'
#' # Decompose per locus
#' result_per_locus <- decompose_signatures(cosmic_counts_per_locus, cosmic_path)
#' print(result_per_locus)
#'
#' # Decompose aggregated counts
#' result_aggregated <- decompose_signatures(cosmic_counts_aggregated, cosmic_path)
#' print(result_aggregated)
#' @export
decompose_signatures <- function(cosmic_counts, cosmic_path) {
  # Load COSMIC signatures
  signatures <- decompTumor2Sig::readAlexandrovSignatures(file = cosmic_path)

  # Check validity of signatures
  if (!decompTumor2Sig::isAlexandrovSet(signatures)) {
    stop("Mutational signatures are not correct.")
  }

  # Initialize the final result
  final_result <- list()

  # Case 1 --> Input is the result of count_triplets_per_locus
  if (is.list(cosmic_counts)) {
    for (locus in names(cosmic_counts)) {
      # Normalization
      counts <- cosmic_counts[[locus]]
      cosmic_counts_normalized <- counts / sum(counts)

      # Check the normalized counts
      if (!decompTumor2Sig::isAlexandrovSet(list(cosmic_counts_normalized))) {
        stop(sprintf("Normalized counts for locus '%s' are not valid for the Alexandrov format.", locus))
      }

      # DECOMPOSITION
      final_result[[locus]] <- decompTumor2Sig::decomposeTumorGenomes(
        genomes = list(cosmic_counts_normalized),
        signatures = signatures
      )
    }
  }
  # Case 2: --> Input is the result of count_triplets_aggregated
  else if (is.vector(cosmic_counts) && is.numeric(cosmic_counts)) {
    # Normalization
    cosmic_counts_normalized <- cosmic_counts / sum(cosmic_counts)

    # Check the normalized counts
    if (!decompTumor2Sig::isAlexandrovSet(list(cosmic_counts_normalized))) {
      stop("Normalized counts are not valid for the Alexandrov format.")
    }

    # DECOMPOSITION
    final_result <- decompTumor2Sig::decomposeTumorGenomes(
      genomes = list(cosmic_counts_normalized),
      signatures = signatures
    )
  } else {
    stop("Input cosmic_counts must be either a list of counts per locus or a single named vector of aggregated counts.")
  }

  return(final_result)
}






