
#' EXTRACT DNA SEQUENCES FROM GENOMIC LOCI
#'
#' This function extracts DNA sequences from a set of genomic loci provided as input,
#' using a reference genome object of class BSgenome.
#'
#' @param loci A character vector where each element represents a genomic locus in the format "chr:start-end".
#'             For example, "chr1:1000-2000".
#' @param genome A reference genome object of class BSgenome (e.g. BSgenome.Hsapiens.UCSC.hg38).
#' #' @return A data frame with the following columns:
#' \itemize{
#'   \item locus: The genomic loci given as input.
#'   \item sequence: The corresponding DNA sequences.
#' }
#' @examples
#' loci <- c("chr1:1000000-1000100", "chr2:2000000-2000100")
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' sequence_df <- extract_sequences(loci, BSgenome.Hsapiens.UCSC.hg38)
#' print(sequence_df)
#' @export
extract_sequences <- function(loci, genome) {

  # Loci check
  if (!is.character(loci)) {
    stop("The 'loci' argument must be a character vector.")
  }

  valid_loci <- grepl("^chr[0-9XY]+:[0-9]+-[0-9]+$", loci)
  if (!all(valid_loci)) {
    stop("Each element in the 'loci' argument must follow the format 'chr:start-end'.")
  }

  # BSgenome check
  if (!inherits(genome, "BSgenome")) {
    stop("The 'genome' argument must be an object of class 'BSgenome'.")
  }

  # Convert loci (e.g., "chr1:1000-2000") into a GRanges object
  genomic_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    data.frame(
      seqnames = sub(":.*", "", loci),  # Extract chr names
      ranges = sub(".*:", "", loci) |> strsplit("-") |>  # Split start and end
        lapply(as.numeric) |>  # Transform start and end to numeric values
        do.call(what = rbind)  # Combine the results into a matrix
    ),
    seqnames.field = "seqnames",
    start.field = "ranges.1",
    end.field = "ranges.2"
  )

  # Extract the sequences (corresponding to the GRanges intervals)
  sequences <- BSgenome::getSeq(genome, genomic_ranges)

  # Store the result using a dataframe
  sequence_df <- data.frame(
    locus = loci,  # Original loci
    sequence = as.character(sequences)  # Sequences
  )

  # Return the result
  return(sequence_df)
}


