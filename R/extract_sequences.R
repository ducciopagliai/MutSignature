# Libraries
library(BSgenome)  # Provides access to genome sequences
library(GenomicRanges)  # Helps work with genomic regions, like chromosomes and their positions


# FUNCTION TO EXTRACT DNA SEQUENCES FROM GENOMIC LOCI
extract_sequences <- function(loci, genome) {

  # Verify that the genome is a BSgenome
  if (!inherits(genome, "BSgenome")) {
    stop("The genome must be a BSgenome object, e.g., BSgenome.Hsapiens.UCSC.hg38")
  }

  # Convert loci (e.g., "chr1:1000-2000") into a GRanges object
  genomic_ranges <- makeGRangesFromDataFrame(
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
  sequences <- getSeq(genome, genomic_ranges)

  # Store the result using a dataframe
  result <- data.frame(
    locus = loci,  # Original loci
    sequence = as.character(sequences)  # Sequences
  )

  # Return the result
  return(result)
}

