
#' FUNCTION TO GENERATE POSSIBLE MUTATIONS FOR EACH TRIPLET
#'
#' This function generates all possible mutations for a given triplet in the COSMIC format,
#' by replacing the central base of the triplet with all other bases except itself.
#' This function
#'
#' @param triplet A character string representing a DNA triplet (e.g., "ACC").
#' @return A vector of all possible mutations in the COSMIC format.
#'
#' @keywords internal
generate_mutations <- function(triplet) {
  # Bases
  first_base <- substr(triplet, 1, 1)
  central_base <- substr(triplet, 2, 2)
  last_base <- substr(triplet, 3, 3)

  # Possible mutations
  possible_mutations <- setdiff(c("A", "C", "G", "T"), central_base)

  # Generate all mutations in the COSMIC format
  mutations <- sapply(possible_mutations, function(i) {
    paste0(first_base, "[", central_base, ">", i, "]", last_base)
  })

  return(mutations)
}

