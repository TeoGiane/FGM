## Useful converters to proper manage BDgraph package output ##

###########################################################################
# str2adj -----------------------------------------------------------------
#' Conversion to adjacency matrix
#'
#' This is a converter from the default string of the sample.graphs vector
#' to its proper adjacency matrix, in the form of a binary matrix.
#'
#' @param p Size of the graph, i.e. the number of nodes.
#' @param graph.str Binary string to convert
#'
#' @return A pxp binary matrix representing the graph in input.

#' @export
str2adj <- function(p, graph.str){
  #library('stringr')
  parsed_vector <- as.numeric(unlist(stringr::str_extract_all(graph.str, "[0-9]")))
  result <- matrix(0, p, p)
  result[upper.tri(result)] <-  parsed_vector
  result <- result + t(result)
  return (result)
}

###########################################################################

###########################################################################
# adj2str -----------------------------------------------------------------
#' Coerce adjacency matrix to a string
#'
#' This is a converter from the the adjacency matrix that represent a graph
#' to the default string of the sample.graphs vector, representing the upper
#' triangular matrix shrinked into a string of length (p-1)(p-2)/2.
#'
#' @param adj.matrix The pxp binary matrix representing the graph.
#'
#' @return The converted string.

#' @export
adj2str <- function(adj.matrix){
  #library('stringr')
  result <- toString(adj.matrix[upper.tri(adj.matrix)])
  result <- stringr::str_remove_all(result, ", ")
  return(result)
}

###########################################################################

###########################################################################
# str2K -------------------------------------------------------------------
#' Conversion to precision matrix
#'
#' This is a converter from a comma-separated string representing the upper
#' triangular precision matrix (diagonal included) to the equivalent pxp positive
#' definite precision matrix. Useful since precision matrices are stored as strings
#' in order to reduce space.
#'
#' @param p Size of the graph, i.e. the number of nodes.
#' @param K.str The comma-separated string containing the values of the upper
#' triangular matrix (diagonal included)
#'
#' @return The converted pxp spd precision matrix.

#' @export
str2K <- function(p, K.str){
  #library('stringr')
  parsed_vector <- as.numeric(unlist(stringr::str_split(K.str, ",")))
  result <- matrix(0, p, p)
  result[upper.tri(result, diag = TRUE)] <- parsed_vector
  result <- (result + t(result)) - diag(diag(result))
  return (result)
}

###########################################################################
