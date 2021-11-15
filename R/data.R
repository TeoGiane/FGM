#' Absorbance mid-IR spectra of 983 fruit purées.
#'
#' This dataset contains the absorbance spectra of two different classes of fruit purées, evaluated
#' in 235 different wavelenghts. Among them, we are able to distinguish data coming from pure strawberry
#' purees and from adulterated purées, obtained through a mixture of different fruits. This dataset is
#' in the proper format for the creation of an FGM object. See \code{\link{FGM}} for details.
#'
#' @usage data("purees")
#'
#' @format A list that contains two elements. \strong{wavelengths} contains the above-mentioned
#' vector of wavelengths, while \strong{data} is a data frame with 983 rows and 236 variables:
#' \describe{
#'   \item{Group}{Factor variable identifying if the specturm comes from a strawberry purees or not.}
#'   \item{W_899.327}{Absorbance value registered at 899.327 nm. All other variables have similar meaning.}
#' }
#' @source \url{https://csr.quadram.ac.uk/example-datasets-for-download/}
"purees"

#' Absorbance mid-IR spectra of 351 strawberry purées.
#'
#' This is a shortcut to load only a part of the dataset described in \code{\link{purees}}.
#' This dataset contains only the 351 absorbance spectra of strawberry purées, measured at 235 different wavelenghts.
#' Differently from \code{\link{purees}}, this file contains only the data.frame of the measurments. The corresponding
#' wavelengths can be loaded using \code{\link{StrawberryWavelengths}}.
#'
#' @usage data("StrawberryPurees")
#'
"StrawberryPurees"

#' Grid of wavelengths where absorbance spectra are measured.
#'
#' This is a shortcut to load only a part of the dataset described in \code{\link{purees}}.
#' This dataset contains a vector of 235 elements representing the grid where the absorbance spectra availabe in \code{\link{StrawberryPurees}}
#' are measured.
#'
#' @usage data("StrawberryWavelengths")
#'
"StrawberryWavelengths"

#' Peterson's Graphs and Precision Matrices
#'
#' This dataset contains a list of four elements, which contains the graph and precision matrices
#' of dimension pxp suggested by Peterson et al. (2015) in order to perform simulations on graphical models.
#'
#' @usage data("Peterson_p16") or data("Peterson_p40")
#'
#' @format A list of four elements, each of the is a list containing:
#' \describe{
#'   \item{G}{The Peterson's graph suggested in the reference paper.}
#'   \item{K}{A pxp precision matrix iduced by G.}
#' }
#' @source \url{https://www.researchgate.net/publication/261020849}
"Peterson_p16"

#' Peterson's Graphs and Precision Matrices
#'
#' This dataset contains a list of four elements, which contains the graph and precision matrices
#' of dimension pxp suggested by Peterson et al. (2015) in order to perform simulations on graphical models.
#'
#' @usage data("Peterson_p16") or data("Peterson_p40")
#'
#' @format A list of four elements, each of the is a list containing:
#' \describe{
#'   \item{G}{The Peterson's graph suggested in the reference paper.}
#'   \item{K}{A pxp precision matrix iduced by G.}
#' }
#' @source \url{https://www.researchgate.net/publication/261020849}
"Peterson_p40"
