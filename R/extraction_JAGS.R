## Little Library to proper extract information from JAGS output ##

###########################################################################
# select_posterior_FM -------------------------------------------------
#' Extract selected JAGS output
#'
#' Useful function that extracts the required information in the form of a dataframe.
#' Use this to extract a specific beta sample, mu, tau ad tau_eps JAGS samplings from
#' the posterior.
#'
#' @param output The JAGS output.
#' @param index Required in param="beta", select all the betas related to the i-th
#' sample. Default vaule is \code{NULL}.
#'
#' @return A dataframe containing the required information as columns.

#' @export
select_posterior_FM <- function(output, param, index=NULL){

  # Necessary conversion of JAGS output
  #library('stringr')
  #library('coda')
  global_df <- data.frame(as.data.frame(coda::as.array.mcmc.list(output)))
  row.names(global_df) <- row.names(data.frame(as.data.frame(coda::as.array.mcmc.list(output))))

  # Check if you provided an index and you didn't want to select beta samples
  if(!is.null(index) & param != "beta"){
    warning("Index is discarded as not necessary. You can remove it.")
  }

  # Pattern selection accrding to input parameter
  if(param=="beta"){
    # Check if the index is missing
    if(is.null(index))
      stop("Index i is required to extract betas of the i'ih sample")
    else
      pattern <- paste0(param,"\\.",index,"\\.[0-9]+")
  } else if (param=="mu" || param=="tau") {
    pattern <- paste0(param,"\\.[0-9]+")
  } else if (param=="tau_eps") {
    pattern <- paste0(param)
  } else {
    stop("Wrong parameter inserted.")
  }

  # Construciton of the output dataframe
  result_df <- data.frame(global_df[,which(stringr::str_detect(names(global_df), pattern))])

  # CHECK
  # cat("Index: ", index, "dims = ", dim(result_df),"\n")

  # Rename in case of tau_eps variable
  if (param=="tau_eps")
    names(result_df) <- param

  # Check if dataframe is not empty
  if(dim(result_df)[2] == 0)
    warning("Empty output, maybe something went wrong.")

  # Return the result
  return(result_df)
}

###########################################################################

###########################################################################
# select_all_beta --------------------------------------------------------
#' Extract the whole beta posterior sampling.
#'
#' Function built specifically to group all the beta posterior samplings coming
#' from JAGS by sample unit. It creates a list in which the [[i]]-th element is
#' a dataframe that collects the beta samplings related to the i-th curves.
#'
#' @param output The JAGS output.
#' @param n.curves The number of curves for which we've executed the sampling
#' strategy with JAGS.
#'
#' @return A list in which the [[i]]-th element is the dataframe.
#' referring to the i-th curve.

#' @export
select_all_beta <- function(output, n.curves){
  # Set up the output to build
  list_result <- list(NULL)

  # Fill the list element-wise
  for (i in 1:n.curves) {
    list_result[[i]] <- select_posterior_FM(output, "beta", i)
  }
  return(list_result)
}

###########################################################################

###########################################################################
# select_beta_per_iteration -----------------------------------------------
#' Select a specific sample iteration of beta
#'
#' Function that extracts the sampling of beta related to the i-th curve at
#' a specific iteration.
#'
#' @param output The JAGS output.
#' @param index This indicates the sample unit which needs to be considered.
#' @param n.iter Iteration at which you want to extract the beta sampling.
#'
#' @return A 1xp dataframe with the desired draw of the simulation.

#' @export
select_beta_per_iteration <- function(output, index, n.iter){
  return(select_posterior_FM(output, "beta", index)[n.iter,])
}

###########################################################################
