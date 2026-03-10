#' Resolve sampler control parameters
#'
#' Internal helper that merges user-supplied sampler controls with defaults.
#'
#' @param control Optional named list of control parameters.
#' @param has_Z Logical; whether a second-stage moderator model is present.
#'
#' @return A named list of resolved control parameters.
#' @keywords internal
resolve_sampler_control <- function(control = NULL, has_Z = FALSE) {

  defaults <- list(
    # gamma prior (only used when Z_block present)
    Sigma_gamma_prior = 10,
    mu_gamma_prior = NULL,

    # sigma^2 prior
    a_sigma_alpha_prior = 2,
    b_sigma_alpha_prior = 2,

    # tau prior
    a_sigma_tau_prior = 2,
    b_sigma_tau_prior = 2
  )

  if (is.null(control)) return(defaults)

  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0) {
    stop(
      sprintf(
        "Unknown control parameter(s): %s",
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  utils::modifyList(defaults, control)
}
