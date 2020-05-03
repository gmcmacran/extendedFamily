#' Additional Binomial Links for glm Models
#'
#' @param link name of link function. One of loglog (Default: loglog)
#' @details
#' family is a generic function with methods for classes "glm" and "lm".
#'
#' @return An object of class "family" (which has a concise print method). This is a list with elements
#'
#'
#' \itemize{
#'   \item family: character: the family name.
#'   \item link: character: the link name.
#'   \item linkfun: function: the link.
#'   \item linkinv: function: the inverse of the link function.
#'   \item variance: function: the variance as a function of the mean.
#'   \item dev.resids function giving the deviance for each observation as a function of (y, mu, wt),
#'         used by the residuals method when computing deviance residuals.
#'   \item aic: function giving the AIC value if appropriate (but NA for the quasi- families).
#'         More precisely, this function returns -2 ll + 2 s, where ll is the log-likelihood and s
#'         is the number of estimated scale parameters. Note that the penalty term for the location
#'         parameters (typically the “regression coefficients”) is added elsewhere, e.g., in glm.fit(),
#'         or AIC(), see the AIC example in glm. See logLik for the assumptions made about the dispersion
#'         parameter.
#'   \item initialize: expression. This needs to set up whatever data objects are needed for the family
#'         as well as n (needed for AIC in the binomial family) and mustart (See
#'         \href{https://stat.ethz.ch/R-manual/R-devel/RHOME/library/stats/html/glm.html}{glm})
#'   \item validmu: logical function. Returns TRUE if a mean vector mu is within the domain of variance.
#'   \item valideta: logical function. Returns TRUE if a linear predictor eta is within the domain of linkinv.
#'         \href{https://stat.ethz.ch/R-manual/R-devel/RHOME/library/stats/html/simulate.html}{simulate}.
#'         It will normally return a matrix with nsim columns and one row for each fitted value, but
#'         it can also return a list of length nsim. Clearly this will be missing for ‘quasi-’ families.
#'         }
#'
#' @examples
#' library(stats)
#' library(extendedFamily)
#'
#' # loglog example
#' data(heart)
#' model <- glm(formula = death ~ anterior + hcabg +
#'                                 kk2 + kk3 + kk4 + age2 + age3 + age4,
#'                       data = heart,
#'                       family = binomialEF(link = "loglog"))
#'
#' @export
binomialEF <- function (link = "loglog")
{
  assertthat::assert_that(length(link) == 1, msg = "Argument link should have length 1.")
  assertthat::assert_that(is.character(link),msg = "Argument link should be a character.")
  assertthat::assert_that(link %in% c("loglog"),msg = "Argument link should be 'loglog'.")

  linktemp <- link
  switch(link,
         "loglog" = {
           linkfun <- function(mu) {-log(-log(mu))}
           linkinv <- function(eta) {
             mu <- pmin(exp(-exp(-eta)),1 - .Machine$double.eps)
             mu <- pmax(mu,.Machine$double.eps)
             return(mu)
           }
           mu.eta <- function(eta) {
             eta <- pmax(eta, -700)
             mu <- exp(-exp(-eta)) * exp(-eta)
             mu <- pmax(mu, .Machine$double.eps)
             return(mu)
           }
           valideta <- function(eta) {TRUE}
         }
  )

  copyMe <- stats::binomial(link = "cloglog")
  variance <- copyMe$variance
  validmu <- copyMe$validmu
  dev.resids <- copyMe$dev.resids
  aic <- copyMe$aic
  initialize <- copyMe$initialize
  simfun <- copyMe$simfun

  structure(list(family = "binomial",
                 link = linktemp,
                 linkfun = linkfun,
                 linkinv = linkinv,
                 variance = variance,
                 dev.resids = dev.resids,
                 aic = aic,
                 mu.eta = mu.eta,
                 initialize = initialize,
                 validmu = validmu,
                 valideta = valideta,
                 simulate = simfun),
            class = "family")
}
