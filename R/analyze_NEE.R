#' Analyze under NEE assumptions
#'
#' Analyze data assuming no early treatment effects
#'
#' @param data Data frame containing the following variables
#' \itemize{
#'    \item Z: indicator of treatment
#'    \item Y: indicator of outcome
#'    \item Y_tau: indicator of early outcome
#'    \item S_star: intermediate biomarker value
#'    \item R: indicator of measurement of intermediate biomarker
#' }
#' @param brange Numeric (2 x 1) vector containing the specified lower and upper bounds of the range for sensitivity parameter \ifelse{html}{\out{&#946;<sub>0</sub>}}{\eqn{\beta_0}}
#' @param design String describing the study design / sampling scheme used. This allows for estimation of sampling weights. Options include "full", "cc" (case-cohort), and "other". When "other" is chosen the weights argument must also be specified
#' @param weights Numeric (n x 1) vector containing pre-estimated sampling weights where n is the number of rows in `data`
#' @param contrast Contrast function for estimand. Options include "logRR", "Difference", and "VE"
#'
#' @return Returns list consisting of 6 vectors corresponding to the ignorance intervals and EUIs of CEP(1, 0), CEP(0, 0), and the difference CEP(1, 0) - CEP(0, 0)
#' @export
#'
#' @examples Z <- rbinom(500, 1, 0.5)
#' S_star <- rbinom(500, 1, 0.2)
#' R <- rep(1, 500)
#' Y_tau <- rep(0, 500)
#' Y <- rbinom(500, 1, 0.1)
#' df <- data.frame(Z, S_star, R, Y_tau, Y)
#' analyze_NEE(df, c(-0.5, 0.5), design = "full", contrast = "VE")
analyze_NEE <- function(data, brange = c(0, 0), design = "full",
                        weights = NULL, contrast = "logRR") {

  # Attach data
  data$S_star[is.na(data$S_star)] <- 0
  Z <- data$Z
  Y <- data$Y
  Y_tau <- data$Y_tau
  S_star <- data$S_star
  R <- data$R
  n <- dim(data)[1]

  # Set or estimate weights
  if (tolower(design) == "full") {
    W <- rep(1, n)
  }
  if (tolower(design) == "cc") {
    # estimate probability control has S_star observed: P(R = 1 | Y = 0)
    solve_pi <- function(x) { sum((1 - Y)*(R - x)) }
    pi_hat <- uniroot(solve_pi, c(0, 1))$root
    W <- 1 / pi_hat*(1 - Y)*R + Y
  }
  if (!is.null(weights)) {
    W <- weights
  }


  # Point estimation

  S_0 <- (1 - Y_tau)*(1 - S_star)
  S_1 <- (1 - Y_tau)*S_star

  # Estimate identifiable parameters risk_z, p(0, 0), p(1, 0), and p(1, 1)
  # risk_0 = P(Y(0) = 1 | Y^tau(0) = Y^tau(1) = 0)
  solve_risk0 <- function(x) { sum((1 - Y_tau)*(1 - Z)*(Y - x)) }
  risk_0 <- uniroot(solve_risk0, c(0, 1))$root

  # risk_1 = P(Y(1) = 1 | Y^tau(0) = Y^tau(1) = 0)
  solve_risk1 <- function(x) { sum((1 - Y_tau)*Z*(Y - x)) }
  risk_1 <- uniroot(solve_risk1, c(0, 1))$root

  # p(0, 0) = P(S^tau(1) = S^tau(0) = 0 | Y^tau(0) = Y^tau(1) = 0)
  # Recall Case CB: P(S^tau(0) = 0 | Y^tau(0) = Y^tau(1) = 0) = 1
  solve_p00 <- function(x) { sum((1 - Y_tau)*Z*(1 - S_star - x)*W) }
  p_00 <- uniroot(solve_p00, c(0, 1))$root
  p_10 <- 1 - p_00
  # p(1, 1) = 0 under NEE scenario

  # Estimate identifiable parameters risk_1(0, 0) and risk_1(1, 0)
  # risk_z(s1, s0) =
  # P(Y(z) = 1 | S^tau(1) = s1, S^tau(0) = s0, Y^tau(1) = Y^tau(0) = 0)
  solve_risk_1_00 <- function(x) { sum(Z*S_0*(Y - x)*W) }
  risk_1_00 <- uniroot(solve_risk_1_00, c(0, 1))$root
  solve_risk_1_10 <- function(x) { sum(Z*S_1*(Y-x)*W) }
  risk_1_10 <- uniroot(solve_risk_1_10, c(0, 1))$root

  # Estimate risk_0(0, 0) using SACE method
  # This is the only partially identifiable term under NEE assumptions
  solve_risk_0_00_min <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(min(brange))*(1 - x)/x )
    risk_0 - x*p_00 - risk_0_10*p_10
  }
  risk_0_00_first <- uniroot(solve_risk_0_00_min, c(0, 1))$root
  risk_0_10_first <- (risk_0 - risk_0_00_first*p_00) / p_10

  solve_risk_0_00_max <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(max(brange))*(1 - x)/x )
    risk_0 - x*p_00 - risk_0_10*p_10
  }
  risk_0_00_second <- uniroot(solve_risk_0_00_max, c(0, 1))$root
  risk_0_10_second <- (risk_0 - risk_0_00_second*p_00) / p_10

  # Calculate CEP(1,0) ignorance interval
  CEP_10_II_low <- min(h(risk_1_10, risk_0_10_first, contrast),
                       h(risk_1_10, risk_0_10_second, contrast))
  CEP_10_II_up <- max(h(risk_1_10, risk_0_10_first, contrast),
                      h(risk_1_10, risk_0_10_second, contrast))
  whichmin_10 <- which.min(c(h(risk_1_10, risk_0_10_first, contrast),
                             h(risk_1_10, risk_0_10_second, contrast)))
  whichmax_10 <- which.max(c(h(risk_1_10, risk_0_10_first, contrast),
                             h(risk_1_10, risk_0_10_second, contrast)))

  # Calculate CEP(0,0) ignorance interval
  CEP_00_II_low <- min(h(risk_1_00, risk_0_00_first, contrast),
                       h(risk_1_00, risk_0_00_second, contrast))
  CEP_00_II_up <- max(h(risk_1_00, risk_0_00_first, contrast),
                      h(risk_1_00, risk_0_00_second, contrast))
  whichmin_00 <- which.min(c(h(risk_1_00, risk_0_00_first, contrast),
                             h(risk_1_00, risk_0_00_second, contrast)))
  whichmax_00 <- which.max(c(h(risk_1_00, risk_0_00_first, contrast),
                             h(risk_1_00, risk_0_00_second, contrast)))

  # Calculate CEP(1,0) - CEP(0,0) ignorance interval
  CEP_diff_II_low <- min(h(risk_1_10, risk_0_10_first, contrast) -
                           h(risk_1_00, risk_0_00_first, contrast),
                         h(risk_1_10, risk_0_10_second, contrast) -
                           h(risk_1_00, risk_0_00_second, contrast))
  CEP_diff_II_up <- max(h(risk_1_10, risk_0_10_first, contrast) -
                          h(risk_1_00, risk_0_00_first, contrast),
                        h(risk_1_10, risk_0_10_second, contrast) -
                          h(risk_1_00, risk_0_00_second, contrast))
  whichmin_diff <- which.min(c(h(risk_1_10, risk_0_10_first, contrast) -
                                 h(risk_1_00, risk_0_00_first, contrast),
                               h(risk_1_10, risk_0_10_second, contrast) -
                                 h(risk_1_00, risk_0_00_second, contrast)))
  whichmax_diff <- which.max(c(h(risk_1_10, risk_0_10_first, contrast) -
                                 h(risk_1_00, risk_0_00_first, contrast),
                               h(risk_1_10, risk_0_10_second, contrast) -
                                 h(risk_1_00, risk_0_00_second, contrast)))




  # Variance estimation

  # Source M-estimation helper functions
  #source("compute.R")
  #source("utilities.R")

  # Source files with estimating equations
  #source("No Early Effects/low_eefun_NEE.R")
  #source("No Early Effects/up_eefun_NEE.R")

  # Define helper function for Imbens-Manski interval computation
  f_cstar <- function(c, low, up, maxsig) {
    pnorm(c + (up - low) / maxsig) - pnorm(-c) - 0.95
  }

  # Get point estimates needed for sandwich variance estimation
  # and choose appropriate stack of estimating equations
  thetahat <- as.numeric(c(risk_0, risk_1_10, p_00, risk_1_00,
                           risk_0_00_first, risk_0_10_first,
                           risk_0_00_second, risk_0_10_second,
                           CEP_00_II_low, CEP_00_II_up, CEP_10_II_low,
                           CEP_10_II_up, CEP_diff_II_low, CEP_diff_II_up))

  eefun <- eefun_NEE

  if (tolower(design) == "cc") {
    thetahat <- c(thetahat, pi_hat)
    eefun <- eefun_NEE_cc
  }

  # Convert data to wide form
  data$W <- W
  data_wide <- plyr::count(data, c('Y', 'Z', 'Y_tau', 'S_star', 'R', 'W'))
  data_wide$Group <- 1:nrow(data_wide)

  # Calculate variance estimate for the lower bound of the ignorance
  # intervals using 'geex' package source code
  mats <- compute_matrices(list(eeFUN = eefun,
                                splitdt = split(data_wide,
                                                f = data_wide$Group)),
                           numDeriv_options = list(method = 'simple'),
                           theta = thetahat, beta0range = brange,
                           contrast, whichmin_00, whichmax_00,
                           whichmin_10, whichmax_10, whichmin_diff,
                           whichmax_diff)

  A <- apply(simplify2array(Map(`*`, mats$A_i,
                                data_wide$freq)), 1:2, sum)
  B <- apply(simplify2array(Map(`*`, mats$B_i,
                                data_wide$freq)), 1:2, sum)

  Sigma <- solve(A) %*% B %*% t(solve(A))

  # Extract variances of interest from sandwich matrix
  var_CEP_00_II_low <- Sigma[9, 9]
  var_CEP_10_II_low <- Sigma[11, 11]
  var_CEP_diff_II_low <- Sigma[13, 13]


  # Same process for upper bound
  if (min(brange) == max(brange)) {
    var_CEP_00_II_up <- var_CEP_00_II_low
    var_CEP_10_II_up <- var_CEP_10_II_low
    var_CEP_diff_II_up <- var_CEP_diff_II_low
  } else {
    mats <-
      compute_matrices(list(eeFUN = up_fun,
                            splitdt = split(data_wide,
                                            f = data_wide$Group)),
                       numDeriv_options = list(method = 'simple'),
                       theta = thetahat_up, beta0range = brange,
                       contrast = contrast)

    A <- apply(simplify2array(Map(`*`, mats$A_i,
                                  data_wide$freq)), 1:2, sum)
    B <- apply(simplify2array(Map(`*`, mats$B_i,
                                  data_wide$freq)), 1:2, sum)
    Sigma <- solve(A) %*% B %*% t(solve(A))

    var_CEP_00_II_up <- Sigma[10, 10]
    var_CEP_10_II_up <- Sigma[12, 12]
    var_CEP_diff_II_up <- Sigma[14, 14]
  } # end of 'else'

  # compute Imbens-Manski intervals for each quantity
  maxsig_10 <- max(sqrt(var_CEP_10_II_low),
                   sqrt(var_CEP_10_II_up))
  cstar_10 <- uniroot(f_cstar, c(1.64, 1.96), low = CEP_10_II_low,
                            up = CEP_10_II_up, maxsig = maxsig_10)$root
  CEP_10_EUI_low <- CEP_10_II_low - cstar_10*sqrt(var_CEP_10_II_low)
  CEP_10_EUI_up <- CEP_10_II_up + cstar_10*sqrt(var_CEP_10_II_up)

  maxsig_00 <- max(sqrt(var_CEP_00_II_low),
                   sqrt(var_CEP_00_II_up))
  cstar_00 <- uniroot(f_cstar, c(1.64, 1.96), low = CEP_00_II_low,
                            up = CEP_00_II_up, maxsig = maxsig_00)$root
  CEP_00_EUI_low <- CEP_00_II_low - cstar_00*sqrt(var_CEP_00_II_low)
  CEP_00_EUI_up <- CEP_00_II_up + cstar_00*sqrt(var_CEP_00_II_up)

  maxsig_diff <- max(sqrt(var_CEP_diff_II_low),
                     sqrt(var_CEP_diff_II_up))
  cstar_diff <- uniroot(f_cstar, c(1.64,1.96), low = CEP_diff_II_low,
                              up = CEP_diff_II_up,
                              maxsig = maxsig_diff)$root
  CEP_diff_EUI_low <- CEP_diff_II_low-cstar_diff*sqrt(var_CEP_diff_II_low)
  CEP_diff_EUI_up <- CEP_diff_II_up + cstar_diff*sqrt(var_CEP_diff_II_up)

  return(list("CEP_10_II" = c(CEP_10_II_low, CEP_10_II_up),
              "CEP_00_II" = c(CEP_00_II_low, CEP_00_II_up),
              "CEP_diff_II" = c(CEP_diff_II_low, CEP_diff_II_up),
              "CEP_10_EUI" = c(CEP_10_EUI_low, CEP_10_EUI_up),
              "CEP_00_EUI" = c(CEP_00_EUI_low, CEP_00_EUI_up),
              "CEP_diff_EUI" = c(CEP_diff_EUI_low, CEP_diff_EUI_up)))

}
