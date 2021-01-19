#' Analyze under NEB assumptions
#'
#' Analyze data assuming no early treatment benefits
#'
#' @param data Data frame containing the following variables
#' \itemize{
#'    \item Z: indicator of treatment
#'    \item Y: indicator of outcome
#'    \item Y_tau: indicator of early outcome
#'    \item S_star: intermediate biomarker value
#'    \item R: indicator of measurement of intermediate biomarker
#' }
#' @param brange0 Numeric (2 x 1) vector containing the specified lower and upper bounds of the range for sensitivity parameter \ifelse{html}{\out{&#946;<sub>0</sub>}}{\eqn{\beta_0}}
#' @param brange1 Numeric (2 x 1) vector containing the specified lower and upper bounds of the range for sensitivity parameter \ifelse{html}{\out{&#946;<sub>1</sub>}}{\eqn{\beta_1}}
#' @param design String describing the study design / sampling scheme used. This allows for estimation of sampling weights. Options include "full", "cc" (case-cohort), and "other". When "other" is chosen the weights argument must also be specified
#' @param weights Numeric (n x 1) vector containing pre-estimated sampling weights where n is the number of rows in `data`
#' @param contrast Contrast function for estimand. Options include "logRR", "Difference", and "VE"
#'
#' @return Returns list consisting of 6 vectors corresponding to the ignorance intervals and EUIs of CEP(1, 0), CEP(0, 0), and the difference CEP(1, 0) - CEP(0, 0)
#' @importFrom stats uniroot
#' @importFrom stats pnorm
#' @export
#'
#' @examples Z <- rbinom(500, 1, 0.5)
#' S_star <- rbinom(500, 1, 0.2)
#' R <- rep(1, 500)
#' Y_tau_0 <- rbinom(500, 1, 0.02)
#' Y_tau_1 <- Y_tau_0 + rbinom(500, 1, (1-Y_tau_0)*Z*0.02)
#' Y_tau <- Y_tau_0*(1-Z) + Y_tau_1*Z
#' Y <- rbinom(500, 1, 0.1)
#' df <- data.frame(Z, S_star, R, Y_tau, Y)
#' analyze_NEB(df, c(-0.5, 0.5), design = "full", contrast = "VE")
analyze_NEB <- function(data, brange0 = c(0, 0), brange1 = c(0, 0),
                        design = "full", weights = NULL,
                        contrast = "logRR") {

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
  # Estimate probability control has S_star observed: P(R = 1 | Y = 0)
  solve_pi <- function(x) { sum((1 - Y)*(R - x)) }
  pi_hat <- uniroot(solve_pi, c(0, 1))$root
  if (tolower(design) == "cc") {
    W <- 1 / pi_hat*(1 - Y)*R + Y
  }
  if (!is.null(weights)) {
    W <- weights
  }


  # Point estimation

  S_0 <- (1 - Y_tau)*(1 - S_star)
  S_1 <- (1 - Y_tau)*S_star

  # p(1, 0) = P(S^tau(1) = 1, S^tau(0) = 0 | Y^tau(0) = Y^tau(1) = 0)
  solve_p10 <- function(x) {
    sum((1 - Y_tau)*Z*(S_star - x)*W)
  }
  p_10 <- uniroot(solve_p10, c(0, 1))$root
  p_00 <- 1 - p_10

  # Estimate risk_1(0, 0), risk_1(1, 0) where risk_z(s1, s0) =
  # P(Y(z) = 1 | S^tau(1) = s1, S^tau(0) = s0, Y^tau(1) = Y^tau(0) = 0)
  solve_risk_1_00 <- function(x) {
    sum((1 - Y_tau)*Z*(!is.na(S_star) & S_star == 0)*(Y - x)*W)
  }
  risk_1_00 <- uniroot(solve_risk_1_00, c(0, 1))$root

  solve_risk_1_10 <- function(x) {
    sum((1 - Y_tau)*Z*(!is.na(S_star) & S_star == 1)*(Y - x)*W)
  }
  risk_1_10 <- uniroot(solve_risk_1_10, c(0, 1))$root

  # Estimate partially identifiable terms using SACE methods
  # Estimate risk_0 = P(Y(0) = 1 | Y^tau(0) = Y^tau(1) = 0)
  mix_1 <- 1 - mean(Y_tau[Z == 1])
  mix_2 <- 1 - mean(Y_tau[Z == 0])
  mix <- mix_1 / mix_2
  prob <- mean(Y[Y_tau == 0 & Z == 0])

  solve_risk0_min <- function(x) {
    risk_new <- 1 /( 1 + exp(min(brange0))*(1 - x) / x )
    prob - x*mix - risk_new*(1 - mix)
  }
  risk_0_first <- uniroot(solve_risk0_min, c(0, 1))$root
  risk_new_first <- (prob - risk_0_first*mix) / (1 - mix)

  solve_risk0_max <- function(x) {
    risk_new <- 1 /( 1 + exp(max(brange0))*(1 - x) / x )
    prob - x*mix - risk_new*(1 - mix)
  }
  risk_0_second <- uniroot(solve_risk0_max, c(0, 1))$root
  risk_new_second <- (prob - risk_0_second*mix) / (1 - mix)

  # Estimate risk_0_00 and risk_0_10
  solve_risk_0_00_minmin <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(min(brange1))*(1 - x) / x )
    risk_0_first - x*p_00 - risk_0_10*p_10
  }
  risk_0_00_first <- uniroot(solve_risk_0_00_minmin, c(0, 1))$root
  risk_0_10_first <- (risk_0_first - risk_0_00_first*p_00) / p_10

  solve_risk_0_00_maxmin <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(max(brange1))*(1 - x) / x )
    risk_0_first - x*p_00 - risk_0_10*p_10
  }
  risk_0_00_second <- uniroot(solve_risk_0_00_maxmin, c(0, 1))$root
  risk_0_10_second <- (risk_0_first - risk_0_00_second*p_00) / p_10

  solve_risk_0_00_minmax <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(min(brange1))*(1-x) / x )
    risk_0_second - x*p_00 - risk_0_10*p_10
  }
  risk_0_00_third <- uniroot(solve_risk_0_00_minmax, c(0, 1))$root
  risk_0_10_third <- (risk_0_second - risk_0_00_third*p_00) / p_10

  solve_risk_0_00_maxmax <- function(x) {
    risk_0_10 <- 1 / ( 1 + exp(max(brange1))*(1-x) / x )
    risk_0_second - x*p_00 - risk_0_10*p_10
  }
  risk_0_00_fourth <- uniroot(solve_risk_0_00_maxmax, c(0, 1))$root
  risk_0_10_fourth <- (risk_0_second - risk_0_00_fourth*p_00) / p_10

  # Calculate CEP(1,0) ignorance interval
  CEP_10_II_low <- min(h(risk_1_10, risk_0_10_first, contrast),
                       h(risk_1_10, risk_0_10_second, contrast),
                       h(risk_1_10, risk_0_10_third, contrast),
                       h(risk_1_10, risk_0_10_fourth, contrast))
  CEP_10_II_up <- max(h(risk_1_10, risk_0_10_first, contrast),
                      h(risk_1_10, risk_0_10_second, contrast),
                      h(risk_1_10, risk_0_10_third, contrast),
                      h(risk_1_10, risk_0_10_fourth, contrast))
  whichmin_10 <- which.min(c(h(risk_1_10, risk_0_10_first, contrast),
                             h(risk_1_10, risk_0_10_second, contrast),
                             h(risk_1_10, risk_0_10_third, contrast),
                             h(risk_1_10, risk_0_10_fourth, contrast)))
  whichmax_10 <- which.max(c(h(risk_1_10, risk_0_10_first, contrast),
                             h(risk_1_10, risk_0_10_second, contrast),
                             h(risk_1_10, risk_0_10_third, contrast),
                             h(risk_1_10, risk_0_10_fourth, contrast)))

  # Calculate CEP(0,0) ignorance interval
  CEP_00_II_low <- min(h(risk_1_00, risk_0_00_first, contrast),
                       h(risk_1_00, risk_0_00_second, contrast),
                       h(risk_1_00, risk_0_00_third, contrast),
                       h(risk_1_00, risk_0_00_fourth, contrast))
  CEP_00_II_up <- max(h(risk_1_00, risk_0_00_first, contrast),
                      h(risk_1_00, risk_0_00_second, contrast),
                      h(risk_1_00, risk_0_00_third, contrast),
                      h(risk_1_00, risk_0_00_fourth, contrast))
  whichmin_00 <- which.min(c(h(risk_1_00, risk_0_00_first, contrast),
                             h(risk_1_00, risk_0_00_second, contrast),
                             h(risk_1_00, risk_0_00_third, contrast),
                             h(risk_1_00, risk_0_00_fourth, contrast)))
  whichmax_00 <- which.max(c(h(risk_1_00, risk_0_00_first, contrast),
                             h(risk_1_00, risk_0_00_second, contrast),
                             h(risk_1_00, risk_0_00_third, contrast),
                             h(risk_1_00, risk_0_00_fourth, contrast)))

  # Calculate CEP(1,0) - CEP(0,0) ignorance interval
  CEP_diff_II_low <- min(h(risk_1_10, risk_0_10_first, contrast) -
                           h(risk_1_00, risk_0_00_first, contrast),
                         h(risk_1_10, risk_0_10_second, contrast) -
                           h(risk_1_00, risk_0_00_second, contrast),
                         h(risk_1_10, risk_0_10_third, contrast) -
                           h(risk_1_00, risk_0_00_third, contrast),
                         h(risk_1_10, risk_0_10_fourth, contrast) -
                           h(risk_1_00, risk_0_00_fourth, contrast))
  CEP_diff_II_up <- max(h(risk_1_10, risk_0_10_first, contrast) -
                          h(risk_1_00, risk_0_00_first, contrast),
                        h(risk_1_10, risk_0_10_second, contrast) -
                          h(risk_1_00, risk_0_00_second, contrast),
                        h(risk_1_10, risk_0_10_third, contrast) -
                          h(risk_1_00, risk_0_00_third, contrast),
                        h(risk_1_10, risk_0_10_fourth, contrast) -
                          h(risk_1_00, risk_0_00_fourth, contrast))
  whichmin_diff <- which.min(c(h(risk_1_10, risk_0_10_first, contrast) -
                                 h(risk_1_00, risk_0_00_first, contrast),
                               h(risk_1_10, risk_0_10_second, contrast) -
                                 h(risk_1_00, risk_0_00_second, contrast),
                               h(risk_1_10, risk_0_10_third, contrast) -
                                 h(risk_1_00, risk_0_00_third, contrast),
                               h(risk_1_10, risk_0_10_fourth, contrast) -
                                 h(risk_1_00, risk_0_00_fourth, contrast)))
  whichmax_diff <- which.max(c(h(risk_1_10, risk_0_10_first, contrast) -
                                 h(risk_1_00, risk_0_00_first, contrast),
                               h(risk_1_10, risk_0_10_second, contrast) -
                                 h(risk_1_00, risk_0_00_second, contrast),
                               h(risk_1_10, risk_0_10_third, contrast) -
                                 h(risk_1_00, risk_0_00_third, contrast),
                               h(risk_1_10, risk_0_10_fourth, contrast) -
                                 h(risk_1_00, risk_0_00_fourth, contrast)))



  # Variance estimation

  # Get point estimate vector needed for sandwich variance estimation
  thetahat <- as.numeric(c(risk_1_00, risk_1_10, mix_1, mix_2, mix, prob,
                           p_10, risk_0_first, risk_new_first,
                           risk_0_second, risk_new_second,
                           risk_0_00_first, risk_0_10_first,
                           risk_0_00_second, risk_0_10_second,
                           risk_0_00_third, risk_0_10_third,
                           risk_0_00_fourth, risk_0_10_fourth,
                           CEP_00_II_low, CEP_00_II_up, CEP_10_II_low,
                           CEP_10_II_up, CEP_diff_II_low, CEP_diff_II_up,
                           pi_hat))

  # Convert data to wide form
  data$W <- W
  data_wide <- plyr::count(data, c('Y', 'Z', 'Y_tau', 'S_star', 'R', 'W'))
  data_wide$Group <- 1:nrow(data_wide)

  # Calculate variance estimate for the ignorance intervals using
  # 'geex' package source code
  mats <- compute_matrices(list(eeFUN = eefun_NEB,
                                splitdt = split(data_wide,
                                                f = data_wide$Group)),
                           numDeriv_options = list(method = 'simple'),
                           theta = thetahat,
                           beta0range = brange0, beta1range = brange1,
                           contrast = contrast, design = design,
                           whichmin_00 = whichmin_00,
                           whichmax_00 = whichmax_00,
                           whichmin_10 = whichmin_10,
                           whichmax_10 = whichmax_10,
                           whichmin_diff = whichmin_diff,
                           whichmax_diff = whichmax_diff)

  A <- apply(simplify2array(Map(`*`, mats$A_i,
                                data_wide$freq)), 1:2, sum)
  B <- apply(simplify2array(Map(`*`, mats$B_i,
                                data_wide$freq)), 1:2, sum)

  Sigma <- solve(A) %*% B %*% t(solve(A))

  # Extract variances of interest from sandwich matrix
  var_CEP_00_II_low <- Sigma[20, 20]
  var_CEP_00_II_up  <- Sigma[21, 21]
  var_CEP_10_II_low <- Sigma[22, 22]
  var_CEP_10_II_up  <- Sigma[23, 23]
  var_CEP_diff_II_low <- Sigma[24, 24]
  var_CEP_diff_II_up  <- Sigma[25, 25]

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
