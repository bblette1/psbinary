# Estimating equation vector for ignorance interval lower bound under NEB
# Used for M estimation variance calculation
eefun_NEB <- function(data, beta0range, beta1range, contrast, design,
                      whichmin_00, whichmax_00, whichmin_10, whichmax_10,
                      whichmin_diff, whichmax_diff) {

  function(theta) {
    with(data,
         c( # Identifiable terms

           # risk_1_00 = theta[1]
           (1-Y_tau)*Z*1*(!is.na(S_star) & S_star==0)*(Y-theta[1]),
           # risk_1_10 = theta[2]
           (1-Y_tau)*Z*1*(!is.na(S_star) & S_star==1)*(Y-theta[2]),
           # P(Ytau(1) = 0) = theta[3]
           Z*(1-Y_tau-theta[3]),
           # P(Ytau(0) = 0) = theta[4]
           (1-Z)*(1-Yt_au-theta[4]),
           # P(Ytau(1) = 0 | Ytau(0) = 0) = theta[5]
           theta[5] - theta[3]/theta[4],
           # P(Y(0) = 1 | Ytau(0) = 0) = theta[6]
           (1-Z)*(1-Y_tau)*(Y-theta[6]),
           # p_10 = theta[7]
           (1-Ytau)*Z*(1*(!is.na(S_star) & S_star==1)-theta[9]),

           # Partially identifiable terms

           # risk_0_first = theta[8] for minimum of beta0
           theta[6] - theta[8]*theta[5] - theta[9]*(1-theta[5]),
           # risk_new_first = theta[9] for minimum of beta0
           exp(min(beta0range))*theta[9]/(1-theta[9]) -
             theta[8]/(1-theta[8]),
           # risk_0_second = theta[10] for maximum of beta0
           theta[6] - theta[10]*theta[5] - theta[11]*(1-theta[5]),
           # risk_new_second = theta[11] for maximum of beta0
           exp(max(beta0range))*theta[11]/(1-theta[11]) -
             theta[10]/(1-theta[10]),

           # risk_0_00_first = theta[12] for min/min of beta0/beta1
           theta[8] - theta[12]*(1-theta[7]) - theta[13]*theta[7],
           # risk_0_10_first = theta[13] for min/min of beta0/beta1
           exp(min(beta1range))*theta[13]/(1-theta[13]) -
             theta[12]/(1-theta[12]),
           # risk_0_00_second = theta[14] for min/max of beta0/beta1
           theta[8] - theta[14]*(1-theta[7]) - theta[15]*theta[7],
           # risk_0_10_second = theta[15] for min/max of beta0/beta1
           exp(max(beta1range))*theta[15]/(1-theta[15]) -
             theta[14]/(1-theta[14]),
           # risk_0_00_third = theta[16] for max/min of beta0/beta1
           theta[10] - theta[16]*(1-theta[7]) - theta[17]*theta[7],
           # risk_0_10_third = theta[17] for max/min of beta0/beta1
           exp(min(beta1range))*theta[17]/(1-theta[17]) -
             theta[16]/(1-theta[16]),
           # risk_0_00_fourth = theta[18] for max/max of beta0/beta1
           theta[10] - theta[18]*(1-theta[7]) - theta[19]*theta[7],
           # risk_0_10_fourth = theta[19] for max/max of beta0/beta1
           exp(max(beta1range))*theta[19]/(1-theta[19]) -
             theta[18]/(1-theta[18]),

           # CEP terms

           # Lower bound of CEP(0, 0) = theta[20]
           ( theta[9] - (1- theta[4]/theta[5]) )*
             (contrast == "VE" & whichmin_00 == 1) +
             ( theta[9] - (1- theta[4]/theta[7]) )*
             (contrast == "VE" & whichmin_00 == 2) +
             ( theta[9] - (theta[4] - theta[5]) )*
             (contrast == "Difference" & whichmin_00 == 1) +
             ( theta[9] - (theta[4] - theta[7]) )*
             (contrast == "Difference" & whichmin_00 == 2) +
             ( theta[9] - (log(theta[4]) - log(theta[5])) )*
             (contrast == "logRR" & whichmin_00 == 1) +
             ( theta[9] - (log(theta[4]) - log(theta[7])) )*
             (contrast == "logRR" & whichmin_00 == 2),

           # Upper bound of CEP(0, 0) = theta[21]
           ( theta[10] - (1- theta[4]/theta[5]) )*
             (contrast == "VE" & whichmax_00 == 1) +
             ( theta[10] - (1- theta[4]/theta[7]) )*
             (contrast == "VE" & whichmax_00 == 2) +
             ( theta[10] - (theta[4] - theta[5]) )*
             (contrast == "Difference" & whichmax_00 == 1) +
             ( theta[10] - (theta[4] - theta[7]) )*
             (contrast == "Difference" & whichmax_00 == 2) +
             ( theta[10] - (log(theta[4]) - log(theta[5])) )*
             (contrast == "logRR" & whichmax_00 == 1) +
             ( theta[10] - (log(theta[4]) - log(theta[7])) )*
             (contrast == "logRR" & whichmax_00 == 2),

           # Lower bound of CEP(1, 0) = theta[22]
           ( theta[11] - (1- theta[2]/theta[6]) )*
             (contrast == "VE" & whichmin_10 == 1) +
             ( theta[11] - (1- theta[2]/theta[8]) )*
             (contrast == "VE" & whichmin_10 == 2) +
             ( theta[11] - (theta[2] - theta[6]) )*
             (contrast == "Difference" & whichmin_10 == 1) +
             ( theta[11] - (theta[2] - theta[8]) )*
             (contrast == "Difference" & whichmin_10 == 2) +
             ( theta[11] - (log(theta[2]) - log(theta[6])) )*
             (contrast == "logRR" & whichmin_10 == 1) +
             ( theta[11] - (log(theta[2]) - log(theta[8])) )*
             (contrast == "logRR" & whichmin_10 == 2),

           # Upper bound of CEP(1, 0) = theta[23]
           ( theta[12] - (1- theta[2]/theta[6]) )*
             (contrast == "VE" & whichmax_10 == 1) +
             ( theta[12] - (1- theta[2]/theta[8]) )*
             (contrast == "VE" & whichmax_10 == 2) +
             ( theta[12] - (theta[2] - theta[6]) )*
             (contrast == "Difference" & whichmax_10 == 1) +
             ( theta[12] - (theta[2] - theta[8]) )*
             (contrast == "Difference" & whichmax_10 == 2) +
             ( theta[12] - (log(theta[2]) - log(theta[6])) )*
             (contrast == "logRR" & whichmax_10 == 1) +
             ( theta[12] - (log(theta[2]) - log(theta[8])) )*
             (contrast == "logRR" & whichmax_10 == 2),

           # Lower bound of CEP(1, 0) - CEP(0, 0) = theta[24]
           ( theta[13] - (1-theta[2]/theta[6]) +(1-theta[4]/theta[5]) )*
             (contrast == "VE" & whichmin_diff == 1) +
             ( theta[13] - (1-theta[2]/theta[8]) +(1-theta[4]/theta[7]) )*
             (contrast == "VE" & whichmin_diff == 2) +
             ( theta[13] - (theta[2]-theta[6]) + (theta[4]-theta[5]) )*
             (contrast == "Difference" & whichmin_diff == 1) +
             ( theta[13] - (theta[2]-theta[8]) + (theta[4]-theta[7]) )*
             (contrast == "Difference" & whichmin_diff == 2) +
             ( theta[13] - (log(theta[2]) - log(theta[6])) +
                 (log(theta[4]) - log(theta[5])) )*
             (contrast == "logRR" & whichmin_diff == 1) +
             ( theta[13] - (log(theta[2]) - log(theta[8])) +
                 (log(theta[4]) - log(theta[7])) )*
             (contrast == "logRR" & whichmin_diff == 2),

           # Upper bound of CEP(1, 0) - CEP(0, 0) = theta[25]
           ( theta[14] - (1-theta[2]/theta[6]) +(1-theta[4]/theta[5]) )*
             (contrast == "VE" & whichmax_diff == 1) +
             ( theta[14] - (1-theta[2]/theta[8]) +(1-theta[4]/theta[7]) )*
             (contrast == "VE" & whichmax_diff == 2) +
             ( theta[14] - (theta[2]-theta[6]) + (theta[4]-theta[5]) )*
             (contrast == "Difference" & whichmax_diff == 1) +
             ( theta[14] - (theta[2]-theta[8]) + (theta[4]-theta[7]) )*
             (contrast == "Difference" & whichmax_diff == 2) +
             ( theta[14] - (log(theta[2]) - log(theta[6])) +
                 (log(theta[4]) - log(theta[5])) )*
             (contrast == "logRR" & whichmax_diff == 1) +
             ( theta[14] - (log(theta[2]) - log(theta[8])) +
                 (log(theta[4]) - log(theta[7])) )*
             (contrast == "logRR" & whichmax_diff == 2),

           # Case-cohort weights, pi = theta[26]
           (1-Y)*(R-theta[15])
         ))
  }

}
