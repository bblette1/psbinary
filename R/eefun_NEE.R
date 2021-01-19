# Estimating equation vector for ignorance interval lower bound under NEE
# Used for M estimation variance calculation
eefun_NEE <- function(data, beta0range, contrast, design, whichmin_00,
                      whichmax_00, whichmin_10, whichmax_10,
                      whichmin_diff, whichmax_diff) {

  function(theta) {
    with(data,
         c( # Identifiable terms

           # risk_0 = theta[1]
           (1-Y_tau)*(1-Z)*(Y-theta[1]),
           # risk_1_10 = theta[2]
           (1-Y_tau)*W*S_star*Z*(Y-theta[2])*(design != "cc") +
             (1-Y_tau)*W*S_star*Z*(1/theta[15]*(1-Y)*R+Y)*
             (Y-theta[2])*(design == "cc"),
           # p_00 = theta[3]
           (1-Y_tau)*W*Z*(1-S_star-theta[3])*(design != "cc") +
             (1-Y_tau)*W*Z*(1/theta[15]*(1-Y)*R+Y)*
             (1-S_star-theta[3])*(design == "cc"),
           # risk_1_00 = theta[4]
           Z*(1-Y_tau)*W*(1-S_star)*(Y-theta[4])*(design != "cc") +
             Z*(1-Y_tau)*W*(1-S_star)*(1/theta[15]*(1-Y)*R+Y)*
             (Y-theta[4])*(design == "cc"),

           # Partially identifiable terms

           # risk_0_00_first = theta[5] for minimum of beta0
           theta[1] - theta[5]*theta[3] - theta[6]*(1-theta[3]),
           # risk_0_10_first = theta[6] for minimum of beta0
           exp(min(beta0range))*theta[6]/(1-theta[6]) -
             theta[5]/(1-theta[5]),
           # risk_0_00_second = theta[7] for maximum of beta0
           theta[1] - theta[7]*theta[3] - theta[8]*(1-theta[3]),
           # risk_0_10_second = theta[8] for maximum of beta0
           exp(max(beta0range))*theta[8]/(1-theta[8]) -
             theta[7]/(1-theta[7]),

           # CEP terms

           # Lower bound of CEP(0, 0) = theta[9]
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

           # Upper bound of CEP(0, 0) = theta[10]
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

           # Lower bound of CEP(1, 0) = theta[11]
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

           # Upper bound of CEP(1, 0) = theta[12]
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

           # Lower bound of CEP(1, 0) - CEP(0, 0) = theta[13]
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

           # Upper bound of CEP(1, 0) - CEP(0, 0) = theta[14]
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

           # Case-cohort weights, pi = theta[15]
           (1-Y)*(R-theta[15])
         ))
  }

}
