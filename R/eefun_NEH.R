# Estimating equation vector for ignorance interval under NEH
# Used for M estimation variance calculation
eefun_NEH <- function(data, beta0range, beta1range, beta2range, beta3range,
                      contrast, design, whichmin_00, whichmax_00,
                      whichmin_10, whichmax_10,
                      whichmin_diff, whichmax_diff) {

  function(theta) {
    with(data,
         c( # Identifiable terms

           # risk_0 = theta[1]
           (1-Y_tau)*(1-Z)*(Y-theta[1]),
           # v = theta[2]
           theta[2]*(1-Y_tau)*Z - (1-Y_tau)*(1-Z),
           # p_1 = theta[3]
           (1-Y_tau)*Z*((!is.na(S_star) & S_star==1)-theta[3])*W*
             (design != "cc") +
             (1-Y_tau)*Z*((!is.na(S_star) & S_star==1)-theta[3])*W*
             (1/theta[47]*(1-Y)*R+Y)*(design == "cc"),

           # Partially identifiable terms and helpers

           # p_10_first = theta[4]
           theta[4]*theta[2] + theta[5]*(1-theta[2]) - theta[3],
           # p_10_mix_first = theta[5]
           exp(min(beta3range))*theta[5]/(1-theta[5]) -
             theta[4]/(1-theta[4]),
           # p_10_second = theta[6]
           theta[6]*theta[2] + theta[7]*(1-theta[2]) - theta[3],
           # p_10_mix_second = theta[7]
           exp(max(beta3range))*theta[7]/(1-theta[7]) -
             theta[6]/(1-theta[6]),

           # risk_0_00_first = theta[8]
           theta[1] - theta[8]*(1-theta[4]) - theta[9]*theta[4],
           # risk_0_10_first = theta[9]
           exp(min(beta0range))*theta[9]/(1-theta[9]) -
             theta[8]/(1-theta[8]),
           # risk_0_00_second = theta[10]
           theta[1] - theta[10]*(1-theta[4]) - theta[11]*theta[4],
           # risk_0_10_second = theta[11]
           exp(max(beta0range))*theta[11]/(1-theta[11]) -
             theta[10]/(1-theta[10]),
           # risk_0_00_third = theta[12]
           theta[1] - theta[12]*(1-theta[6]) - theta[13]*theta[6],
           # risk_0_10_third = theta[13]
           exp(min(beta0range))*theta[13]/(1-theta[13]) -
             theta[12]/(1-theta[12]),
           # risk_0_00_fourth = theta[14]
           theta[1] - theta[14]*(1-theta[6]) - theta[15]*theta[6],
           # risk_0_10_fourth = theta[15]
           exp(max(beta0range))*theta[15]/(1-theta[15]) -
             theta[14]/(1-theta[14]),

           # row1 = theta[16]
           Z*(1-Y_tau)*(!is.na(S_star) & S_star==0)*W*(Y-theta[16])*
             (design != "cc") +
             Z*(1-Y_tau)*(!is.na(S_star) & S_star==0)*W*(Y-theta[16])*
             (1/theta[47]*(1-Y)*R+Y)*(design == "cc"),
           # Ytau_prob = theta[17]
           (1-Z)*(1-Y_tau-theta[17]),
           # r1_00_helper = theta[18]
           (1-Y_tau)*(!is.na(S_star) & S_star==0)*Z*W - Z*theta[18]*
             (design != "cc") +
             (1-Y_tau)*(!is.na(S_star) & S_star==0)*Z*W - Z*theta[18]*
             (1/theta[47]*(1-Y)*R+Y)*(design == "cc"),
           # r1_00_mix_first = theta[19]
           (1-theta[4])*theta[17] - theta[19]*theta[18],										           # r1_00_mix_second = theta[20]
           (1-theta[6])*theta[17] - theta[20]*theta[18],

           # risk_1_00_first = theta[21]
           theta[16] - theta[21]*theta[19] - theta[22]*(1-theta[19]),
           # risk_1_0star_first = theta[22]
           exp(min(beta1range))*theta[22]/(1-theta[22]) -
             theta[21]/(1-theta[21]),
           # risk_1_00_second = theta[23]
           theta[16] - theta[23]*theta[19] - theta[24]*(1-theta[19]),
           # risk_1_0star_second = theta[24]
           exp(max(beta1range))*theta[24]/(1-theta[24]) -
             theta[23]/(1-theta[23]),
           # risk_1_00_third = theta[25]
           theta[16] - theta[25]*theta[20] - theta[26]*(1-theta[20]),
           # risk_1_0star_third = theta[26]
           exp(min(beta1range))*theta[26]/(1-theta[26]) -
             theta[25]/(1-theta[25]),
           # risk_1_00_fourth = theta[27]
           theta[16] - theta[27]*theta[20] - theta[28]*(1-theta[20]),
           # risk_1_0star_fourth = theta[28]
           exp(max(beta1range))*theta[28]/(1-theta[28]) -
             theta[27]/(1-theta[27]),

           # row2 = theta[29]
           Z*(1-Y_tau)*(!is.na(S_star) & S_star==1)*W*(Y-theta[29])*
             (design != "cc") +
             Z*(1-Y_tau)*(!is.na(S_star) & S_star==1)*W*(Y-theta[29])*
             (1/theta[47]*(1-Y)*R+Y)*(design == "cc"),
           # Ytau_prob1 = theta[30]
           Z*(1-Y_tau-theta[30]),
           # r_10_mix_first = theta[31]
           theta[4]*theta[19] - theta[31]*theta[30]*theta[3],
           # r_10_mix_second = theta[32]
           theta[6]*theta[20] - theta[32]*theta[30]*theta[3],

           # risk_1_10_first = theta[33]
           theta[29] - theta[33]*theta[31] - theta[34]*(1-theta[31]),
           # risk_1_1star_first = theta[34]
           exp(min(beta2range))*theta[34]/(1-theta[34]) -
             theta[33]/(1-theta[33]),
           # risk_1_10_second = theta[35]
           theta[29] - theta[35]*theta[31] - theta[36]*(1-theta[31]),
           # risk_1_1star_second = theta[36]
           exp(max(beta2range))*theta[36]/(1-theta[36]) -
             theta[35]/(1-theta[35]),
           # risk_1_10_third = theta[37]
           theta[29] - theta[37]*theta[32] - theta[38]*(1-theta[32]),
           # risk_1_1star_third = theta[38]
           exp(min(beta2range))*theta[38]/(1-theta[38]) -
             theta[37]/(1-theta[37]),
           # risk_1_10_fourth = theta[39]
           theta[29] - theta[39]*theta[32] - theta[40]*(1-theta[32]),
           # risk_1_1star_fourth = theta[40]
           exp(max(beta2range))*theta[40]/(1-theta[40]) -
             theta[39]/(1-theta[39]),

           # CEP terms

           # Lower bound of CEP(0, 0) = theta[41]
           ( theta[41] - (1- theta[21]/theta[8]) )*
             (contrast == "VE" & whichmin_00 == 1) +
             ( theta[41] - (1- theta[21]/theta[10]) )*
             (contrast == "VE" & whichmin_00 == 2) +
             ( theta[41] - (1- theta[23]/theta[8]) )*
             (contrast == "VE" & whichmin_00 == 3) +
             ( theta[41] - (1- theta[23]/theta[10]) )*
             (contrast == "VE" & whichmin_00 == 4) +
             ( theta[41] - (1- theta[25]/theta[12]) )*
             (contrast == "VE" & whichmin_00 == 5) +
             ( theta[41] - (1- theta[25]/theta[14]) )*
             (contrast == "VE" & whichmin_00 == 6) +
             ( theta[41] - (1- theta[27]/theta[12]) )*
             (contrast == "VE" & whichmin_00 == 7) +
             ( theta[41] - (1- theta[27]/theta[14]) )*
             (contrast == "VE" & whichmin_00 == 8) +

             ( theta[41] - (theta[21] - theta[8]) )*
             (contrast == "Difference" & whichmin_00 == 1) +
             ( theta[41] - (theta[21] - theta[10]) )*
             (contrast == "Difference" & whichmin_00 == 2) +
             ( theta[41] - (theta[23] - theta[8]) )*
             (contrast == "Difference" & whichmin_00 == 3) +
             ( theta[41] - (theta[23] - theta[10]) )*
             (contrast == "Difference" & whichmin_00 == 4) +
             ( theta[41] - (theta[25] - theta[12]) )*
             (contrast == "Difference" & whichmin_00 == 5) +
             ( theta[41] - (theta[25] - theta[14]) )*
             (contrast == "Difference" & whichmin_00 == 6) +
             ( theta[41] - (theta[27] - theta[12]) )*
             (contrast == "Difference" & whichmin_00 == 7) +
             ( theta[41] - (theta[27] - theta[14]) )*
             (contrast == "Difference" & whichmin_00 == 8) +

             ( theta[41] - (log(theta[21]) - log(theta[8])) )*
             (contrast == "logRR" & whichmin_00 == 1) +
             ( theta[41] - (log(theta[21]) - log(theta[10])) )*
             (contrast == "logRR" & whichmin_00 == 2) +
             ( theta[41] - (log(theta[23]) - log(theta[8])) )*
             (contrast == "logRR" & whichmin_00 == 3) +
             ( theta[41] - (log(theta[23]) - log(theta[10])) )*
             (contrast == "logRR" & whichmin_00 == 4) +
             ( theta[41] - (log(theta[25]) - log(theta[12])) )*
             (contrast == "logRR" & whichmin_00 == 5) +
             ( theta[41] - (log(theta[25]) - log(theta[14])) )*
             (contrast == "logRR" & whichmin_00 == 6) +
             ( theta[41] - (log(theta[27]) - log(theta[12])) )*
             (contrast == "logRR" & whichmin_00 == 7) +
             ( theta[41] - (log(theta[27]) - log(theta[14])) )*
             (contrast == "logRR" & whichmin_00 == 8),

           # Upper bound of CEP(0, 0) = theta[42]
           ( theta[42] - (1- theta[21]/theta[8]) )*
             (contrast == "VE" & whichmax_00 == 1) +
             ( theta[42] - (1- theta[21]/theta[10]) )*
             (contrast == "VE" & whichmax_00 == 2) +
             ( theta[42] - (1- theta[23]/theta[8]) )*
             (contrast == "VE" & whichmax_00 == 3) +
             ( theta[42] - (1- theta[23]/theta[10]) )*
             (contrast == "VE" & whichmax_00 == 4) +
             ( theta[42] - (1- theta[25]/theta[12]) )*
             (contrast == "VE" & whichmax_00 == 5) +
             ( theta[42] - (1- theta[25]/theta[14]) )*
             (contrast == "VE" & whichmax_00 == 6) +
             ( theta[42] - (1- theta[27]/theta[12]) )*
             (contrast == "VE" & whichmax_00 == 7) +
             ( theta[42] - (1- theta[27]/theta[14]) )*
             (contrast == "VE" & whichmax_00 == 8) +

             ( theta[42] - (theta[21] - theta[8]) )*
             (contrast == "Difference" & whichmax_00 == 1) +
             ( theta[42] - (theta[21] - theta[10]) )*
             (contrast == "Difference" & whichmax_00 == 2) +
             ( theta[42] - (theta[23] - theta[8]) )*
             (contrast == "Difference" & whichmax_00 == 3) +
             ( theta[42] - (theta[23] - theta[10]) )*
             (contrast == "Difference" & whichmax_00 == 4) +
             ( theta[42] - (theta[25] - theta[12]) )*
             (contrast == "Difference" & whichmax_00 == 5) +
             ( theta[42] - (theta[25] - theta[14]) )*
             (contrast == "Difference" & whichmax_00 == 6) +
             ( theta[42] - (theta[27] - theta[12]) )*
             (contrast == "Difference" & whichmax_00 == 7) +
             ( theta[42] - (theta[27] - theta[14]) )*
             (contrast == "Difference" & whichmax_00 == 8) +

             ( theta[42] - (log(theta[21]) - log(theta[8])) )*
             (contrast == "logRR" & whichmax_00 == 1) +
             ( theta[42] - (log(theta[21]) - log(theta[10])) )*
             (contrast == "logRR" & whichmax_00 == 2) +
             ( theta[42] - (log(theta[23]) - log(theta[8])) )*
             (contrast == "logRR" & whichmax_00 == 3) +
             ( theta[42] - (log(theta[23]) - log(theta[10])) )*
             (contrast == "logRR" & whichmax_00 == 4) +
             ( theta[42] - (log(theta[25]) - log(theta[12])) )*
             (contrast == "logRR" & whichmax_00 == 5) +
             ( theta[42] - (log(theta[25]) - log(theta[14])) )*
             (contrast == "logRR" & whichmax_00 == 6) +
             ( theta[42] - (log(theta[27]) - log(theta[12])) )*
             (contrast == "logRR" & whichmax_00 == 7) +
             ( theta[42] - (log(theta[27]) - log(theta[14])) )*
             (contrast == "logRR" & whichmax_00 == 8),

           # Lower bound of CEP(1, 0) = theta[43]
           ( theta[43] - (1- theta[33]/theta[9]) )*
             (contrast == "VE" & whichmin_10 == 1) +
             ( theta[43] - (1- theta[33]/theta[11]) )*
             (contrast == "VE" & whichmin_10 == 2) +
             ( theta[43] - (1- theta[35]/theta[9]) )*
             (contrast == "VE" & whichmin_10 == 3) +
             ( theta[43] - (1- theta[35]/theta[11]) )*
             (contrast == "VE" & whichmin_10 == 4) +
             ( theta[43] - (1- theta[37]/theta[13]) )*
             (contrast == "VE" & whichmin_10 == 5) +
             ( theta[43] - (1- theta[37]/theta[15]) )*
             (contrast == "VE" & whichmin_10 == 6) +
             ( theta[43] - (1- theta[39]/theta[13]) )*
             (contrast == "VE" & whichmin_10 == 7) +
             ( theta[43] - (1- theta[39]/theta[15]) )*
             (contrast == "VE" & whichmin_10 == 8) +

             ( theta[43] - (theta[33] - theta[9]) )*
             (contrast == "Difference" & whichmin_10 == 1) +
             ( theta[43] - (theta[33] - theta[11]) )*
             (contrast == "Difference" & whichmin_10 == 2) +
             ( theta[43] - (theta[35] - theta[9]) )*
             (contrast == "Difference" & whichmin_10 == 3) +
             ( theta[43] - (theta[35] - theta[11]) )*
             (contrast == "Difference" & whichmin_10 == 4) +
             ( theta[43] - (theta[37] - theta[13]) )*
             (contrast == "Difference" & whichmin_10 == 5) +
             ( theta[43] - (theta[37] - theta[15]) )*
             (contrast == "Difference" & whichmin_10 == 6) +
             ( theta[43] - (theta[39] - theta[13]) )*
             (contrast == "Difference" & whichmin_10 == 7) +
             ( theta[43] - (theta[39] - theta[15]) )*
             (contrast == "Difference" & whichmin_10 == 8) +

             ( theta[43] - (log(theta[33]) - log(theta[9])) )*
             (contrast == "logRR" & whichmin_10 == 1) +
             ( theta[43] - (log(theta[33]) - log(theta[11])) )*
             (contrast == "logRR" & whichmin_10 == 2) +
             ( theta[43] - (log(theta[35]) - log(theta[9])) )*
             (contrast == "logRR" & whichmin_10 == 3) +
             ( theta[43] - (log(theta[35]) - log(theta[11])) )*
             (contrast == "logRR" & whichmin_10 == 4) +
             ( theta[43] - (log(theta[37]) - log(theta[13])) )*
             (contrast == "logRR" & whichmin_10 == 5) +
             ( theta[43] - (log(theta[37]) - log(theta[15])) )*
             (contrast == "logRR" & whichmin_10 == 6) +
             ( theta[43] - (log(theta[39]) - log(theta[13])) )*
             (contrast == "logRR" & whichmin_10 == 7) +
             ( theta[43] - (log(theta[39]) - log(theta[15])) )*
             (contrast == "logRR" & whichmin_10 == 8),

           # Upper bound of CEP(1, 0) = theta[44]
           ( theta[44] - (1- theta[33]/theta[9]) )*
             (contrast == "VE" & whichmax_10 == 1) +
             ( theta[44] - (1- theta[33]/theta[11]) )*
             (contrast == "VE" & whichmax_10 == 2) +
             ( theta[44] - (1- theta[35]/theta[9]) )*
             (contrast == "VE" & whichmax_10 == 3) +
             ( theta[44] - (1- theta[35]/theta[11]) )*
             (contrast == "VE" & whichmax_10 == 4) +
             ( theta[44] - (1- theta[37]/theta[13]) )*
             (contrast == "VE" & whichmax_10 == 5) +
             ( theta[44] - (1- theta[37]/theta[15]) )*
             (contrast == "VE" & whichmax_10 == 6) +
             ( theta[44] - (1- theta[39]/theta[13]) )*
             (contrast == "VE" & whichmax_10 == 7) +
             ( theta[44] - (1- theta[39]/theta[15]) )*
             (contrast == "VE" & whichmax_10 == 8) +

             ( theta[44] - (theta[33] - theta[9]) )*
             (contrast == "Difference" & whichmax_10 == 1) +
             ( theta[44] - (theta[33] - theta[11]) )*
             (contrast == "Difference" & whichmax_10 == 2) +
             ( theta[44] - (theta[35] - theta[9]) )*
             (contrast == "Difference" & whichmax_10 == 3) +
             ( theta[44] - (theta[35] - theta[11]) )*
             (contrast == "Difference" & whichmax_10 == 4) +
             ( theta[44] - (theta[37] - theta[13]) )*
             (contrast == "Difference" & whichmax_10 == 5) +
             ( theta[44] - (theta[37] - theta[15]) )*
             (contrast == "Difference" & whichmax_10 == 6) +
             ( theta[44] - (theta[39] - theta[13]) )*
             (contrast == "Difference" & whichmax_10 == 7) +
             ( theta[44] - (theta[39] - theta[15]) )*
             (contrast == "Difference" & whichmax_10 == 8) +

             ( theta[44] - (log(theta[33]) - log(theta[9])) )*
             (contrast == "logRR" & whichmax_10 == 1) +
             ( theta[44] - (log(theta[33]) - log(theta[11])) )*
             (contrast == "logRR" & whichmax_10 == 2) +
             ( theta[44] - (log(theta[35]) - log(theta[9])) )*
             (contrast == "logRR" & whichmax_10 == 3) +
             ( theta[44] - (log(theta[35]) - log(theta[11])) )*
             (contrast == "logRR" & whichmax_10 == 4) +
             ( theta[44] - (log(theta[37]) - log(theta[13])) )*
             (contrast == "logRR" & whichmax_10 == 5) +
             ( theta[44] - (log(theta[37]) - log(theta[15])) )*
             (contrast == "logRR" & whichmax_10 == 6) +
             ( theta[44] - (log(theta[39]) - log(theta[13])) )*
             (contrast == "logRR" & whichmax_10 == 7) +
             ( theta[44] - (log(theta[39]) - log(theta[15])) )*
             (contrast == "logRR" & whichmax_10 == 8),

           # Lower bound of CEP(1, 0) - CEP(0, 0) = theta[45]
           ( theta[45] - (1 - theta[33]/theta[9]) +
               (1 - theta[21]/theta[8]) )*
               (contrast == "VE" & whichmin_diff == 1) +
             ( theta[45] - (1 - theta[33]/theta[9]) +
               (1 - theta[23]/theta[8]) )*
               (contrast == "VE" & whichmin_diff == 2) +
             ( theta[45] - (1 - theta[33]/theta[11]) +
               (1 - theta[21]/theta[10]) )*
               (contrast == "VE" & whichmin_diff == 3) +
             ( theta[45] - (1 - theta[33]/theta[11]) +
               (1 - theta[23]/theta[10]) )*
               (contrast == "VE" & whichmin_diff == 4) +
             ( theta[45] - (1 - theta[35]/theta[9]) +
               (1 - theta[21]/theta[8]) )*
               (contrast == "VE" & whichmin_diff == 5) +
             ( theta[45] - (1 - theta[35]/theta[9]) +
               (1 - theta[23]/theta[8]) )*
               (contrast == "VE" & whichmin_diff == 6) +
             ( theta[45] - (1 - theta[35]/theta[11]) +
               (1 - theta[21]/theta[10]) )*
               (contrast == "VE" & whichmin_diff == 7) +
             ( theta[45] - (1 - theta[35]/theta[11]) +
               (1 - theta[23]/theta[10]) )*
               (contrast == "VE" & whichmin_diff == 8) +
             ( theta[45] - (1 - theta[37]/theta[13]) +
               (1 - theta[25]/theta[12]) )*
               (contrast == "VE" & whichmin_diff == 9) +
             ( theta[45] - (1 - theta[37]/theta[13]) +
               (1 - theta[27]/theta[12]) )*
               (contrast == "VE" & whichmin_diff == 10) +
             ( theta[45] - (1 - theta[37]/theta[15]) +
               (1 - theta[25]/theta[14]) )*
               (contrast == "VE" & whichmin_diff == 11) +
             ( theta[45] - (1 - theta[37]/theta[15]) +
               (1 - theta[27]/theta[14]) )*
               (contrast == "VE" & whichmin_diff == 12) +
             ( theta[45] - (1 - theta[39]/theta[13]) +
               (1 - theta[25]/theta[12]) )*
               (contrast == "VE" & whichmin_diff == 13) +
             ( theta[45] - (1 - theta[39]/theta[13]) +
               (1 - theta[27]/theta[12]) )*
               (contrast == "VE" & whichmin_diff == 14) +
             ( theta[45] - (1 - theta[39]/theta[15]) +
               (1 - theta[25]/theta[14]) )*
               (contrast == "VE" & whichmin_diff == 15) +
             ( theta[45] - (1 - theta[39]/theta[15]) +
               (1 - theta[27]/theta[14]) )*
               (contrast == "VE" & whichmin_diff == 16) +

             ( theta[45] - (theta[33] - theta[9]) +
               (theta[21] - theta[8]) )*
               (contrast == "Difference" & whichmin_diff == 1) +
             ( theta[45] - (theta[33] - theta[9]) +
               (theta[23] - theta[8]) )*
               (contrast == "Difference" & whichmin_diff == 2) +
             ( theta[45] - (theta[33] - theta[11]) +
               (theta[21] - theta[10]) )*
               (contrast == "Difference" & whichmin_diff == 3) +
             ( theta[45] - (theta[33] - theta[11]) +
               (theta[23] - theta[10]) )*
               (contrast == "Difference" & whichmin_diff == 4) +
             ( theta[45] - (theta[35] - theta[9]) +
               (theta[21] - theta[8]) )*
               (contrast == "Difference" & whichmin_diff == 5) +
             ( theta[45] - (theta[35] - theta[9]) +
               (theta[23] - theta[8]) )*
               (contrast == "Difference" & whichmin_diff == 6) +
             ( theta[45] - (theta[35] - theta[11]) +
               (theta[21] - theta[10]) )*
               (contrast == "Difference" & whichmin_diff == 7) +
             ( theta[45] - (theta[35] - theta[11]) +
               (theta[23] - theta[10]) )*
               (contrast == "Difference" & whichmin_diff == 8) +
             ( theta[45] - (theta[37] - theta[13]) +
               (theta[25] - theta[12]) )*
               (contrast == "Difference" & whichmin_diff == 9) +
             ( theta[45] - (theta[37] - theta[13]) +
               (theta[27] - theta[12]) )*
               (contrast == "Difference" & whichmin_diff == 10) +
             ( theta[45] - (theta[37] - theta[15]) +
               (theta[25] - theta[14]) )*
               (contrast == "Difference" & whichmin_diff == 11) +
             ( theta[45] - (theta[37] - theta[15]) +
               (theta[27] - theta[14]) )*
               (contrast == "Difference" & whichmin_diff == 12) +
             ( theta[45] - (theta[39] - theta[13]) +
               (theta[25] - theta[12]) )*
               (contrast == "Difference" & whichmin_diff == 13) +
             ( theta[45] - (theta[39] - theta[13]) +
               (theta[27] - theta[12]) )*
               (contrast == "Difference" & whichmin_diff == 14) +
             ( theta[45] - (theta[39] - theta[15]) +
               (theta[25] - theta[14]) )*
               (contrast == "Difference" & whichmin_diff == 15) +
             ( theta[45] - (theta[39] - theta[15]) +
               (theta[27] - theta[14]) )*
               (contrast == "Difference" & whichmin_diff == 16) +

             ( theta[45] - (log(theta[33]) - log(theta[9])) +
               (log(theta[21]) - log(theta[8])) )*
               (contrast == "logRR" & whichmin_diff == 1) +
             ( theta[45] - (log(theta[33]) - log(theta[9])) +
               (log(theta[23]) - log(theta[8])) )*
               (contrast == "logRR" & whichmin_diff == 2) +
             ( theta[45] - (log(theta[33]) - log(theta[11])) +
               (log(theta[21]) - log(theta[10])) )*
               (contrast == "logRR" & whichmin_diff == 3) +
             ( theta[45] - (log(theta[33]) - log(theta[11])) +
               (log(theta[23]) - log(theta[10])) )*
               (contrast == "logRR" & whichmin_diff == 4) +
             ( theta[45] - (log(theta[35]) - log(theta[9])) +
               (log(theta[21]) - log(theta[8])) )*
               (contrast == "logRR" & whichmin_diff == 5) +
             ( theta[45] - (log(theta[35]) - log(theta[9])) +
               (log(theta[23]) - log(theta[8])) )*
               (contrast == "logRR" & whichmin_diff == 6) +
             ( theta[45] - (log(theta[35]) - log(theta[11])) +
               (log(theta[21]) - log(theta[10])) )*
               (contrast == "logRR" & whichmin_diff == 7) +
             ( theta[45] - (log(theta[35]) - log(theta[11])) +
               (log(theta[23]) - log(theta[10])) )*
               (contrast == "logRR" & whichmin_diff == 8) +
             ( theta[45] - (log(theta[37]) - log(theta[13])) +
               (log(theta[25]) - log(theta[12])) )*
               (contrast == "logRR" & whichmin_diff == 9) +
             ( theta[45] - (log(theta[37]) - log(theta[13])) +
               (log(theta[27]) - log(theta[12])) )*
               (contrast == "logRR" & whichmin_diff == 10) +
             ( theta[45] - (log(theta[37]) - log(theta[15])) +
               (log(theta[25]) - log(theta[14])) )*
               (contrast == "logRR" & whichmin_diff == 11) +
             ( theta[45] - (log(theta[37]) - log(theta[15])) +
               (log(theta[27]) - log(theta[14])) )*
               (contrast == "logRR" & whichmin_diff == 12) +
             ( theta[45] - (log(theta[39]) - log(theta[13])) +
               (log(theta[25]) - log(theta[12])) )*
               (contrast == "logRR" & whichmin_diff == 13) +
             ( theta[45] - (log(theta[39]) - log(theta[13])) +
               (log(theta[27]) - log(theta[12])) )*
               (contrast == "logRR" & whichmin_diff == 14) +
             ( theta[45] - (log(theta[39]) - log(theta[15])) +
               (log(theta[25]) - log(theta[14])) )*
               (contrast == "logRR" & whichmin_diff == 15) +
             ( theta[45] - (log(theta[39]) - log(theta[15])) +
               (log(theta[27]) - log(theta[14])) )*
               (contrast == "logRR" & whichmin_diff == 16),

           # Upper bound of CEP(1, 0) - CEP(0, 0) = theta[45]
           ( theta[46] - (1 - theta[33]/theta[9]) +
               (1 - theta[21]/theta[8]) )*
               (contrast == "VE" & whichmax_diff == 1) +
             ( theta[46] - (1 - theta[33]/theta[9]) +
               (1 - theta[23]/theta[8]) )*
               (contrast == "VE" & whichmax_diff == 2) +
             ( theta[46] - (1 - theta[33]/theta[11]) +
               (1 - theta[21]/theta[10]) )*
               (contrast == "VE" & whichmax_diff == 3) +
             ( theta[46] - (1 - theta[33]/theta[11]) +
               (1 - theta[23]/theta[10]) )*
               (contrast == "VE" & whichmax_diff == 4) +
             ( theta[46] - (1 - theta[35]/theta[9]) +
               (1 - theta[21]/theta[8]) )*
               (contrast == "VE" & whichmax_diff == 5) +
             ( theta[46] - (1 - theta[35]/theta[9]) +
               (1 - theta[23]/theta[8]) )*
               (contrast == "VE" & whichmax_diff == 6) +
             ( theta[46] - (1 - theta[35]/theta[11]) +
               (1 - theta[21]/theta[10]) )*
               (contrast == "VE" & whichmax_diff == 7) +
             ( theta[46] - (1 - theta[35]/theta[11]) +
               (1 - theta[23]/theta[10]) )*
               (contrast == "VE" & whichmax_diff == 8) +
             ( theta[46] - (1 - theta[37]/theta[13]) +
               (1 - theta[25]/theta[12]) )*
               (contrast == "VE" & whichmax_diff == 9) +
             ( theta[46] - (1 - theta[37]/theta[13]) +
               (1 - theta[27]/theta[12]) )*
               (contrast == "VE" & whichmax_diff == 10) +
             ( theta[46] - (1 - theta[37]/theta[15]) +
               (1 - theta[25]/theta[14]) )*
               (contrast == "VE" & whichmax_diff == 11) +
             ( theta[46] - (1 - theta[37]/theta[15]) +
               (1 - theta[27]/theta[14]) )*
               (contrast == "VE" & whichmax_diff == 12) +
             ( theta[46] - (1 - theta[39]/theta[13]) +
               (1 - theta[25]/theta[12]) )*
               (contrast == "VE" & whichmax_diff == 13) +
             ( theta[46] - (1 - theta[39]/theta[13]) +
               (1 - theta[27]/theta[12]) )*
               (contrast == "VE" & whichmax_diff == 14) +
             ( theta[46] - (1 - theta[39]/theta[15]) +
               (1 - theta[25]/theta[14]) )*
               (contrast == "VE" & whichmax_diff == 15) +
             ( theta[46] - (1 - theta[39]/theta[15]) +
               (1 - theta[27]/theta[14]) )*
               (contrast == "VE" & whichmax_diff == 16) +

             ( theta[46] - (theta[33] - theta[9]) +
               (theta[21] - theta[8]) )*
               (contrast == "Difference" & whichmax_diff == 1) +
             ( theta[46] - (theta[33] - theta[9]) +
               (theta[23] - theta[8]) )*
               (contrast == "Difference" & whichmax_diff == 2) +
             ( theta[46] - (theta[33] - theta[11]) +
               (theta[21] - theta[10]) )*
               (contrast == "Difference" & whichmax_diff == 3) +
             ( theta[46] - (theta[33] - theta[11]) +
               (theta[23] - theta[10]) )*
               (contrast == "Difference" & whichmax_diff == 4) +
             ( theta[46] - (theta[35] - theta[9]) +
               (theta[21] - theta[8]) )*
               (contrast == "Difference" & whichmax_diff == 5) +
             ( theta[46] - (theta[35] - theta[9]) +
               (theta[23] - theta[8]) )*
               (contrast == "Difference" & whichmax_diff == 6) +
             ( theta[46] - (theta[35] - theta[11]) +
               (theta[21] - theta[10]) )*
               (contrast == "Difference" & whichmax_diff == 7) +
             ( theta[46] - (theta[35] - theta[11]) +
               (theta[23] - theta[10]) )*
               (contrast == "Difference" & whichmax_diff == 8) +
             ( theta[46] - (theta[37] - theta[13]) +
               (theta[25] - theta[12]) )*
               (contrast == "Difference" & whichmax_diff == 9) +
             ( theta[46] - (theta[37] - theta[13]) +
               (theta[27] - theta[12]) )*
               (contrast == "Difference" & whichmax_diff == 10) +
             ( theta[46] - (theta[37] - theta[15]) +
               (theta[25] - theta[14]) )*
               (contrast == "Difference" & whichmax_diff == 11) +
             ( theta[46] - (theta[37] - theta[15]) +
               (theta[27] - theta[14]) )*
               (contrast == "Difference" & whichmax_diff == 12) +
             ( theta[46] - (theta[39] - theta[13]) +
               (theta[25] - theta[12]) )*
               (contrast == "Difference" & whichmax_diff == 13) +
             ( theta[46] - (theta[39] - theta[13]) +
               (theta[27] - theta[12]) )*
               (contrast == "Difference" & whichmax_diff == 14) +
             ( theta[46] - (theta[39] - theta[15]) +
               (theta[25] - theta[14]) )*
               (contrast == "Difference" & whichmax_diff == 15) +
             ( theta[46] - (theta[39] - theta[15]) +
               (theta[27] - theta[14]) )*
               (contrast == "Difference" & whichmax_diff == 16) +

             ( theta[46] - (log(theta[33]) - log(theta[9])) +
               (log(theta[21]) - log(theta[8])) )*
               (contrast == "logRR" & whichmax_diff == 1) +
             ( theta[46] - (log(theta[33]) - log(theta[9])) +
               (log(theta[23]) - log(theta[8])) )*
               (contrast == "logRR" & whichmax_diff == 2) +
             ( theta[46] - (log(theta[33]) - log(theta[11])) +
               (log(theta[21]) - log(theta[10])) )*
               (contrast == "logRR" & whichmax_diff == 3) +
             ( theta[46] - (log(theta[33]) - log(theta[11])) +
               (log(theta[23]) - log(theta[10])) )*
               (contrast == "logRR" & whichmax_diff == 4) +
             ( theta[46] - (log(theta[35]) - log(theta[9])) +
               (log(theta[21]) - log(theta[8])) )*
               (contrast == "logRR" & whichmax_diff == 5) +
             ( theta[46] - (log(theta[35]) - log(theta[9])) +
               (log(theta[23]) - log(theta[8])) )*
               (contrast == "logRR" & whichmax_diff == 6) +
             ( theta[46] - (log(theta[35]) - log(theta[11])) +
               (log(theta[21]) - log(theta[10])) )*
               (contrast == "logRR" & whichmax_diff == 7) +
             ( theta[46] - (log(theta[35]) - log(theta[11])) +
               (log(theta[23]) - log(theta[10])) )*
               (contrast == "logRR" & whichmax_diff == 8) +
             ( theta[46] - (log(theta[37]) - log(theta[13])) +
               (log(theta[25]) - log(theta[12])) )*
               (contrast == "logRR" & whichmax_diff == 9) +
             ( theta[46] - (log(theta[37]) - log(theta[13])) +
               (log(theta[27]) - log(theta[12])) )*
               (contrast == "logRR" & whichmax_diff == 10) +
             ( theta[46] - (log(theta[37]) - log(theta[15])) +
               (log(theta[25]) - log(theta[14])) )*
               (contrast == "logRR" & whichmax_diff == 11) +
             ( theta[46] - (log(theta[37]) - log(theta[15])) +
               (log(theta[27]) - log(theta[14])) )*
               (contrast == "logRR" & whichmax_diff == 12) +
             ( theta[46] - (log(theta[39]) - log(theta[13])) +
               (log(theta[25]) - log(theta[12])) )*
               (contrast == "logRR" & whichmax_diff == 13) +
             ( theta[46] - (log(theta[39]) - log(theta[13])) +
               (log(theta[27]) - log(theta[12])) )*
               (contrast == "logRR" & whichmax_diff == 14) +
             ( theta[46] - (log(theta[39]) - log(theta[15])) +
               (log(theta[25]) - log(theta[14])) )*
               (contrast == "logRR" & whichmax_diff == 15) +
             ( theta[46] - (log(theta[39]) - log(theta[15])) +
               (log(theta[27]) - log(theta[14])) )*
               (contrast == "logRR" & whichmax_diff == 16),

           # pi = theta[47]
           (1-Y)*(R-theta[47])

    ))
  }

}
