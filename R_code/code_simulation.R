# Code for Part III: Simulation

# Master thesis 
# Sven Jacobs
# Winter 2022, M.Sc. Economics, Bonn University
# Supervisor: Prof. Dr. Dominik Liebl

################################################################################

rm(list = ls())

library(tidyverse)

library(rdrobust)
if (!requireNamespace("remotes")) {
    install.packages("remotes")
}
remotes::install_github("kolesarm/RDHonest")
library(RDHonest)

################################################################################

##################
#    Figure 8    #
##################

x_beta <- seq(0, 1, 0.01)
pdf_beta <- dbeta(x_beta, 2, 4)
df_beta <- data.frame(x_beta, pdf_beta)

ggplot(df_beta) +
    geom_line(aes(x = x_beta, y = pdf_beta)) +
    labs(x = "", y = "PDF") +
    theme_bw(base_size = 20) +
    theme(axis.title.x = element_blank())

ggsave(filename = "../figures/figure_08.pdf",
       width = 15, height = 8, units = "in", dpi = 600)

################################################################################
#                              Monte Carlo Study                               #
################################################################################

model1 <- function(x) {
    
    ifelse(x < 0, 0.48 + 1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5,
           0.52 + 0.84*x - 3.00*x^2 + 7.99*x^3 - 9.01*x^4 + 3.56*x^5)
    
}

model2 <- function(x) {
    
    ifelse(x < 0, 3.71 + 2.30*x + 3.28*x^2 + 1.45*x^3 + 0.23*x^4 + 0.03*x^5,
           0.26 + 18.49*x - 54.81*x^2 + 74.30*x^3 - 45.02*x^4 + 9.83*x^5)
    
}

model3 <- function(x) {
    
    ifelse(x < 0, 0.48 + 1.27*x - 0.5*7.18*x^2 + 0.7*20.21*x^3 + 1.1*21.54*x^4 + 1.5*7.33*x^5,
           0.52 + 0.84*x - 0.1*3.00*x^2 - 0.3*7.99*x^3 - 0.1*9.01*x^4 + 3.56*x^5)
    
}


generate_data <- function(model) {
    
    X <- 2*rbeta(500, 2, 4) - 1
    epsilon <- rnorm(500, 0, 0.1295)
    Y <- model(X) + epsilon
    
    return(data.frame(X, Y))
    
}


# Population RD ATE
tau_model1 <- 0.04
tau_model2 <- -3.45
tau_model3 <- 0.04

# True smoothness constants M (2nd-order Hölder class)
# Computed with WolframAlpha
M_model1 <- 14.36
M_model2 <- 109.62
M_model3 <- 45.406

# Number of Monte Carlo replications
sim_rep <- 5000    

################################################################################

# Note:

# The simulation code is divided into two parts.
# Part one for the approaches "Robust bias-correction" and
# "Armstrong and Kolesár", part two for "Undersmoothing".
# The reason is that undersmoothing (due to a small bandwidth)
# can cause numerical instability during the simulation.
# For Model 2 (Ludwig and Miller, 2007) results cannot be obtained
# for undersmoothing.

################################################################################

################################
#    Without undersmoothing    #
################################
    
simulate_model <- function(model, sim_rep) {
        
    if(identical(model, model1)) {
        
        tau <- tau_model1
        M_true <- M_model1
        
    } else if(identical(model, model2)) {
        
        tau <- tau_model2
        M_true <- M_model2 
        
    } else if(identical(model, model3)) {
        
        tau <- tau_model3
        M_true <- M_model3
        
    }
    
    
    h_MSE_vec <- rep(NA, sim_rep)
    h_CE_vec <- rep(NA, sim_rep)
    h_MSE_M_ROT <- rep(NA, sim_rep)
    h_MSE_M <- rep(NA, sim_rep)
    h_CIL_M_ROT <- rep(NA, sim_rep)
    h_CIL_M <- rep(NA, sim_rep)
    
    tau_MSE_vec <- rep(NA, sim_rep)
    tau_CE_vec <- rep(NA, sim_rep)
    tau_MSE_M_ROT <- rep(NA, sim_rep)
    tau_MSE_M <- rep(NA, sim_rep)
    tau_CIL_M_ROT <- rep(NA, sim_rep)
    tau_CIL_M <- rep(NA, sim_rep)
    
    SErr_MSE_vec <- rep(NA, sim_rep)
    SErr_CE_vec <- rep(NA, sim_rep)
    SErr_MSE_M_ROT <- rep(NA, sim_rep)
    SErr_MSE_M <- rep(NA, sim_rep)
    SErr_CIL_M_ROT <- rep(NA, sim_rep)
    SErr_CIL_M <- rep(NA, sim_rep)
    
    cov_MSE_vec <- rep(NA, sim_rep)
    cov_CE_vec <- rep(NA, sim_rep)
    cov_MSE_M_ROT <- rep(NA, sim_rep)
    cov_MSE_M <- rep(NA, sim_rep)
    cov_CIL_M_ROT <- rep(NA, sim_rep)
    cov_CIL_M <- rep(NA, sim_rep)
    
    length_MSE_vec <- rep(NA, sim_rep)
    length_CE_vec <- rep(NA, sim_rep)
    length_MSE_M_ROT <- rep(NA, sim_rep)
    length_MSE_M <- rep(NA, sim_rep)
    length_CIL_M_ROT <- rep(NA, sim_rep)
    length_CIL_M <- rep(NA, sim_rep)
    
    M_ROT_vec <- rep(NA, sim_rep)
    

    for (i in 1:sim_rep) {
        
        data <- generate_data(model = model)
        X <- data$X
        Y <- data$Y
    
        
        MSE_out <- rdrobust(Y, X)
        h_MSE_vec[i] <- MSE_out$bws[1, 1]
        tau_MSE_vec[i] <- MSE_out$coef[1]
        SErr_MSE_vec[i] <- (MSE_out$coef[1] - tau)^2
        MSE_ci_low <- MSE_out$ci[3, 1]
        MSE_ci_high <- MSE_out$ci[3, 2]
        cov_MSE_vec[i] <- tau >= MSE_ci_low & tau <= MSE_ci_high
        length_MSE_vec[i] <- MSE_ci_high - MSE_ci_low
        
        CE_out <- rdrobust(Y, X, bwselect = "cerrd")
        h_CE_vec[i] <- CE_out$bws[1, 1]
        tau_CE_vec[i] <- CE_out$coef[1]
        SErr_CE_vec[i] <- (CE_out$coef[1] - tau)^2
        CE_ci_low <- CE_out$ci[3, 1]
        CE_ci_high <- CE_out$ci[3, 2]
        cov_CE_vec[i] <- tau >= CE_ci_low & tau <= CE_ci_high
        length_CE_vec[i] <- CE_ci_high - CE_ci_low
        
        
        MSE_M_ROT_out <- RDHonest(Y ~ X, opt.criterion = "MSE")$coefficients
        M_ROT_vec[i] <- MSE_M_ROT_out$M
        h_MSE_M_ROT[i] <- MSE_M_ROT_out$bandwidth 
        tau_MSE_M_ROT[i] <- MSE_M_ROT_out$estimate
        SErr_MSE_M_ROT[i] <- (MSE_M_ROT_out$estimate - tau)^2
        MSE_M_ROT_ci_low <- MSE_M_ROT_out$conf.low
        MSE_M_ROT_ci_high <- MSE_M_ROT_out$conf.high
        cov_MSE_M_ROT[i] <- tau >= MSE_M_ROT_ci_low & tau <= MSE_M_ROT_ci_high
        length_MSE_M_ROT[i] <- MSE_M_ROT_ci_high - MSE_M_ROT_ci_low
        
        MSE_M_out <- RDHonest(Y ~ X, opt.criterion = "MSE", M = M_true)$coefficients
        h_MSE_M[i] <- MSE_M_out$bandwidth 
        tau_MSE_M[i] <- MSE_M_out$estimate
        SErr_MSE_M[i] <- (MSE_M_out$estimate - tau)^2
        MSE_M_ci_low <- MSE_M_out$conf.low
        MSE_M_ci_high <- MSE_M_out$conf.high
        cov_MSE_M[i] <- tau >= MSE_M_ci_low & tau <= MSE_M_ci_high
        length_MSE_M[i] <- MSE_M_ci_high - MSE_M_ci_low
        
        CIL_M_ROT_out <- RDHonest(Y ~ X, opt.criterion = "FLCI")$coefficients
        h_CIL_M_ROT[i] <- CIL_M_ROT_out$bandwidth 
        tau_CIL_M_ROT[i] <- CIL_M_ROT_out$estimate
        SErr_CIL_M_ROT[i] <- (CIL_M_ROT_out$estimate - tau)^2
        CIL_M_ROT_ci_low <- CIL_M_ROT_out$conf.low
        CIL_M_ROT_ci_high <- CIL_M_ROT_out$conf.high
        cov_CIL_M_ROT[i] <- tau >= CIL_M_ROT_ci_low & tau <= CIL_M_ROT_ci_high
        length_CIL_M_ROT[i] <- CIL_M_ROT_ci_high - CIL_M_ROT_ci_low
        
        CIL_M_out <- RDHonest(Y ~ X, opt.criterion = "FLCI", M = M_true)$coefficients
        h_CIL_M[i] <- CIL_M_out$bandwidth 
        tau_CIL_M[i] <- CIL_M_out$estimate
        SErr_CIL_M[i] <- (CIL_M_out$estimate - tau)^2
        CIL_M_ci_low <- CIL_M_out$conf.low
        CIL_M_ci_high <- CIL_M_out$conf.high
        cov_CIL_M[i] <- tau >= CIL_M_ci_low & tau <= CIL_M_ci_high
        length_CIL_M[i] <- CIL_M_ci_high - CIL_M_ci_low

    }
    
    
    df <- data.frame(matrix(nrow = 6, ncol = 6))
    df_names <- c("Approach", "Bandwidth", "Bias", "RMSE", "CI_coverage", "CI_length")
    colnames(df) <- df_names
    
    df$Approach <- c("RBC MSE", "RBC CE",
                     "AK MSE_M_ROT", "AK MSE_M", "AK CIL_M_ROT", "AK CIL_M")
    
    df$Bandwidth <- c(mean(h_MSE_vec), mean(h_CE_vec),
                      mean(h_MSE_M_ROT), mean(h_MSE_M), mean(h_CIL_M_ROT), mean(h_CIL_M))
    
    df$Bias <- c(mean(tau_MSE_vec), mean(tau_CE_vec),
                 mean(tau_MSE_M_ROT), mean(tau_MSE_M), mean(tau_CIL_M_ROT), mean(tau_CIL_M))
    df$Bias <- df$Bias - tau

    df$RMSE <- c(mean(SErr_MSE_vec), mean(SErr_CE_vec),
                 mean(SErr_MSE_M_ROT), mean(SErr_MSE_M), mean(SErr_CIL_M_ROT), mean(SErr_CIL_M))
    df$RMSE <- sqrt(df$RMSE)

    df$CI_coverage <- c(mean(cov_MSE_vec), mean(cov_CE_vec),
                        mean(cov_MSE_M_ROT), mean(cov_MSE_M), mean(cov_CIL_M_ROT), mean(cov_CIL_M))
    df$CI_coverage <- df$CI_coverage * 100

    df$CI_length <- c(mean(length_MSE_vec), mean(length_CE_vec),
                      mean(length_MSE_M_ROT), mean(length_MSE_M), mean(length_CIL_M_ROT), mean(length_CIL_M))
    
    
    cat("The average ROT M over the", sim_rep, "replications is:", mean(M_ROT_vec))
    
    return(df)
    
}
    

##########################################
#    Table 3 – Without undersmoothing    #
##########################################

set.seed(123)
simulate_model(model1, sim_rep)

##########################################
#    Table 4 – Without undersmoothing    #
##########################################

set.seed(123)
simulate_model(model2, sim_rep)

##########################################
#    Table 5 – Without undersmoothing    #
##########################################

set.seed(123)
simulate_model(model3, sim_rep)


########################
#    Undersmoothing    #
########################

simulate_model_US <- function(model, sim_rep) {
    
    if(identical(model, model1)) {
        
        tau <- tau_model1
        M_true <- M_model1
        
    } else if(identical(model, model2)) {
        
        tau <- tau_model2
        M_true <- M_model2 
        
    } else if(identical(model, model3)) {
        
        tau <- tau_model3
        M_true <- M_model3
        
    }
    
    
    h_MSE_vec <- rep(NA, sim_rep)
    
    tau_US2_vec <- rep(NA, sim_rep)
    tau_US3_vec <- rep(NA, sim_rep)
    tau_MSE_vec <- rep(NA, sim_rep)
    
    SErr_US2_vec <- rep(NA, sim_rep)
    SErr_US3_vec <- rep(NA, sim_rep)
    SErr_MSE_vec <- rep(NA, sim_rep)
    
    cov_US2_vec <- rep(NA, sim_rep)
    cov_US3_vec <- rep(NA, sim_rep)
    cov_MSE_vec <- rep(NA, sim_rep)
    
    length_US2_vec <- rep(NA, sim_rep)
    length_US3_vec <- rep(NA, sim_rep)
    length_MSE_vec <- rep(NA, sim_rep)
    
    
    for (i in 1:sim_rep) {
        
        data <- generate_data(model = model)
        X <- data$X
        Y <- data$Y
        
        
        MSE_out <- rdrobust(Y, X)
        h_MSE_vec[i] <- MSE_out$bws[1, 1]
        tau_MSE_vec[i] <- MSE_out$coef[1]
        SErr_MSE_vec[i] <- (MSE_out$coef[1] - tau)^2
        MSE_ci_low <- MSE_out$ci[3, 1]
        MSE_ci_high <- MSE_out$ci[3, 2]
        cov_MSE_vec[i] <- tau >= MSE_ci_low & tau <= MSE_ci_high
        length_MSE_vec[i] <- MSE_ci_high - MSE_ci_low
        
        
        US2_out <- rdrobust(Y, X, h = MSE_out$bws[1, 1]/2)
        tau_US2_vec[i] <- US2_out$coef[1]
        SErr_US2_vec[i] <- (US2_out$coef[1] - tau)^2
        US2_ci_low <- US2_out$ci[1, 1]
        US2_ci_high <- US2_out$ci[1, 2]
        cov_US2_vec[i] <- tau >= US2_ci_low & tau <= US2_ci_high
        length_US2_vec[i] <- US2_ci_high - US2_ci_low
        
        US3_out <- rdrobust(Y, X, h = MSE_out$bws[1, 1]/3)
        tau_US3_vec[i] <- US3_out$coef[1]
        SErr_US3_vec[i] <- (US3_out$coef[1] - tau)^2
        US3_ci_low <- US3_out$ci[1, 1]
        US3_ci_high <- US3_out$ci[1, 2]
        cov_US3_vec[i] <- tau >= US3_ci_low & tau <= US3_ci_high
        length_US3_vec[i] <- US3_ci_high - US3_ci_low
        
    }
    
    
    df <- data.frame(matrix(nrow = 3, ncol = 6))
    df_names <- c("Approach", "Bandwidth", "Bias", "RMSE", "CI_coverage", "CI_length")
    colnames(df) <- df_names
    
    df$Approach <- c("US 1/2", "US 1/3", "RBC MSE")
    
    df$Bandwidth <- c(mean(h_MSE_vec/2), mean(h_MSE_vec/3), mean(h_MSE_vec))
    
    df$Bias <- c(mean(tau_US2_vec), mean(tau_US3_vec), mean(tau_MSE_vec))
    df$Bias <- df$Bias - tau
    
    df$RMSE <- c(mean(SErr_US2_vec), mean(SErr_US3_vec), mean(SErr_MSE_vec))
    df$RMSE <- sqrt(df$RMSE)
    
    df$CI_coverage <- c(mean(cov_US2_vec), mean(cov_US3_vec), mean(cov_MSE_vec))
    df$CI_coverage <- df$CI_coverage * 100
    
    df$CI_length <- c(mean(length_US2_vec), mean(length_US3_vec), mean(length_MSE_vec))
    
    
    return(df)
    
}


##################################
#    Table 3 – Undersmoothing    #
##################################

set.seed(456)
simulate_model_US(model1, sim_rep)

##################################
#    Table 4 – Undersmoothing    #
##################################

# set.seed(123)
# simulate_model_US(model2, sim_rep)

##################################
#    Table 5 – Undersmoothing    #
##################################

set.seed(456)
simulate_model_US(model3, sim_rep)