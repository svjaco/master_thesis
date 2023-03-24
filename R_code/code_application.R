# Code for Part II: Application

# Master thesis 
# Sven Jacobs
# Winter 2022, M.Sc. Economics, Bonn University
# Supervisor: Prof. Dr. Dominik Liebl

################################################################################

rm(list = ls())

library(haven)
library(tidyverse)

library(rddensity)
library(rdrobust)
if (!requireNamespace("remotes")) {
    install.packages("remotes")
}
remotes::install_github("kolesarm/RDHonest")
library(RDHonest)

################################################################################

##############
#    Data    #
##############

# Loading data
data_original <- read_stata("../data/LM2007.dta")

# Outliers
# Excluding two counties with distinctly higher mortality rate
data <- data_original[data_original$mort_age59_related_postHS <= 100, ]

# Change column to percentage (like others)
data$census1960_pctsch534 <- data$census1960_pctsch534 * 100

# Renaming
cutoff <- 59.1984
Y <- data$mort_age59_related_postHS
X <- data$povrate60
D <- X >= cutoff

###################
#    Figure A3    #
###################

ggplot(data.frame(X, Y), aes(X, Y)) +
    geom_point(col = "darkblue") +
    geom_vline(xintercept = cutoff) +
    xlab("Poverty rate") +
    ylab("Mortality rate per 100000, HS related causes, Ages 5–9, 1973–83") +
    theme_bw(base_size = 18)

ggsave(filename = "../figures/figure_A03.pdf",
       width = 15, height = 8, units = "in", dpi = 600)

##################
#    Figure 4    #
##################

rdplot_Y <- rdplot(Y, X, c = cutoff, p = 3,
                   title = "", x.label = "Poverty rate",
                   y.label = "Mortality rate per 100000, HS related causes, Ages 5–9, 1973–83")$rdplot
rdplot_Y + theme_bw(base_size = 18)

ggsave(filename = "../figures/figure_04.pdf",
       width = 15, height = 8, units = "in", dpi = 600)

###################
#    Figure A4    #
###################

ggplot(data.frame(X, D), aes(x = X, fill = D)) +
    geom_histogram(boundary = cutoff, binwidth = 1.5, col = "black") +
    geom_vline(xintercept = cutoff) +
    labs(x = "Poverty rate", y = "Number of counties") +
    scale_fill_manual(values = c("blue", "red")) +
    guides(fill = "none") +
    theme_bw(base_size = 18)

ggsave(filename = "../figures/figure_A04.pdf",
       width = 15, height = 8, units = "in", dpi = 600)

##################
#    Figure 5    #
##################

density_test <- rddensity(X, c = cutoff)
summary(density_test)

density_plot <- rdplotdensity(rdd = density_test, X = X, type = "both", xlabel = "Poverty rate",
                              lcol = c(4, 2))$Estplot
density_plot + theme_bw(base_size = 18) +
    theme(plot.title = element_blank(), axis.title.y = element_blank(), legend.position = "none",
          axis.text.y = element_text(angle = 90, hjust = 0.5))

ggsave(filename = "../figures/figure_05.pdf",
       width = 15, height = 8, units = "in", dpi = 600)

####################
#    Figure A5a    #
####################

injury <- data$mort_age59_injury_postHS

rdplot_injury <- rdplot(injury, X, c = cutoff, p = 3,
                        title = "", x.label = "Poverty rate",
                        y.label = "Mortality rate per 100000, Injuries, Ages 5–9, 1973–83")$rdplot
rdplot_injury + theme_bw(base_size = 18) + theme(plot.title = element_blank())

ggsave(filename = "../figures/figure_A05a.pdf",
       width = 8, height = 8, units = "in", dpi = 600)

####################
#    Figure A6a    #
####################

h_CCT_injury <- rdrobust(injury, X, c = cutoff)$bws[1, 1]
llplot_injury <- rdplot(injury[abs(X - cutoff) <= h_CCT_injury], X[abs(X - cutoff) <= h_CCT_injury],
                        c = cutoff, p = 1, kernel = "triangular",
                        title = "", x.label = "Poverty rate",
                        y.label = "Mortality rate per 100000, Injuries, Ages 5–9, 1973–83")$rdplot
llplot_injury + theme_bw(base_size = 18) + theme(plot.title = element_blank())

ggsave(filename = "../figures/figure_A06a.pdf",
       width = 8, height = 8, units = "in", dpi = 600)

####################
#    Figure A5b    #
####################

Y_pre <- data$mort_age59_related_preHS

rdplot_Y_pre <- rdplot(Y_pre, X, c = cutoff, p = 3,
                       title = "", x.label = "Poverty rate",
                       y.label = "Mortality rate per 100000, HS related causes, Ages 5–9, 1959–64")$rdplot
rdplot_Y_pre + theme_bw(base_size = 18) + theme(plot.title = element_blank())

ggsave(filename = "../figures/figure_A05b.pdf",
       width = 8, height = 8, units = "in", dpi = 600)

####################
#    Figure A6b    #
####################

h_CCT_Y_pre <- rdrobust(Y_pre, X, c = cutoff)$bws[1, 1]
llplot_Y_pre <- rdplot(Y_pre[abs(X - cutoff) <= h_CCT_Y_pre], X[abs(X - cutoff) <= h_CCT_Y_pre],
                       c = cutoff, p = 1, kernel = "triangular",
                       title = "", x.label = "Poverty rate",
                       y.label = "Mortality rate per 100000, HS related causes, Ages 5–9, 1959–64")$rdplot
llplot_Y_pre + theme_bw(base_size = 18) + theme(plot.title = element_blank())

ggsave(filename = "../figures/figure_A06b.pdf",
       width = 8, height = 8, units = "in", dpi = 600)

##################
#    Table A2    #
##################

summary(rdrobust(data$mort_age59_injury_postHS, X, c = cutoff))
summary(rdrobust(data$mort_age59_all_postHS, X, c = cutoff))
summary(rdrobust(data$mort_age25plus_related_postHS, X, c = cutoff))
summary(rdrobust(data$mort_age25plus_injuries_postHS, X, c = cutoff))
summary(rdrobust(data$mort_age59_related_preHS, X, c = cutoff))
summary(rdrobust(data$census1960_pop, X, c = cutoff))
summary(rdrobust(data$census1960_pctsch1417, X, c = cutoff))
summary(rdrobust(data$census1960_pctsch534, X, c = cutoff))
summary(rdrobust(data$census1960_pctsch25plus, X, c = cutoff))
summary(rdrobust(data$census1960_pop1417, X, c = cutoff))
summary(rdrobust(data$census1960_pop534, X, c = cutoff))
summary(rdrobust(data$census1960_pop25plus, X, c = cutoff))
summary(rdrobust(data$census1960_pcturban, X, c = cutoff))
summary(rdrobust(data$census1960_pctblack, X, c = cutoff))

##################
#    Table A3    #
##################

summary(rdrobust(Y[X < cutoff], X[X < cutoff], c = median(X[X < cutoff])))
summary(rdrobust(Y[X < cutoff], X[X < cutoff], c = 40))
summary(rdrobust(Y[X < cutoff], X[X < cutoff], c = 50))

summary(rdrobust(Y, X, c = cutoff))

summary(rdrobust(Y[X >= cutoff], X[X >= cutoff], c = median(X[X >= cutoff])))
summary(rdrobust(Y[X >= cutoff], X[X >= cutoff], c = 70))

##################
#    Table A4    #
##################

summary(rdrobust(Y[abs(X - cutoff) >= 0], X[abs(X - cutoff) >= 0], c = cutoff))
summary(rdrobust(Y[abs(X - cutoff) >= 0.1], X[abs(X - cutoff) >= 0.1], c = cutoff))
summary(rdrobust(Y[abs(X - cutoff) >= 0.2], X[abs(X - cutoff) >= 0.2], c = cutoff))
summary(rdrobust(Y[abs(X - cutoff) >= 0.3], X[abs(X - cutoff) >= 0.3], c = cutoff))
summary(rdrobust(Y[abs(X - cutoff) >= 0.4], X[abs(X - cutoff) >= 0.4], c = cutoff))
summary(rdrobust(Y[abs(X - cutoff) >= 0.5], X[abs(X - cutoff) >= 0.5], c = cutoff))

#################
#    Table 2    #
#################

summary(rdrobust(Y, X, c = cutoff))
summary(rdrobust(Y, X, c = cutoff, bwselect = "cerrd"))


US_half <- rdrobust(Y, X, c = cutoff)$bws[1, 1] / 2
US_third <- rdrobust(Y, X, c = cutoff)$bws[1, 1] / 3
summary(rdrobust(Y, X, c = cutoff, h = US_half))
summary(rdrobust(Y, X, c = cutoff, h = US_third))


census_1960 <- cbind(data$census1960_pop,
                     data$census1960_pctsch1417,
                     data$census1960_pctsch534, 
                     data$census1960_pctsch25plus,
                     data$census1960_pop1417,
                     data$census1960_pop534,
                     data$census1960_pop25plus,
                     data$census1960_pcturban,
                     data$census1960_pctblack)
summary(rdrobust(Y, X, c = cutoff, covs = census_1960))


p_value_AK <- function(RDHonest_out) {
    
    estimate_se <- RDHonest_out$coefficients$estimate / RDHonest_out$coefficients$std.error
    maxbias_se <- RDHonest_out$coefficients$maximum.bias / RDHonest_out$coefficients$std.error
    
    fun_search <- function(alpha) {
        abs(estimate_se) - RDHonest::CVb(maxbias_se, alpha)
    }
    
    uniroot(fun_search, interval = c(0.001, 0.999))$root
    
}

HS_AK_MSE <- RDHonest(Y ~ X, cutoff = cutoff, opt.criterion = "MSE")
p_value_AK(HS_AK_MSE)
HS_AK_FLCI <- RDHonest(Y ~ X, cutoff = cutoff, opt.criterion = "FLCI")
p_value_AK(HS_AK_FLCI)


summary(rdrobust(Y, X, c = cutoff, h = 9))
summary(rdrobust(Y, X, c = cutoff, h = 18))

HS_AK_h9 <- RDHonest(Y ~ X, cutoff = cutoff, h = 9)
p_value_AK(HS_AK_h9)
HS_AK_h18 <- RDHonest(Y ~ X, cutoff = cutoff, h = 18)
# p_value_AK(HS_AK_h18)

##################
#    Figure 6    #
##################

h_vector <- c(3.457, 2.304, 6.913, 4.650, 7.115, 9.000, 4.551, 4.670)
estimate_vector <- c(-3.428, -2.590, -2.389, -3.248, -2.445, -2.182, -3.283, -3.239)
ci_lower_vector <- c(-5.856, -5.488, -5.426, -6.092, -5.155, -5.722, -6.039, -6.022)
ci_upper_vector <- c(-1.001, 0.307, -0.083, -0.749, -0.339, -0.350, -0.526, -0.456)
labels <- c("US, 1/2", "US, 1/3", "RBC, MSE", "RBC, CE", "RBC, MSE, w/ cov.",
            "RBC, LM1", "AK, minimax MSE", "AK, CIL")
df <- data.frame(h_vector, estimate_vector, ci_lower_vector, ci_upper_vector, labels)

ggplot(df) +
    geom_errorbar(aes(x = h_vector, ymin = ci_lower_vector, ymax = ci_upper_vector),
                  width = 0.1, col = "darkblue", linewidth = 1) +
    geom_point(aes(x = h_vector, y = estimate_vector), size = 3, col = "red") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
    labs(x = "Bandwidth", y = "RD inference") +
    geom_text(data = df[-c(4, 8), ], aes(x = h_vector - 0.1, y = ci_lower_vector + 0.2, label = labels),
              angle = 90, hjust = 0, size = 6) +
    geom_text(data = df[4, ], aes(x = h_vector + 0.22, y = ci_lower_vector - 0.2, label = labels), size = 6) +
    geom_text(data = df[8, ], aes(x = h_vector + 0.32, y = ci_upper_vector, label = labels), size = 6) +
    theme_bw(base_size = 18)

ggsave(filename = "../figures/figure_06.pdf",
       width = 15, height = 8, units = "in", dpi = 600)