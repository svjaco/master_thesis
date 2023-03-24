# Code for Part I: Theory  

# Master thesis 
# Sven Jacobs
# Winter 2022, M.Sc. Economics, Bonn University
# Supervisor: Prof. Dr. Dominik Liebl

################################################################################

rm(list = ls())

library(tidyverse)

################################################################################

##################
#    Figure 3    #
##################

set.seed(123)

control_fun <- function(x) {-0.5*x^3 + 0.75*x^2 + 0.3}
treatment_fun <- function(x) {-0.5*x^3 + 0.9*x^2 + 0.45}

n_control <- 50
n_treatment <- 50

cutoff <- 0.5

X_control <- seq(0, cutoff - 0.01, length.out = n_control)
X_treatment <- seq(cutoff, 1, length.out = n_treatment)

epsilon_control <- rnorm(n_control, sd = 0.03)
epsilon_treatment <- rnorm(n_treatment, sd = 0.03)

Y_control <- control_fun(X_control) + epsilon_control
Y_treatment <- treatment_fun(X_treatment) + epsilon_treatment

X_tilde_control <- X_control - cutoff
X_tilde_treatment <- X_treatment - cutoff

h <- 0.15

weights_control <- rep(NA, n_control)
weights_treatment <- rep(NA, n_treatment)
weights_control <- ifelse(abs(X_tilde_control/h) <= 1, 1 - abs(X_tilde_control/h), 0)
weights_treatment <- ifelse(abs(X_tilde_treatment/h) <= 1, 1 - abs(X_tilde_treatment/h), 0)

fit_control <- lm(Y_control ~ X_tilde_control, weights = weights_control)
fit_treatment <- lm(Y_treatment ~ X_tilde_treatment, weights = weights_treatment)

intercept_control <- fit_control$coefficients[[1]]
slope_control <- fit_control$coefficients[[2]]
intercept_treatment <- fit_treatment$coefficients[[1]]
slope_treatment <- fit_treatment$coefficients[[2]]
coordinates_fits <- tibble(
    x = c(0, 0),
    y = c(intercept_control, intercept_treatment),
    xend = c(-h, h),
    yend = c(intercept_control - h*slope_control, intercept_treatment + h*slope_treatment)
)

X_tilde <- c(X_tilde_control, X_tilde_treatment)
Y <- c(Y_control, Y_treatment) 
treatment <- ifelse(X_tilde >= 0, "treatment", "control")
weights <- c(weights_control, weights_treatment)

data <- tibble(treatment, X_tilde, Y, weights)

control_fun_plot <- function(x) {-0.5*(x + cutoff)^3 + 0.75*(x + cutoff)^2 + 0.3}
treatment_fun_plot <- function(x) {-0.5*(x + cutoff)^3 + 0.9*(x + cutoff)^2 + 0.45}

ll_estimation_plot <- ggplot(data, aes(x = X_tilde, y = Y)) +
    geom_point(data = filter(data, abs(X_tilde) <= h), aes(color = treatment, size = weights)) +
    scale_size(range = c(1, 5), guide = "none") +
    geom_point(data = filter(data, abs(X_tilde) > h), aes(color = treatment), size = 2, alpha = 0.2) +
    scale_color_manual(labels = c("Control observations", "Treatment observations"),
                       values = c("control" = "blue", "treatment" = "red")) +
    stat_function(fun = control_fun_plot, color = "blue",
                  xlim = c(-cutoff, 0), linewidth = 1) +
    stat_function(fun = control_fun_plot, color = "blue",
                  xlim = c(0, cutoff), linewidth = 1, linetype = "dashed") +
    stat_function(fun = treatment_fun_plot, color = "red",
                  xlim = c(0, cutoff), linewidth = 1) +
    stat_function(fun = treatment_fun_plot, color = "red",
                  xlim = c(-cutoff, 0), linewidth = 1, linetype = "dashed") +
    annotate("rect", xmin = -h, xmax = 0, ymin = -Inf, ymax = Inf, fill = "blue", alpha = 0.1) +
    annotate("rect", xmin = 0, xmax = h, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
    geom_segment(data = coordinates_fits,
                 aes(x = x, y = y, xend = xend, yend = yend, linetype = "solid"),
                 linewidth = 1) +
    scale_linetype_manual(labels = "Local linear fit", values = c("solid" = "solid")) +
    geom_segment(x = 0, xend = 0, y = control_fun_plot(0), yend = treatment_fun_plot(0),
                 linetype = "dashed") +
    annotate("text", x = 0.03, y = 0.525, label = "tau[SRD]", parse = TRUE, size = 8) +
    annotate("text", x = 0.4, y = 0.5, label = "E[Y(0) | X]", color = "blue", size = 8) +
    annotate("text", x = 0.4, y = 0.75, label = "E[Y(1) | X]", color = "red", size = 8) +
    xlab("Assignment variable, X") + 
    scale_x_continuous(breaks = c(0, -h, h), labels = c("c", "c - h", "c + h")) +
    guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
    theme_light(base_size = 20) +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = c(0.13, 0.85), legend.title = element_blank(),
          legend.spacing.y = unit(-1, "mm"), legend.text = element_text(size = 20))

ggsave(filename = "../figures/figure_03.pdf", plot = ll_estimation_plot,
       width = 15, height = 8, units = "in", dpi = 600)