
#20bases functions
load(file="../../06INLA/Poisson/Basis_C_phi3_Gen25k_1_INLA.RData")
ySeq_INLA20<-dnorm(x = xSeq , mean = linear_predictor_mean[(0.8*n)+ktheta]  , sd = linear_predictor_sd[(0.8*n)+ktheta] )

#100bases functions
load(file="../../06INLA/Poisson/100Basis_C_phi3_Gen25k_1_INLA.RData")
# INLA
ySeq_INLA100<-dnorm(x = xSeq , mean = linear_predictor_mean[(0.8*n)+ktheta]  , sd = linear_predictor_sd[(0.8*n)+ktheta] )

##############################################################

#50bases functions
##############################################################
plot(x = xSeq, y = ySeq_INLA, type = "l", col = "blue", 
     xlim = c(-1.5, 0.5), ylim = c(0, 5))

lines(x = xSeq, y = ySeq_INLA100, col = "orange")
lines(x = xSeq, y = ySeq_INLA20, col = "red")
abline(v = eta_testdata_real[ktheta], col = "green", lty = 2, lwd = 2)

legend("topright", 
       legend = c("INLA (20 basis)", "INLA (50 basis)", "INLA (100 basis)", "True Eta"), 
       col = c("red","blue", "orange", "green"),
       lty = c(1, 1, 1, 2), 
       lwd = c(2, 2, 2, 2)) 

##############################################################
##############################################################
##############################################################
library(ggplot2)

# Create a data frame for the main lines
df <- data.frame(
  x = rep(xSeq, 3),
  y = c(ySeq_INLA100, ySeq_INLA, ySeq_INLA20),
  Method = factor(rep(c("INLA (100 basis)", "INLA (50 basis)", "INLA (20 basis)"), 
                      each = length(xSeq)),
                  levels = c("INLA (100 basis)", "INLA (50 basis)", "INLA (20 basis)"))
)

# Create a data frame for the vertical line (True η)
vline_df <- data.frame(
  x = eta_testdata_real[ktheta],
  Method = "True η"
)

# Combine levels so legend includes "True η"
all_levels <- c("INLA (100 basis)", "INLA (50 basis)", "INLA (20 basis)", "True η")

ggplot() +
  geom_line(data = df, aes(x = x, y = y, color = Method, linetype = Method), size = 1.2) +
  geom_vline(data = vline_df, aes(xintercept = x, color = Method, linetype = Method), size = 1) +
  scale_color_manual(
    values = c("INLA (100 basis)" = "red", 
               "INLA (50 basis)" = "blue", 
               "INLA (20 basis)" = "green", 
               "True η" = "black"),
    breaks = all_levels
  ) +
  scale_linetype_manual(
    values = c("INLA (100 basis)" = "solid", 
               "INLA (50 basis)" = "solid", 
               "INLA (20 basis)" = "solid", 
               "True η" = "dashed"),
    breaks = all_levels
  ) +
  labs(
    x = expression(eta[i]),
    y = "Density",
    color = "Method",
    linetype = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_cartesian(xlim = c(-1.5, 0.5), ylim = c(0, 5))

