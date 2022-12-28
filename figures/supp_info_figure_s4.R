library(ggplot2)
library(ggbiplot)

# load data

yearly.dat <- readRDS(file = "data/for_plotting/yearlydat.RDS")

# compare with PCA
ord_predictor_PCA <-
  princomp(x = yearly.dat[, -1], scores = T, cor = T)

# make plot
ggbiplot(ord_predictor_PCA,
         varname.size = 5) + 
  xlim(c(-3, 3))  +
  theme_bw() + 
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )
