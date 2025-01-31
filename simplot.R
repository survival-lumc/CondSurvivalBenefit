library(tidyverse)
library(ggforce)
library(data.table)
#library(ggExtra) 
library(cowplot) 
library(ggpubr)

summ_cs <- readRDS("Simulation_livers/resultsallfinal.rds")

# Put summary data in longitudinal format --------------------------------------
# ------------------------------------------------------------------------------
summ_cs$ben_naive <- summ_cs$esrmst1_naive - summ_cs$esrmst0_naive
summ_cs$ben_cs <- summ_cs$esrmst1_cs - summ_cs$esrmst0_cs
summ_cs$ben_w <- summ_cs$esrmst1_w - summ_cs$esrmst0_w

l_df <- data.frame(
  nrun = rep(summ_cs$nrun, each = 9),
  cs = rep(summ_cs$cs, each = 9),
  rmst0 = rep(summ_cs$rmst0, each = 9),
  rmst1 = rep(summ_cs$rmst1, each = 9),
  Estimand = rep(c(rep("esrmst0", 3), rep("esrmst1", 3), rep("ben", 3)), nrow(summ_cs)),
  Method = rep(c("naive", "cs", "w"), 3*nrow(summ_cs)),
  Estimate = c(t(matrix(c(summ_cs$esrmst0_naive, summ_cs$esrmst0_cs, summ_cs$esrmst0_w,
                          summ_cs$esrmst1_naive, summ_cs$esrmst1_cs, summ_cs$esrmst1_w,
                          summ_cs$ben_naive, summ_cs$ben_cs, summ_cs$ben_w), ncol = 9)))
)

rm(summ_cs)

l_df$Truth <- case_when(
  l_df$Estimand =="esrmst0" ~  l_df$rmst0,
  l_df$Estimand =="esrmst1" ~  l_df$rmst1,
  l_df$Estimand =="ben" ~ (l_df$rmst1 - l_df$rmst0)
)

l_df$Estimand <- factor(l_df$Estimand, levels = c("esrmst0", "esrmst1", "ben"),
                        labels = c("RMST0", "RMST1", "Benefit"))
l_df$Method <- factor(l_df$Method, levels = c("naive", "cs", "w"),
                      labels = c("Naive", "Cross-sections unweighted", "Cross-sections weighted"))

l_df <- select(l_df, -c("rmst0", "rmst1", "nrun"))


### Loess ----------------------------------------------------------------------
### ----------------------------------------------------------------------------
my_gam <- function (method, estimand) {
  gam_fit <- mgcv::gam(Truth ~ s(Estimate, bs = "cs"), data = l_df[l_df$Method == method & l_df$Estimand == estimand,])
  Truth_cal <- predict(gam_fit, l_df[l_df$Method == method & l_df$Estimand == estimand,])
  return(Truth_cal)
}
# Root mean square bias
l_df$Truth_cal <- NA_real_
l_df$Truth_cal[l_df$Method == "Naive" & l_df$Estimand == "RMST0"] <- my_gam("Naive", "RMST0")
l_df$Truth_cal[l_df$Method == "Cross-sections unweighted" & l_df$Estimand == "RMST0"] <- my_gam("Cross-sections unweighted", "RMST0")
l_df$Truth_cal[l_df$Method == "Cross-sections weighted" & l_df$Estimand == "RMST0"] <- my_gam("Cross-sections weighted", "RMST0")
l_df$Truth_cal[l_df$Method == "Naive" & l_df$Estimand == "RMST1"] <- my_gam("Naive", "RMST1")
l_df$Truth_cal[l_df$Method == "Cross-sections unweighted" & l_df$Estimand == "RMST1"] <- my_gam("Cross-sections unweighted", "RMST1")
l_df$Truth_cal[l_df$Method == "Cross-sections weighted" & l_df$Estimand == "RMST1"] <- my_gam("Cross-sections weighted", "RMST1")
l_df$Truth_cal[l_df$Method == "Naive" & l_df$Estimand == "Benefit"] <- my_gam("Naive", "Benefit")
l_df$Truth_cal[l_df$Method == "Cross-sections unweighted" & l_df$Estimand == "Benefit"] <- my_gam("Cross-sections unweighted", "Benefit")
l_df$Truth_cal[l_df$Method == "Cross-sections weighted" & l_df$Estimand == "Benefit"] <- my_gam("Cross-sections weighted", "Benefit")


### Plot function --------------------------------------------------------------
### ----------------------------------------------------------------------------
cols <- c("Naive" = "#D41159", "Cross-sections unweighted" = "#FFC20A", 
          "Cross-sections weighted" = "#0C7BDC", "Perfect calibration" = "darkgrey")
lline <- c("Naive" = "longdash", "Cross-sections unweighted" = "dashed", 
           "Cross-sections weighted" = "solid", "Perfect calibration" = "solid")


myplot <- function(data, xlims = c(0,3)) {
  
  df <- data
  x1 <- round(xlims[[1]])
  x2 <- round(xlims[[2]])
  len <- x2-x1
  
  rmsb <- df |>
    group_by(Method) |>
    summarise(RMSB = sqrt(sum((Truth_cal - Estimate)^2)/n()))
  
  p <- ggplot(df, aes(x = Estimate, y = Truth_cal, color = Method, fill = Method, linetype = Method)) +
    geom_line(linewidth = 1) +
    #geom_smooth(alpha = 0.3, linewidth = 1.25,
    #            method = 'gam', formula = formula(y ~ s(x, bs = "cs"))) +
    geom_abline(slope = 1, color = "grey") +
    scale_x_continuous(breaks = seq(x1, x2, by = len/6))+
    scale_y_continuous(breaks = seq(x1, x2, by = len/6)) +
    scale_linetype_manual(
      name = "Method",
      values = lline
    ) +
    scale_colour_manual(
      name = "Method",
      values = cols,
      aesthetics = c("colour", "fill")
    ) +
    ylab("Truth") +
    theme_bw() + 
    theme(legend.position="none") +
    coord_cartesian(xlim = xlims) 
  
  xmin <- ggplot_build(p)$layout$panel_params[[1]]$x$continuous_range[[1]]
  xmax <- ggplot_build(p)$layout$panel_params[[1]]$x$continuous_range[[2]]
  rangex <- xmax - xmin
  ymin <- ggplot_build(p)$layout$panel_params[[1]]$y$continuous_range[[1]]
  ymax <- ggplot_build(p)$layout$panel_params[[1]]$y$continuous_range[[2]]
  rangey <- ymax - ymin
  
  p <- p + 
    annotate(
      geom = "rect", 
      xmin = xmin + 0.1 * rangex, xmax = xmin + 0.5 * rangex, 
      ymin = ymin + 0.7 * rangey, ymax = ymin + 0.95 * rangey,
      alpha = 1, 
      fill = "white", 
      color = "darkgrey", 
      linewidth = 1
    ) +
    annotate(
      geom = "text", 
      label="Bias", 
      x = xmin + 0.115 * rangex,
      y = ymin + 0.9 * rangey, 
      hjust = 0, 
      col = "black", 
      size=4,
      fontface = 2
    ) +
    annotate(
      geom = "text", 
      label=paste("Naive:", formatC(round(rmsb$RMSB[rmsb$Method == "Naive"],3),3,format="f")),
      x = xmin + 0.115 * rangex,
      y = ymin + 0.83 * rangey, 
      hjust = 0, 
      col = cols["Naive"], 
      size=3.5
    ) +
    annotate(
      geom = "text", 
      label=paste("CS unweighted:", 
                  formatC(round(rmsb$RMSB[rmsb$Method == "Cross-sections unweighted"],3),3,format="f")), 
      x = xmin + 0.115 * rangex,
      y = ymin + 0.785 * rangey, 
      hjust = 0, 
      col = cols["Cross-sections unweighted"], 
      size=3.5
    ) +
    annotate(
      geom = "text", 
      label=paste("CS weighted:", 
                  formatC(round(rmsb$RMSB[rmsb$Method == "Cross-sections weighted"],3),3,format="f")), 
      x = xmin + 0.115 * rangex,
      y = ymin + 0.740 * rangey, 
      hjust = 0, 
      col = cols["Cross-sections weighted"], 
      size=3.5
    ) 
  
  ph <- ggplot(df, aes(x = Truth)) +
    geom_histogram(alpha = 0.9, colour = "black", fill = "lightblue", binwidth = 0.05) +
    theme_void() + 
    xlim(xlims) +
    coord_flip()
  
  plot_grid(
    p, 
    plot_grid(ph, NULL, ncol =1, nrow=2, rel_heights = c(0.9, 0.1)), 
    rel_widths = c(0.85, 0.15)
  )
  
}

### Plot -----------------------------------------------------------------------
### ----------------------------------------------------------------------------
p0 <- myplot(l_df[l_df$Estimand == "RMST0", c("Truth","Truth_cal", "Estimate", "Method")])
p1 <- myplot(l_df[l_df$Estimand == "RMST1", c("Truth","Truth_cal", "Estimate", "Method")])
pb <- myplot(l_df[l_df$Estimand == "Benefit", c("Truth","Truth_cal", "Estimate", "Method")], c(-2.85, 3))
mylegend <- get_legend(
  ggplot(l_df, aes(x = Estimate, y = Truth, colour = Method, linetype = Method)) + 
    geom_line(linewidth = 1.25) +
    scale_linetype_manual(
      name = "Method",
      values = lline
    ) +
    scale_colour_manual(
      name = "Method",
      values = cols,
      aesthetics = c("colour", "fill")
    )
)
mygglegend <- as_ggplot(mylegend)


pg <- plot_grid(
  plot_grid(NULL),
  plot_grid(p0 , p1, nrow = 1, ncol = 2, labels = c("RMST0", "RMST1"), 
            vjust = 0, hjust = -3),
  plot_grid(NULL, pb, mygglegend, nrow = 1, rel_widths = c(0.25, 0.5, 0.25),
            labels = c("", "Benefit", ""), vjust = 0, hjust = -3),
  nrow = 3,
  rel_heights = c(0.05, 0.475, 0.475)
)









