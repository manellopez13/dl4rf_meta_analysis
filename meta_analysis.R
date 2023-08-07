library(metafor)

# FUNCTION DEFINITION - - - - - -

preprocess_data <- function(df) {
  
  # Get lower and upper bounds intervals.
  df$diff_1 <- df$mean - df$bottom
  df$diff_2 <- df$top - df$mean
  
  # Assume: normal dist, and n big, so t student approx normal.
  # Diff 1 = 1.96 * SE.
  df$se <- (df$diff_1 + df$diff_2) / (2 * 1.96)
  
  # SE = SD/sqrt(n).
  df$sd <- df$se * sqrt(df$n)
  
  return(df)
  
}

# A little helper function to add Q-test, I^2, and tau^2 estimate info.
mlabfun <- function(text, res) {
  list(bquote(paste(.(text),
                    " (",
                    I^2, " = ", .(formatC(res$I2, digits=1, format="f")), "%, ",
                    tau^2, " = ", .(formatC(res$tau2, digits=2, format="f")), ")")))}

make_forest_plot <- function(df, xlabel, x_low, x_high) {
  
  # Reorder rows.
  df <- df[order(df$rob,df$study,decreasing=TRUE),]
  
  # Get effect size and random effects model results.
  es <- escalc(measure="MN", mi=mean, sdi=sd, ni=n, data=df, slab=study)
  res <- rma(yi=yi, vi=vi, data=es)
  
  n_pos = x_low + 20
  
  n_rows = nrow(df)
  n_rob = sum(df$rob == 1)
  n_low_rob = n_rows - n_rob
  
  row_1 = 3
  row_2 = row_1 + n_rob - 1
  row_3 = row_2 + 4
  row_4 = row_3 + n_low_rob - 1
  
  # Forest plot.
  forest(res, xlim=c(x_low, x_high), ylim=c(-2, n_rows + 9),
         header="Study",
         ilab=n, ilab.xpos=n_pos, # xpos is consistent with x axis
         mlab=mlabfun("RE model for all studies", res),
         rows=c(row_1:row_2, row_3:row_4),
         refline=res$b,
         xlab=xlabel,
  )
  
  op <- par(font=2) # bold
  text(n_pos, n_rows + 8, "n")
  
  ### switch to bold italic font
  par(font=4)
  
  ### add text for the subgroups
  text(x_low, c(row_4+1, row_2+1), pos=4, c("Low ROB", "ROB"))
  
  ### fit random-effects model in the three subgroups
  res.group0 <- rma(es[["yi"]], es[["vi"]], subset=(es$rob==0))
  res.group1 <- rma(es[["yi"]], es[["vi"]], subset=(es$rob==1))
  
  print(res.group0)
  
  ### add summary polygons for the three subgroups
  addpoly(res.group0, row=row_3-1.5, mlab=mlabfun("RE model for subgroup", res.group0))
  addpoly(res.group1, row=row_1-1.5, mlab=mlabfun("RE model for subgroup", res.group1))
  
  ### fit meta-regression model to test for subgroup differences
  resdiff <- rma(es[["yi"]], es[["vi"]], mods = ~ es$rob)
  
  ### add text for the test of subgroup differences
  text(x_low, -2, pos=4, cex=0.8, bquote(paste("Test for subgroup differences: ",
                                              Q[M], " = ", .(formatC(resdiff$QM, digits=2, format="f")), ", df = ", .(resdiff$p - 1),
                                              ", p = ", .(formatC(resdiff$QMp, digits=2, format="f")))))
  
}

# DATA DEFINITION - - - - - - - -

data <- read.csv("meta_analysis_data_with_sim_f1.csv")

sensitivity <- data[c("study", "n", "sens", "sens_b", "sens_t", "rob")]
colnames(sensitivity)[3] <- "mean"
colnames(sensitivity)[4] <- "bottom"
colnames(sensitivity)[5] <- "top"

# Removal of row 2, which has no report of precision.
precision <- data[c("study", "n", "prec", "prec_b", "prec_t", "rob")]
precision <- precision[-c(2),]
colnames(precision)[3] <- "mean"
colnames(precision)[4] <- "bottom"
colnames(precision)[5] <- "top"

# Removal of row 2, which has no report of precision nor F1-score.
f1 <- data[c("study", "n", "f1", "f1_b", "f1_t", "rob")]
f1 <- f1[-c(2),]
colnames(f1)[3] <- "mean"
colnames(f1)[4] <- "bottom"
colnames(f1)[5] <- "top"

# MAIN - - - - - - - -

sensitivity <- preprocess_data(sensitivity)
precision <- preprocess_data(precision)
f1 <- preprocess_data(f1)

make_forest_plot(sensitivity, "Sensitivity", -20, 140)
make_forest_plot(precision, "Precision", -10, 130)

# Add asterisks in F1-score plot to signal simulated data.
x_low = 25
x_high = 115
make_forest_plot(f1, "F1-score", x_low, x_high)
text(x_high-1.75, 4, pos=4, "*")
text(x_high-1.75, 18, pos=4, "*")
