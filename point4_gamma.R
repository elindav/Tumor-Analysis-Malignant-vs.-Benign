
# Sections 1-4, consistent with point4_gamma.Rmd.
# Run from the RStudio project directory.
# Descriptive R quantiles are used for outlier screening; lecture order-statistic quantiles are used for Q-Q plots.
if (!file.exists("data/wdbc_project.csv")) stop("Open the project directory containing data/wdbc_project.csv")
if (!dir.exists("outputs")) dir.create("outputs")
if (!dir.exists("figures")) dir.create("figures")

# ---- data-import ----
bc <- read.csv("data/wdbc_project.csv", stringsAsFactors = FALSE)

stopifnot(
  nrow(bc) == 569,
  ncol(bc) == 31,
  all(c("diagnosis", "area_mean") %in% names(bc)),
  all(bc$diagnosis %in% c("M", "B"))
)

numerical <- bc[, names(bc) != "diagnosis", drop = FALSE]
stopifnot(all(vapply(numerical, is.numeric, logical(1))))

missing_count <- sum(is.na(bc))
nonfinite_count <- sum(!is.finite(as.matrix(numerical)))
duplicate_count <- sum(duplicated(bc))
nonpositive_count <- sum(bc$area_mean <= 0, na.rm = TRUE)

if (missing_count > 0 || nonfinite_count > 0 || nonpositive_count > 0) {
  stop("Check data quality and positive Gamma support before analysis.")
}

x_M <- bc$area_mean[bc$diagnosis == "M"]
x_B <- bc$area_mean[bc$diagnosis == "B"]


# ---- descriptive-statistics ----
count_iqr_outliers <- function(x) {
  q <- as.numeric(quantile(x, c(0.25, 0.75)))
  iqr <- q[2] - q[1]
  sum(x < q[1] - 1.5 * iqr | x > q[2] + 1.5 * iqr)
}

area_summary <- data.frame(
  Diagnosis = c("Malignant", "Benign"),
  N = c(length(x_M), length(x_B)),
  Mean = c(mean(x_M), mean(x_B)),
  Median = c(median(x_M), median(x_B)),
  SD = c(sd(x_M), sd(x_B)),
  Minimum = c(min(x_M), min(x_B)),
  Maximum = c(max(x_M), max(x_B)),
  IQR_Outliers = c(count_iqr_outliers(x_M), count_iqr_outliers(x_B))
)

print(area_summary, row.names = FALSE, digits = 5)


# ---- comparative-boxplot ----
boxplot(
  x_B, x_M,
  names = c("Benign", "Malignant"),
  col = c("lightblue", "lightpink"),
  main = "Mean Nuclear Area by Diagnosis",
  ylab = "Mean nuclear area (dataset units)"
)

png("figures/01_comparative_boxplot.png", width=1400, height=700, res=130)
boxplot(x_B, x_M, names=c("Benign", "Malignant"), col=c("lightblue", "lightpink"), main="Mean Nuclear Area by Diagnosis", ylab="Mean nuclear area (dataset units)")
dev.off()

# ---- gamma-method-of-moments ----
fit_gamma_mom <- function(x) {
  m1 <- mean(x)
  m2 <- mean(x^2)
  v <- m2 - m1^2
  
  if (v <= 0) stop("Gamma fitting requires positive empirical variance.")
  
  list(
    shape = m1^2 / v,
    rate = m1 / v
  )
}

fit_M <- fit_gamma_mom(x_M)
fit_B <- fit_gamma_mom(x_B)

fits_table <- data.frame(
  Diagnosis = c("Malignant", "Benign"),
  Shape = c(fit_M$shape, fit_B$shape),
  Rate = c(fit_M$rate, fit_B$rate),
  Fitted_Mean = c(fit_M$shape / fit_M$rate,
                  fit_B$shape / fit_B$rate)
)

print(fits_table, row.names = FALSE, digits = 7)


# ---- gamma-histograms ----
plot_hist_gamma <- function(x, fit, label) {
  grid_x <- seq(min(x), max(x), length.out = 400)
  gamma_y <- dgamma(grid_x, shape = fit$shape, rate = fit$rate)
  h <- hist(x, breaks = "FD", freq = FALSE,
            right = FALSE, plot = FALSE)
  
  hist(x, breaks = h$breaks, freq = FALSE, right = FALSE,
       col = "grey85", border = "white",
       ylim = c(0, 1.1 * max(c(h$density, gamma_y))),
       main = label, xlab = "Mean nuclear area", ylab = "Density")
  lines(grid_x, gamma_y, col = "blue", lwd = 2)
}

par(mfrow = c(1, 2))
plot_hist_gamma(x_M, fit_M, "Malignant")
plot_hist_gamma(x_B, fit_B, "Benign")
par(mfrow = c(1, 1))

png("figures/02_gamma_histograms.png", width=1400, height=700, res=130)
par(mfrow=c(1,2)); plot_hist_gamma(x_M, fit_M, "Malignant"); plot_hist_gamma(x_B, fit_B, "Benign")
dev.off()

# ---- gamma-qq ----
course_sample_quantiles <- function(x, alpha) {
  ordered <- sort(x)
  indices <- pmax(1, floor(alpha * length(ordered)))
  ordered[indices]
}

plot_gamma_qq <- function(x, fit, label) {
  alpha <- seq(0.05, 0.95, by = 0.05)
  sample_q <- course_sample_quantiles(x, alpha)
  theoretical_q <- qgamma(alpha, shape = fit$shape,
                          rate = fit$rate)
  limits <- range(c(sample_q, theoretical_q))
  
  plot(theoretical_q, sample_q,
       xlim = limits, ylim = limits, asp = 1,
       pch = 19, cex = 0.8, col = "blue",
       main = label,
       xlab = "Theoretical Gamma quantiles",
       ylab = "Sample quantiles")
  abline(0, 1, col = "red", lwd = 2)
}

par(mfrow = c(1, 2))
plot_gamma_qq(x_M, fit_M, "Malignant")
plot_gamma_qq(x_B, fit_B, "Benign")
par(mfrow = c(1, 1))

png("figures/03_gamma_qq.png", width=1400, height=700, res=130)
par(mfrow=c(1,2)); plot_gamma_qq(x_M, fit_M, "Malignant"); plot_gamma_qq(x_B, fit_B, "Benign")
dev.off()

# ---- gamma-cdf ----
plot_gamma_cdf <- function(x, fit, label) {
  plot(ecdf(x), verticals = TRUE, do.points = FALSE,
       col = "grey40", main = label,
       xlab = "Mean nuclear area", ylab = "Cumulative probability")
  
  grid_x <- seq(min(x), max(x), length.out = 400)
  lines(grid_x,
        pgamma(grid_x, shape = fit$shape, rate = fit$rate),
        col = "blue", lwd = 2)
  legend("bottomright", c("Empirical CDF", "Fitted Gamma CDF"),
         col = c("grey40", "blue"), lwd = 2, bty = "n", cex = 0.8)
}

par(mfrow = c(1, 2))
plot_gamma_cdf(x_M, fit_M, "Malignant")
plot_gamma_cdf(x_B, fit_B, "Benign")
par(mfrow = c(1, 1))

png("figures/04_gamma_cdf.png", width=1400, height=700, res=130)
par(mfrow=c(1,2)); plot_gamma_cdf(x_M, fit_M, "Malignant"); plot_gamma_cdf(x_B, fit_B, "Benign")
dev.off()

# Saved summary tables for the group
write.csv(area_summary, "outputs/area_summary.csv", row.names = FALSE)
write.csv(fits_table, "outputs/gamma_mom_parameters.csv", row.names = FALSE)
cat("Updated Gamma MoM analysis complete. Figures saved in figures/.\n")
