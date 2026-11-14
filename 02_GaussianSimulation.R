# ===========================================
# Trans-GLM Gaussian simulation
# 6 sources, 0-6 transferable, h = 1, 10, 20
# 10 repetitions for each setting
# ===========================================

# Packages and functions
source("00_requirements.R")
source("01_SimulationFunctions.R")


# Variables
reps = 10
h_values = c(1, 10, 20)
K = 6                 # Total source datasets
p = 500               # Number of predictors
s = 10                # Number of nonzero coefficients
n_target = 100
n_source = 100         # Observations in each source
n_test = 1000

# Results
results = data.frame()
detection = data.frame()

# Simulation
for (r in 1:reps) 
{
  set.seed(1000 + r)
  target = generate_target(n_target, K, p, s)

  test = generate_target(n_test, K, p, s)

  # Method 1: Fit target-only Lasso once
  set.seed(2000 + r)
  fit_target = cv.glmnet(target$x, target$y, family = "gaussian", alpha = 1, nfolds = 3)
  pred_target = as.numeric(predict(fit_target, newx = test$x, s = "lambda.min"))

  for (h in h_values) 
  {
    for (Ka in 0:K) 
    {
      set.seed(3000 + 100 * r + Ka)
      source_data = generate_sources(h, K, Ka, n_source, p, s)

      # Method 2: Trans-GLM using every source.
      set.seed(4000 + 100 * r + Ka)
      fit_all = glmtrans(target = target, source = source_data,
                         family = "gaussian", transfer.source.id = "all",
                         alpha = 1, nfolds = 3, valid.nfolds = 3,
                         detection.info = FALSE)
      pred_all = as.numeric(predict(fit_all, newx = test$x))

      # Method 3: Trans-GLM automatically choosing sources.
      set.seed(5000 + 100 * r + Ka)
      fit_auto = glmtrans(target = target, source = source_data,
                          family = "gaussian", transfer.source.id = "auto",
                          alpha = 1, nfolds = 3, valid.nfolds = 3,
                          detection.info = FALSE)

      pred_auto = as.numeric(predict(fit_auto, newx = test$x))
      selected = fit_auto$transfer.source.id

      # Compare predictions with the independent test outcomes.
      results = rbind(results,
        data.frame(h = h, Ka = Ka, repetition = r,
                   method = "Target-only Lasso",
                   mse = mse(test$y, pred_target),
                   mae = mae(test$y, pred_target)),
        data.frame(h = h, Ka = Ka, repetition = r,
                   method = "Trans-GLM: all",
                   mse = mse(test$y, pred_all),
                   mae = mae(test$y, pred_all)),
        data.frame(h = h, Ka = Ka, repetition = r,
                   method = "Trans-GLM: auto",
                   mse = mse(test$y, pred_auto),
                   mae = mae(test$y, pred_auto)))

      # Record which sources automatic Trans-GLM chose.
      for (k in 1:K) 
      {
        detection = rbind(detection, data.frame(
          h = h, Ka = Ka, repetition = r, source = k,
          transferable = (k <= Ka), selected = (k %in% selected)))
      }

      cat("Repetition", r, "of", reps,
          "| h =", h, "| Ka =", Ka, "\n")
    }
  }
}

# Results
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

write.csv(results, "results/simulation_results.csv", row.names = FALSE)
write.csv(detection, "results/source_detection.csv", row.names = FALSE)

# Average MSE and MAE
summary_results = aggregate(cbind(mse, mae) ~ h + Ka + method, data = results, FUN = mean)
write.csv(summary_results, "results/average_results.csv", row.names = FALSE)
print(summary_results)

# Frequency each source was selected
selection_rates = aggregate(selected ~ h + Ka + transferable, data = detection, FUN = mean)
write.csv(selection_rates, "results/selection_rates.csv", row.names = FALSE)

# Plots
# Figure 1: MSE by number of transferable sources.
plot_mse = ggplot(summary_results, aes(x = Ka, y = mse, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ h, labeller = label_both) +
  scale_x_continuous(breaks = 0:6) +
  labs(
    title = "Prediction error by number of transferable sources",
    x = "Number of transferable sources (out of 6)",
    y = "Average test MSE", color = "Method") + theme_minimal()

ggsave("figures/mse_lineplot.png", plot_mse, width = 12, height = 5)

# Figure 2: How often each source was selected
plot_detection = ggplot(selection_rates,
  aes(x = Ka, y = selected,
      color = factor(transferable, levels = c(FALSE, TRUE),
      labels = c("Not transferable", "Transferable")),
      group = transferable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ h, labeller = label_both) +
  scale_x_continuous(breaks = 0:6) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "How often automatic Trans-GLM selects sources",
    x = "Number of transferable sources (out of 6)",
    y = "Average selection rate", color = "Source type") +
  theme_minimal()

ggsave("figures/source_selection_lineplot.png", plot_detection, width = 12, height = 5)
