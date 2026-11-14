# Simple functions used in the Gaussian simulation.
# All datasets have the same p predictors.

# Generate target data using the Gaussian glmtrans simulation model.
# n: sample size; K: total number of sources; p: predictors; s: nonzero coefficients.
# Returns a list containing x (predictor matrix) and y (response).
generate_target = function(n, K, p, s) 
{
  data = models(family = "gaussian", type = "target", cov.type = 2,
                n.target = n, K = K, s = s, p = p)$target
  return(data)
}

# Generate K source datasets, the first Ka of which are transferable.
# h: coefficient mismatch; n_source: observations per source.
# Returns a list of K source datasets with x and y.
generate_sources = function(h, K, Ka, n_source, p, s) 
{
  data = models(family = "gaussian", type = "source", cov.type = 2,
                h = h, K = K, Ka = Ka, n.source = rep(n_source, K),
                s = s, p = p)$source
  return(data)
}

# Mean squared prediction error.
mse = function(observed, predicted)
{
  return(mean((observed - predicted)^2))
}

# Mean absolute prediction error.
mae = function(observed, predicted)
{
  return(mean(abs(observed - predicted)))
}
