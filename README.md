# fisher_deriv.py runs CLASS many times for each parameter and fits a polynomial at each l to compute a polynomial derivative to be fed into a fisher matrix (all inside run_fisher function). It subsequently inverts to find the covariance matrix, extracts the forecasted constraints and computes the ellipse parameters.

# after the first run the CLASS outputs are imported from "data"
