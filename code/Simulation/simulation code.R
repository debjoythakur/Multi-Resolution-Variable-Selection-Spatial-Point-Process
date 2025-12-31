## ==========================================================
## OVERVIEW (what this script does)
## ==========================================================
## 1) Create a spatial grid on [0,1]^2 and simulate:
##    - spatial covariate fields X_p(s) from a Gaussian process (GP)
##    - spatially varying true coefficient surfaces beta_p(s) (piecewise patterns)
## 2) Define the *true* Poisson intensity on the grid:
##       log lambda(s) = b0 + contrast * z(s)
##    where z(s) is a standardized linear combination of the *active* covariates.
##    b0 is calibrated so that E[#points] ≈ mu_target.
## 3) Simulate a point pattern (events) from the intensity:
##    - IPPP: inhomogeneous Poisson process
##    - Thomas: clustered process = baseline intensity × cluster field S(s)
## 4) Convert point pattern likelihood into a Poisson regression problem
##    using quadrature (data + dummy points) with offset log(w_i).
## 5) Fit five methods (variable selection):
##    - Local LASSO  (Haar-localized features) with BIC selection along glmnet path
##    - Local SCAD   (Haar-localized features) with BIC selection along ncvreg path
##    - Global LASSO (plain covariates) with BIC selection along glmnet path
##    - Global SCAD  (plain covariates) with BIC selection along ncvreg path
##    - Adaptive LASSO (plain covariates) using ridge-init + weighted LASSO, both via BIC
## 6) Evaluate performance:
##    - beta_rmspe: RMS error between estimated beta surfaces and true beta surfaces
##    - TPR_global / FPR_global: recovery of active covariates {1,...,p0}
##    - local_TPR: location-wise set-based recovery of nonzero covariates
##    - runtime_fit: time for fitting + selection
## 7) Repeat over mu_target (large-sample behavior) and/or clustering regimes.

rm(list = ls(all = TRUE))
suppressPackageStartupMessages({
  library(spatstat)
  library(glmnet)   # LASSO / ridge / adaptive LASSO
  library(ncvreg)   # SCAD
})

## ---------- Helpers ----------
normalize_to_unit_square <- function(s) {
  s <- as.matrix(s)
  r1 <- range(s[,1]); r2 <- range(s[,2])
  cbind((s[,1]-r1[1])/max(r1[2]-r1[1], 1e-12),
        (s[,2]-r2[1])/max(r2[2]-r2[1], 1e-12))
}

haar_phi_1d <- function(j,k,x){ a <- k/2^j; b <- (k+1)/2^j; as.numeric(x>=a & x<b)*2^(j/2) }
haar_psi_1d <- function(j,k,x){
  a <- k/2^j; m <- (k+0.5)/2^j; b <- (k+1)/2^j
  (as.numeric(x>=a & x<m)-as.numeric(x>=m & x<b))*2^(j/2)
}
build_haar2d_meta <- function(s01, j0=0, J=1){
  x <- s01[,1]; y <- s01[,2]
  Phi <- list(); meta <- list(); cnt <- 0L
  n_k <- max(2^j0, 1)
  for (kx in 0:(n_k-1)) for (ky in 0:(n_k-1)) {
    cnt <- cnt+1L; Phi[[cnt]] <- haar_phi_1d(j0,kx,x)*haar_phi_1d(j0,ky,y)
    meta[[cnt]] <- list(j=j0, type="phi", kx=kx, ky=ky)
  }
  if (J > j0) {
    for (j in j0:(J-1)) for (kx in 0:(2^j-1)) for (ky in 0:(2^j-1)) {
      cnt <- cnt+1L; Phi[[cnt]] <- haar_psi_1d(j,kx,x)*haar_phi_1d(j,ky,y)
      meta[[cnt]] <- list(j=j, type="psi_x", kx=kx, ky=ky)
      cnt <- cnt+1L; Phi[[cnt]] <- haar_phi_1d(j,kx,x)*haar_psi_1d(j,ky,y)
      meta[[cnt]] <- list(j=j, type="psi_y", kx=kx, ky=ky)
      cnt <- cnt+1L; Phi[[cnt]] <- haar_psi_1d(j,kx,x)*haar_psi_1d(j,ky,y)
      meta[[cnt]] <- list(j=j, type="psi_xy", kx=kx, ky=ky)
    }
  }
  B <- do.call(cbind, Phi)
  meta_df <- do.call(rbind, lapply(meta, as.data.frame))
  list(B=B, meta=meta_df)
}
build_localized_features_safe <- function(X, basis){
  X <- as.matrix(X); basis <- as.matrix(basis)
  stopifnot(nrow(X) == nrow(basis))
  n <- nrow(X); P <- ncol(X); R <- ncol(basis)
  Z <- matrix(0.0, n, P*R)
  for (p in seq_len(P)) {
    cols <- ((p-1)*R + 1):(p*R)
    Z[, cols] <- basis * X[, p]
  }
  Z
}
simulate_gp_fields <- function(coords, P, sigma2=1, rho=0.25) {
  n <- nrow(coords)
  Dx <- as.matrix(dist(coords))
  K  <- sigma2 * exp(-Dx / rho)
  diag(K) <- diag(K) + 1e-8
  L  <- chol(K)
  Z0 <- matrix(rnorm(n*P), n, P)
  L %*% Z0
}
beta_truth_mat <- function(coords01, pn = 40, p0 = 5) {
  stopifnot(is.matrix(coords01) && ncol(coords01) == 2)
  x <- coords01[, 1]; y <- coords01[, 2]
  n <- nrow(coords01)
  B <- matrix(0, n, pn)
  b1 <- numeric(n); b1[(x < .5)&(y < .5)] <- 1; b1[(x >= .5)&(y >= .5)] <- -1; B[,1] <- b1
  b2 <- rep(-1,n); b2[(x < .5)&(y < .5)] <- 1; b2[(x >= .5)&(x < .75)&(y >= .5)&(y < .75)] <- .5
  b2[(x >= .75)&(y >= .75)] <- 0; B[,2] <- b2
  B[,3] <- ifelse(x <= .5, 1, 0)
  B[,4] <- ifelse(y <= .5, 1, 0)
  b5 <- numeric(n)
  b5[(x <= .5)&(y >  .5)] <- 1;       b5[(x > .5)&(y >  .5)] <- 0
  b5[(x <= .5)&(y <= .5)] <- sqrt(2); b5[(x > .5)&(y <= .5)] <- 0
  B[,5] <- b5
  B[,6] <- ifelse(x + y <= .5, sqrt(2), 0)
  B[,7] <- ifelse(x - y <= .5, sqrt(2), 0)
  B[,8] <- ifelse(x + y <= 1.5, sqrt(2), 0)
  B[,9] <- ifelse(x - y >= -.5, sqrt(2), 0)
  b10 <- numeric(n)
  in_band <- (x + y >= .5) & (x + y <= 1.5)
  left   <- (x - y <= -.5) & in_band
  right  <- (x - y >=  .5) & in_band
  center <- in_band & (abs(x - y) <= .5)
  bottom <- (x + y <= .5) & (abs(x - y) <= .5)
  top    <- (x + y >= 1.5) & (abs(x - y) <= .5)
  b10[left]  <- 1; b10[right] <- 1; b10[center] <- .5; b10[bottom] <- sqrt(2); b10[top] <- -1
  B[,10] <- b10
  B
}
rmse <- function(a,b) sqrt(mean((a-b)^2))



## ---------- Common simulation + quadrature prep (shared by all methods) ----------
prep_ippp <- function(seed, pn, p0, qtiles, J, j0, mu_target, gp_sigma2, gp_rho, contrast_s){
  set.seed(seed)
  W  <- owin(c(0,1), c(0,1))
  gx <- qtiles[1]; gy <- qtiles[2]
  xs <- seq(0.5/gx, 1-0.5/gx, length.out = gx)
  ys <- seq(0.5/gy, 1-0.5/gy, length.out = gy)
  grid  <- as.matrix(expand.grid(x=xs, y=ys))
  grid01 <- normalize_to_unit_square(grid)
  
  # Truth + covariates on grid
  Xgrid1 <- simulate_gp_fields(grid, P = p0, sigma2 = gp_sigma2, rho = gp_rho)
  Xgrid2 <- matrix(0, nrow = nrow(Xgrid1), ncol = pn - p0)
  Xgrid  <- cbind(Xgrid1, Xgrid2); colnames(Xgrid) <- paste0("x", seq_len(ncol(Xgrid)))
  Beta   <- beta_truth_mat(grid01, pn = pn, p0 = p0)  # n x pn
  
  base_linear <- rowSums(Beta[,1:p0, drop=FALSE] * Xgrid[,1:p0, drop=FALSE])
  z   <- scale(base_linear)[,1]
  b0  <- 0
  tmp <- exp(b0 + contrast_s * z)
  b0  <- b0 + log(mu_target / mean(tmp))
  lambda_true_grid <- exp(b0 + contrast_s * z)
  
  # IPPP simulation, quadrature (data+dummy) with offsets
  lambda_im <- im(xcol = xs, yrow = ys,
                  mat = matrix(lambda_true_grid, nrow=gx, ncol=gy, byrow=FALSE))
  Xpp <- rpoispp(lambda = lambda_im, win = W)
  
  Q  <- quadscheme(Xpp, method="grid", ntile=qtiles)
  isD_int <- as.logical(is.data(Q))
  wi <- try(as.numeric(weights(Q)), silent=TRUE)
  if (inherits(wi,"try-error") || length(wi)==0L) wi <- if (!is.null(Q$w)) as.numeric(Q$w) else numeric(0)
  if (length(wi)==0L && exists("w.quad", where=asNamespace("spatstat.model"), inherits=FALSE))
    wi <- as.numeric(getFromNamespace("w.quad","spatstat.model")(Q))
  Upp <- try(spatstat.geom::as.ppp(Q), silent=TRUE)
  if (!inherits(Upp,"try-error")) {
    U  <- data.frame(x=Upp$x, y=Upp$y)
    yi <- as.integer(isD_int)
  } else {
    Xd  <- if (!is.null(Q$data)) Q$data else spatstat.model::data.quad(Q)
    Xdm <- if (!is.null(Q$dummy)) Q$dummy else spatstat.model::dummy.quad(Q)
    U   <- data.frame(x=c(Xd$x, Xdm$x), y=c(Xd$y, Xdm$y))
    yi  <- c(rep(1L, npoints(Xd)), rep(0L, npoints(Xdm)))
  }
  M <- nrow(U)
  Wow <- as.owin(W)
  cell_area <- (diff(Wow$xrange)/gx) * (diff(Wow$yrange)/gy)
  if (length(wi)!=M) wi <- rep(NA_real_, M)
  bad <- !is.finite(wi) | wi<=0
  if (all(bad)) wi <- ifelse(yi==1L, 0, cell_area) else wi[bad] <- cell_area
  off <- log(pmax(wi, 1e-12))
  
  # Map quadrature points to nearest grid cell and pull covariates
  ix <- pmin(pmax(1, round(U$x * gx + 0.5)), gx)
  iy <- pmin(pmax(1, round(U$y * gy + 0.5)), gy)
  Xq <- Xgrid[(iy - 1) * gx + ix, , drop=FALSE]
  
  # Haar bases (for localized methods)
  U01 <- normalize_to_unit_square(U)
  B_qu_U   <- build_haar2d_meta(U01, j0=j0, J=J)
  B_U      <- B_qu_U$B
  B_qu_G   <- build_haar2d_meta(grid01, j0=j0, J=J)
  B_grid   <- B_qu_G$B
  R        <- ncol(B_U)
  
  list(
    W=W, gx=gx, gy=gy,
    grid=grid, grid01=grid01,
    Xgrid=Xgrid, Beta=Beta,
    lambda_true_grid=lambda_true_grid,
    U=U, yi=yi, off=off,
    ix=ix, iy=iy, Xq=Xq,
    B_U=B_U, B_grid=B_grid, R=R
  )
}
```


## Local Metrices
```{r}
## ---------- Local metrics helper (atom-level, existing) ----------
compute_local_metrics <- function(B_U, B_grid, Beta, Gamma, U01, pn, p0){
  support_hat_atom <- abs(Gamma) > 1e-8
  sel_cov_hat <- colSums(support_hat_atom) > 0
  S_true <- 1:p0
  S_hat  <- which(sel_cov_hat)
  TP_ids <- intersect(S_hat, S_true)
  FP_ids <- setdiff(S_hat, S_true)
  TPR_cov <- length(TP_ids) / p0
  FPR_cov <- length(FP_ids) / (ncol(Beta) - p0)
  
  # Local atom-level
  Beta_U <- beta_truth_mat(U01, pn = ncol(Beta), p0 = p0)
  BtB   <- crossprod(B_U)
  BtY   <- crossprod(B_U, Beta_U)
  Theta_true <- solve(BtB, BtY)
  support_true_atom <- abs(Theta_true) > 1e-8
  
  active_atoms_per_i <- lapply(seq_len(nrow(B_U)), function(i) which(B_U[i, ] != 0))
  is_true_nz_loc <- abs(Beta_U) > 1e-12
  
  local_results <- data.frame(p = integer(0), i = integer(0), TPR = numeric(0), FPR = numeric(0))
  for (i in seq_len(nrow(B_U))) {
    Ai <- active_atoms_per_i[[i]]
    if (!length(Ai)) next
    for (p in S_true) {
      if (!is_true_nz_loc[i, p]) next
      true_set <- which(support_true_atom[Ai, p])
      hat_set  <- which(support_hat_atom [Ai, p])
      
      TPc <- length(intersect(true_set, hat_set))
      FNc <- length(setdiff(true_set, hat_set))
      FPc <- length(setdiff(hat_set, true_set))
      TNc <- length(setdiff(seq_along(Ai), union(true_set, hat_set)))
      
      TPR_loc <- if ((TPc + FNc) > 0) TPc / (TPc + FNc) else NA_real_
      FPR_loc <- if ((FPc + TNc) > 0) FPc / (FPc + TNc) else NA_real_
      local_results <- rbind(local_results, data.frame(p = p, i = i, TPR = TPR_loc, FPR = FPR_loc))
    }
  }
  
  list(TPR_cov=TPR_cov, FPR_cov=FPR_cov)
}

## ---------- NEW helper: set-based local TPR (your definition) ----------
##  J_true(i) := { p : beta_p(s_i) != 0 }
##  J_hat(i)  := { p : exists active atom r at i with Gamma[r,p] != 0 }
compute_local_TPR_sets <- function(B_U, U01, Gamma, Beta, tol = 1e-8) {
  stopifnot(nrow(B_U) == nrow(U01), ncol(Gamma) == ncol(Beta))
  nU  <- nrow(B_U)
  pn  <- ncol(Beta)
  p0_infer <- sum(colSums(abs(Beta) > tol) > 0)
  Beta_U   <- beta_truth_mat(U01, pn = pn, p0 = p0_infer)
  active_atoms_per_i <- lapply(seq_len(nU), function(i) which(abs(B_U[i, ]) != 0))
  ratios <- vector("numeric", nU); ratios[] <- NA_real_
  for (i in seq_len(nU)) {
    Ai <- active_atoms_per_i[[i]]
    J_true_i <- which(abs(Beta_U[i, ]) > tol)
    if (length(J_true_i) == 0L) next
    if (length(Ai)) {
      any_nonzero_per_p <- colSums(abs(Gamma[Ai, , drop = FALSE]) > tol) > 0
      J_hat_i <- which(any_nonzero_per_p)
    } else {
      J_hat_i <- integer(0)
    }
    inter_sz <- length(intersect(J_true_i, J_hat_i))
    ratios[i] <- inter_sz / length(J_true_i)
  }
  list(
    local_TPR = mean(ratios, na.rm = TRUE),
    per_location  = ratios,
    n_effective   = sum(!is.na(ratios))
  )
}

bic_poisson_offset <- function(y, X, beta, offset, n_eff = length(y), weights = NULL) {
  eta <- as.numeric(X %*% beta + offset)
  ll  <- loglik_poisson_offset(y, eta, weights)
  s   <- sum(abs(beta) > 0)
  -2*ll + s * log(n_eff)
}
loglik_poisson_offset <- function(y, eta, weights = NULL) {
  mu <- exp(eta)
  if (is.null(weights)) sum(y * eta - mu) else sum(weights * (y * eta - mu))
}
## ========================= Robust path selectors (NO STANDARDIZATION) =========================
select_by_bic_glmnet <- function(X, y, offset,
                                 penalty.factor = NULL,
                                 family = "poisson") {
  n_eff = length(y)
  X <- as.matrix(X)
  pf <- if (is.null(penalty.factor)) rep(1, ncol(X)) else as.numeric(penalty.factor)
  fit <- glmnet::glmnet(
    x = X, y = y, family = family, alpha = 1,
    standardize = FALSE,                # <-- NO STANDARDIZATION
    offset = offset,
    penalty.factor = pf
  )
  lambdas <- fit$lambda
  coefs   <- stats::coef(fit)
  BICs <- numeric(length(lambdas))
  betas <- vector("list", length(lambdas))
  for (k in seq_along(lambdas)) {
    beta_full <- as.numeric(coefs[, k])   # (Intercept + p)
    beta <- beta_full[-1]
    BICs[k] <- bic_poisson_offset(y, X, beta, offset, n_eff = n_eff)
    betas[[k]] <- beta
  }
  kmin <- which.min(BICs)
  list(beta = betas[[kmin]], lambda = lambdas[kmin],
       bic_min = BICs[kmin],
       bic_path = data.frame(lambda = lambdas, bic = BICs))
}

select_by_bic_ncvreg <- function(X, y, offset ) {
  n_eff = length(y)
  fit <- ncvreg::ncvreg(X = X, y = y, family = "poisson", penalty = "SCAD", offset = offset)
  lambdas <- fit$lambda
  BICs <- numeric(length(lambdas))
  betas <- vector("list", length(lambdas))
  for (k in seq_along(lambdas)) {
    beta_full <- as.numeric(fit$beta[, k]) # (Intercept + p)
    beta <- beta_full[-1]
    BICs[k] <- bic_poisson_offset(y, X, beta, offset, n_eff = n_eff)
    betas[[k]] <- beta
  }
  kmin <- which.min(BICs)
  list(beta = betas[[kmin]], lambda = lambdas[kmin],
       bic_min = BICs[kmin],
       bic_path = data.frame(lambda = lambdas, bic = BICs))
}

## ========================= Design prep (NO DROPS / NO IMPUTATION) =========================
# Pass-through: keep all columns; only ensure matrix form. No imputation, no dropping.
.prep_design <- function(X, keep_always = NULL, impute = "none",
                         drop_all_na = FALSE, drop_constant = FALSE, sd_tol = 0) {
  X <- as.matrix(X)
  # Do not touch NAs here; upstream checks already enforce finite values
  list(X = X, keep = seq_len(ncol(X)), p_full = ncol(X),
       dropped = data.frame(name=character(), reason=character()),
       kept_always = character(0))
}
nonzero_map_local <- function(Gamma, meta, cov_names, tol = 1e-8) {
  R <- nrow(Gamma); P <- ncol(Gamma)
  stopifnot(nrow(meta) == R, length(cov_names) == P)
  nz <- which(abs(Gamma) > tol, arr.ind = TRUE)
  if (!nrow(nz)) {
    return(list(
      selected_covariates = character(0),
      selection_by_resolution = meta[0, c("atom_index","j","type","kx","ky",
                                          "x_lo01","x_hi01","y_lo01","y_hi01",
                                          "x_lo","x_hi","y_lo","y_hi","desc")]
    ))
  }
  df <- cbind(
    covariate = cov_names[nz[,2]],
    meta[nz[,1], c("atom_index","j","type","kx","ky",
                   "x_lo01","x_hi01","y_lo01","y_hi01",
                   "x_lo","x_hi","y_lo","y_hi","desc")],
    stringsAsFactors = FALSE
  )
  list(
    selected_covariates = unique(df$covariate),
    selection_by_resolution = df[order(df$covariate, df$j, df$type, df$kx, df$ky), ]
  )
}


## ===================== METHODS =====================

## ---- 1) Localized LASSO (uses Haar-localized covariates) ----
local_lasso <- function(prep){
  t0  <- proc.time()[3]
  
  # Design: localized features Z = [B_U ⊙ Xq]
  Z   <- build_localized_features_safe(prep$Xq, prep$B_U)
  yi  <- prep$yi
  off <- prep$off
  
  # Optional (if you want weights in likelihood/BIC later; your select_by_bic_glmnet ignores weights)
  # wi  <- if (!is.null(prep$wi)) prep$wi else NULL
  
  # No drops/no impute
  cz <- .prep_design(Z)
  Zc <- cz$X
  
  # ----- Select lambda by BIC on glmnet path -----
  sel <- select_by_bic_glmnet(X = Zc, y = yi, offset = off, penalty.factor = NULL, family = "poisson")
  beta_sel <- sel$beta          # length = ncol(Zc)
  runtime_fit <- proc.time()[3] - t0
  
  # ----- Reconstruct Gamma (R x P) from beta_sel -----
  R <- ncol(prep$B_U)           # number of Haar atoms
  P <- ncol(prep$Xq)            # number of covariates (this equals pn)
  
  # beta_sel corresponds to vec(Gamma) with your construction:
  # columns grouped by covariate: for p=1..P, block of size R
  Gamma_vec <- numeric(R * P)
  Gamma_vec[cz$keep] <- beta_sel
  Gamma <- matrix(Gamma_vec, nrow = R, ncol = P, byrow = FALSE)
  
  # ----- β̂(s) on grid -----
  beta_hat_mat <- prep$B_grid %*% Gamma   # (n_grid x R) %*% (R x P) = n_grid x P
  
  # ----- Metrics -----
  beta_rmspe <- rmse(as.vector(beta_hat_mat), as.vector(prep$Beta))
  
  loc <- compute_local_metrics(prep$B_U, prep$B_grid, prep$Beta, Gamma,
                               normalize_to_unit_square(prep$U), pn = P,
                               p0 = sum(colSums(prep$Beta != 0) > 0))
  
  tpr2 <- compute_local_TPR_sets(
    B_U   = prep$B_U,
    U01   = normalize_to_unit_square(prep$U),
    Gamma = Gamma,
    Beta  = prep$Beta,
    tol   = 1e-8
  )
  
  list(
    beta_rmspe = beta_rmspe,
    TPR_global = loc$TPR_cov,
    FPR_global = loc$FPR_cov,
    runtime_fit = runtime_fit,
    local_TPR     = tpr2$local_TPR,
    lambda_bic        = sel$lambda,
    bic_min           = sel$bic_min
  )
}

## ============================================================
## Localized SCAD (Haar-localized covariates) -- BIC/WQBIC-style
##   - replaces cv.ncvreg with your select_by_bic_ncvreg()
##   - constructs Gamma correctly as (R x P), where P = ncol(Xq)
##   - avoids any "fit not found" and avoids wrong pn dimension
## ============================================================
local_scad <- function(prep){
  t0  <- proc.time()[3]
  
  # Design: localized features Z = [B_U ⊙ Xq]
  Z   <- build_localized_features_safe(prep$Xq, prep$B_U)
  yi  <- prep$yi
  off <- prep$off
  
  # No drops/no impute (pass-through)
  cz <- .prep_design(Z)
  Zc <- cz$X
  
  # ----- Select lambda by BIC on ncvreg path -----
  sel <- select_by_bic_ncvreg(X = Zc, y = yi, offset = off)
  beta_sel <- sel$beta          # length = ncol(Zc)
  runtime_fit <- proc.time()[3] - t0
  
  # ----- Reconstruct Gamma (R x P) from beta_sel -----
  R <- ncol(prep$B_U)           # number of Haar atoms
  P <- ncol(prep$Xq)            # number of covariates (this equals pn)
  
  Gamma_vec <- numeric(R * P)
  Gamma_vec[cz$keep] <- beta_sel
  Gamma <- matrix(Gamma_vec, nrow = R, ncol = P, byrow = FALSE)
  
  # ----- β̂(s) on grid -----
  beta_hat_mat <- prep$B_grid %*% Gamma   # (n_grid x R) %*% (R x P) = n_grid x P
  
  # ----- Metrics -----
  beta_rmspe <- rmse(as.vector(beta_hat_mat), as.vector(prep$Beta))
  
  loc <- compute_local_metrics(prep$B_U, prep$B_grid, prep$Beta, Gamma,
                               normalize_to_unit_square(prep$U), pn = P,
                               p0 = sum(colSums(prep$Beta != 0) > 0))
  
  tpr2 <- compute_local_TPR_sets(
    B_U   = prep$B_U,
    U01   = normalize_to_unit_square(prep$U),
    Gamma = Gamma,
    Beta  = prep$Beta,
    tol   = 1e-8
  )
  
  list(
    beta_rmspe = beta_rmspe,
    TPR_global = loc$TPR_cov,
    FPR_global = loc$FPR_cov,
    runtime_fit = runtime_fit,
    local_TPR     = tpr2$local_TPR,
    lambda_bic        = sel$lambda,
    bic_min           = sel$bic_min
  )
}

## ============================================================
## 3) Global LASSO (uses Xq only) -- BIC/WQBIC-style (no CV)
## ============================================================
lasso <- function(prep){
  t0  <- proc.time()[3]
  Z   <- as.matrix(prep$Xq)
  yi  <- prep$yi
  off <- prep$off
  
  cz <- .prep_design(Z)
  Zc <- cz$X
  
  sel <- select_by_bic_glmnet(X = Zc, y = yi, offset = off,
                              penalty.factor = NULL, family = "poisson")
  beta_hat <- sel$beta
  runtime_fit <- proc.time()[3] - t0
  
  p0 <- sum(colSums(prep$Beta != 0) > 0)
  S_true <- 1:p0
  S_hat  <- which(abs(beta_hat) > 1e-8)
  
  TPR <- length(intersect(S_hat, S_true)) / p0
  FPR <- length(setdiff(S_hat, S_true)) / (ncol(prep$Xgrid) - p0)
  
  beta_rmspe <- rmse(as.vector(beta_hat), as.vector(colMeans(prep$Beta)))
  
  list(
    beta_rmspe = beta_rmspe,
    TPR_global = TPR,
    FPR_global = FPR,
    runtime_fit = runtime_fit,
    local_TPR     = NA_real_,
    lambda_bic        = sel$lambda,
    bic_min           = sel$bic_min
  )
}


## ============================================================
## 4) Global SCAD (uses Xq only) -- BIC/WQBIC-style (no CV)
## ============================================================
scad <- function(prep){
  t0  <- proc.time()[3]
  Z   <- as.matrix(prep$Xq)
  yi  <- prep$yi
  off <- prep$off
  
  cz <- .prep_design(Z)
  Zc <- cz$X
  
  sel <- select_by_bic_ncvreg(X = Zc, y = yi, offset = off)
  beta_hat <- sel$beta
  runtime_fit <- proc.time()[3] - t0
  
  p0 <- sum(colSums(prep$Beta != 0) > 0)
  S_true <- 1:p0
  S_hat  <- which(abs(beta_hat) > 1e-8)
  
  TPR <- length(intersect(S_hat, S_true)) / p0
  FPR <- length(setdiff(S_hat, S_true)) / (ncol(prep$Xgrid) - p0)
  
  beta_rmspe <- rmse(as.vector(beta_hat), as.vector(colMeans(prep$Beta)))
  
  list(
    beta_rmspe = beta_rmspe,
    TPR_global = TPR,
    FPR_global = FPR,
    runtime_fit = runtime_fit,
    local_TPR     = NA_real_,
    lambda_bic        = sel$lambda,
    bic_min           = sel$bic_min
  )
}


## ============================================================
## 5) Global Adaptive LASSO (weighted LASSO on Xq) -- BIC (no CV)
##   Step 0: init ridge/lasso by BIC to get b_init
##   Step 1: weighted lasso by BIC with penalty.factor
## ============================================================
adaptive_lasso <- function(prep, gamma_adapt = 1.0, eps_w = 1e-6,
                           init = c("ridge","lasso")){
  init <- match.arg(init)
  
  t0  <- proc.time()[3]
  Z   <- as.matrix(prep$Xq)
  yi  <- prep$yi
  off <- prep$off
  
  cz <- .prep_design(Z)
  Zc <- cz$X
  
  ## ---- Step 0: initial estimator (ridge or lasso) via BIC ----
  alpha0 <- if (init == "ridge") 0 else 1
  
  # NOTE: select_by_bic_glmnet() currently hardcodes alpha=1 inside glmnet().
  # For ridge init (alpha0=0), we do a small local modification here.
  fit0 <- glmnet::glmnet(
    x = Zc, y = yi, family = "poisson", alpha = alpha0,
    standardize = FALSE, offset = off
  )
  lambdas0 <- fit0$lambda
  coefs0   <- stats::coef(fit0)
  
  BICs0 <- numeric(length(lambdas0))
  betas0 <- vector("list", length(lambdas0))
  for (k in seq_along(lambdas0)) {
    beta_full <- as.numeric(coefs0[, k])
    beta_k <- beta_full[-1]
    BICs0[k] <- bic_poisson_offset(yi, Zc, beta_k, off, n_eff = length(yi))
    betas0[[k]] <- beta_k
  }
  k0 <- which.min(BICs0)
  b_init <- betas0[[k0]]
  lam0   <- lambdas0[k0]
  
  ## adaptive weights
  w_pf <- 1 / (abs(b_init) + eps_w)^gamma_adapt
  
  ## ---- Step 1: weighted lasso via BIC ----
  sel <- select_by_bic_glmnet(X = Zc, y = yi, offset = off,
                              penalty.factor = w_pf, family = "poisson")
  beta_hat <- sel$beta
  runtime_fit <- proc.time()[3] - t0
  
  p0 <- sum(colSums(prep$Beta != 0) > 0)
  S_true <- 1:p0
  S_hat  <- which(abs(beta_hat) > 1e-8)
  
  TPR <- length(intersect(S_hat, S_true)) / p0
  FPR <- length(setdiff(S_hat, S_true)) / (ncol(prep$Xgrid) - p0)
  
  beta_rmspe <- rmse(as.vector(beta_hat), as.vector(colMeans(prep$Beta)))
  
  list(
    beta_rmspe = beta_rmspe,
    TPR_global = TPR,
    FPR_global = FPR,
    runtime_fit = runtime_fit,
    local_TPR     = NA_real_,
    lambda0_bic       = lam0,
    lambda_bic        = sel$lambda
  )
}


## ===================== BIG SIM WRAPPER =====================
bind_rows_fill <- function(...) {
  rows <- list(...)
  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    miss <- setdiff(all_names, names(x))
    if (length(miss)) x[miss] <- NA
    x[all_names]
  })
  out <- do.call(rbind, lapply(rows, as.data.frame))
  rownames(out) <- NULL
  out
}
run_all_methods <- function(
    seed        = 1,
    pn          = 1000,
    p0          = 10,
    qtiles      = c(16,16),
    J           = 2, j0 = 0,
    mu_target   = 900,
    gp_sigma2   = 1, gp_rho = 0.25,
    contrast_s  = 0.7
){
  prep <- prep_ippp(seed, pn, p0, qtiles, J, j0, mu_target, gp_sigma2, gp_rho, contrast_s)
  
  res_local_lasso <- c(list(method="local_lasso"),    local_lasso(prep))
  res_local_scad  <- c(list(method="local_scad"),     local_scad(prep))
  res_lasso       <- c(list(method="lasso"),          lasso(prep))
  res_scad        <- c(list(method="scad"),           scad(prep))
  res_adalasso    <- c(list(method="adaptive_lasso"), adaptive_lasso(prep, gamma_adapt=1.0, eps_w=1e-6, init="ridge"))
  
  out <- bind_rows_fill(res_local_lasso, res_local_scad, res_lasso, res_scad, res_adalasso)
  
  for (cn in setdiff(names(out), "method")) out[[cn]] <- as.numeric(out[[cn]])
  out
}

## ---------- Heuristic for quadrature resolution ----------
choose_qtiles <- function(n_expected, k_per_point = 5, min_dummy = 1000, max_side = 128) {
  D <- max(k_per_point * n_expected, min_dummy)
  side <- ceiling(sqrt(D))
  side <- min(side, max_side)
  c(side, side)
}
## ===================== Example run =====================
res_table <- run_all_methods(
  seed        = 1,
  pn          = 500,
  p0          = 10,
  qtiles      = c(16,16),
  J           = 2, j0 = 0,
  mu_target   = 500,
  gp_sigma2   = 1, gp_rho = 0.25,
  contrast_s  = 0.7
)
print(res_table)
```


## Large Sample Behavior
```{r}
set.seed(1)
pn <- 1000
p0 <- 10
mu_grid <- c(100, 200, 500, 900)

results <- list()
row_id <- 1L
for (mu in mu_grid) {
  qt <- c(16,16)
  cat(sprintf("\n=== mu_target = %d | qtiles = c(%d,%d) ===\n", mu, qt[1], qt[2]))
  
  tab <- run_all_methods(
    seed        = 1,
    pn          = pn,
    p0          = p0,
    qtiles      = c(16,16),
    J           = 2, j0 = 0,
    mu_target   = mu,
    gp_sigma2   = 1, gp_rho = 0.25,
    contrast_s  = 0.7
  )
  tab$mu_target   <- mu
  tab$qtiles_side <- qt[1]
  
  results[[row_id]] <- tab
  row_id <- row_id + 1L
}

results_df <- do.call(rbind, results)
# Keep key columns; include new local_TPR_new fields
keep_cols <- c("mu_target","qtiles_side","method",
               "beta_rmspe","TPR_global","FPR_global","runtime_fit",
               "local_TPR")
keep_cols <- intersect(keep_cols, names(results_df))
results_df <- results_df[, keep_cols]
row.names(results_df) <- NULL

print(results_df)


## ============================================================
## Inhomogeneous Thomas process simulator used by all methods
## Parents: Poisson(kappa_parent) on W; cluster field S(u) from Gaussian kernels.
## Intensity: lambda(u) = exp(b0 + contrast_s*z(u)) * S(u)
## ============================================================
prep_thomas <- function(seed, pn, p0, qtiles, J, j0, mu_target,
                        gp_sigma2, gp_rho, contrast_s,
                        kappa_parent = 60,
                        sigma_kernel = 0.10){
  set.seed(seed)
  W  <- owin(c(0,1), c(0,1))
  gx <- qtiles[1]; gy <- qtiles[2]
  xs <- seq(0.5/gx, 1-0.5/gx, length.out = gx)
  ys <- seq(0.5/gy, 1-0.5/gy, length.out = gy)
  grid  <- as.matrix(expand.grid(x=xs, y=ys))
  grid01 <- normalize_to_unit_square(grid)
  
  ## Truth + covariates on grid
  Xgrid1 <- simulate_gp_fields(grid, P = p0, sigma2 = gp_sigma2, rho = gp_rho)
  Xgrid2 <- matrix(0, nrow = nrow(Xgrid1), ncol = pn - p0)
  Xgrid  <- cbind(Xgrid1, Xgrid2); colnames(Xgrid) <- paste0("x", seq_len(ncol(Xgrid)))
  Beta   <- beta_truth_mat(grid01, pn = pn, p0 = p0)
  
  base_linear <- rowSums(Beta[,1:p0, drop=FALSE] * Xgrid[,1:p0, drop=FALSE])
  z   <- scale(base_linear)[,1]
  b0  <- 0
  tmp <- exp(b0 + contrast_s * z)
  b0  <- b0 + log(mu_target / mean(tmp))
  lambda_baseline_grid <- exp(b0 + contrast_s * z)
  
  ## Parents
  parents <- rpoispp(lambda = kappa_parent, win = W)
  
  ## Cluster field S(u) ≈ kernel density / kappa_parent (normalize)
  dens_im <- density.ppp(parents, sigma = sigma_kernel, at = "pixels",
                         edge = TRUE, dimyx = c(gx, gy))
  S_grid <- as.vector(dens_im$v) / kappa_parent
  
  ## Final intensity on grid
  lambda_true_grid <- pmax(lambda_baseline_grid * S_grid, 0)
  lambda_im <- im(xcol = xs, yrow = ys,
                  mat = matrix(lambda_true_grid, nrow = gx, ncol = gy, byrow = FALSE))
  
  ## Simulate realizations
  Xpp <- rpoispp(lambda = lambda_im, win = W)
  
  ## Quadrature & offsets
  Q  <- quadscheme(Xpp, method="grid", ntile=qtiles)
  isD_int <- as.logical(is.data(Q))
  wi <- try(as.numeric(weights(Q)), silent=TRUE)
  if (inherits(wi,"try-error") || length(wi)==0L) wi <- if (!is.null(Q$w)) as.numeric(Q$w) else numeric(0)
  if (length(wi)==0L && exists("w.quad", where=asNamespace("spatstat.model"), inherits=FALSE))
    wi <- as.numeric(getFromNamespace("w.quad","spatstat.model")(Q))
  Upp <- try(spatstat.geom::as.ppp(Q), silent=TRUE)
  if (!inherits(Upp,"try-error")) {
    U  <- data.frame(x=Upp$x, y=Upp$y)
    yi <- as.integer(isD_int)
  } else {
    Xd  <- if (!is.null(Q$data)) Q$data else spatstat.model::data.quad(Q)
    Xdm <- if (!is.null(Q$dummy)) Q$dummy else spatstat.model::dummy.quad(Q)
    U   <- data.frame(x=c(Xd$x, Xdm$x), y=c(Xd$y, Xdm$y))
    yi  <- c(rep(1L, npoints(Xd)), rep(0L, npoints(Xdm)))
  }
  M <- nrow(U)
  Wow <- as.owin(W)
  cell_area <- (diff(Wow$xrange)/gx) * (diff(Wow$yrange)/gy)
  if (length(wi)!=M) wi <- rep(NA_real_, M)
  bad <- !is.finite(wi) | wi<=0
  if (all(bad)) wi <- ifelse(yi==1L, 0, cell_area) else wi[bad] <- cell_area
  off <- log(pmax(wi, 1e-12))
  
  ## Map quadrature to nearest grid cell
  ix <- pmin(pmax(1, round(U$x * gx + 0.5)), gx)
  iy <- pmin(pmax(1, round(U$y * gy + 0.5)), gy)
  Xq <- Xgrid[(iy - 1) * gx + ix, , drop=FALSE]
  
  ## Haar bases
  U01 <- normalize_to_unit_square(U)
  B_qu_U   <- build_haar2d_meta(U01, j0=j0, J=J)
  B_U      <- B_qu_U$B
  B_qu_G   <- build_haar2d_meta(grid01, j0=j0, J=J)
  B_grid   <- B_qu_G$B
  R        <- ncol(B_U)
  
  list(
    W=W, gx=gx, gy=gy,
    grid=grid, grid01=grid01,
    Xgrid=Xgrid, Beta=Beta,
    lambda_true_grid=lambda_true_grid,
    U=U, yi=yi, off=off,
    ix=ix, iy=iy, Xq=Xq,
    B_U=B_U, B_grid=B_grid, R=R,
    parents=parents, sigma_kernel=sigma_kernel, kappa_parent=kappa_parent
  )
}

## ===================== BIG SIM WRAPPER =====================
bind_rows_fill <- function(...) {
  rows <- list(...)
  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    miss <- setdiff(all_names, names(x))
    if (length(miss)) x[miss] <- NA
    x[all_names]
  })
  out <- do.call(rbind, lapply(rows, as.data.frame))
  rownames(out) <- NULL
  out
}
run_all_methods <- function(
    seed        = 1,
    pn          = 1000,
    p0          = 10,
    qtiles      = c(16,16),
    J           = 2, j0 = 0,
    mu_target   = 900,
    gp_sigma2   = 1, gp_rho = 0.25,
    contrast_s  = 0.7, kappa_parent = 60,
    sigma_kernel = 0.10
){
  prep <- prep_thomas(seed, pn, p0, qtiles, J, j0, mu_target,
                      gp_sigma2, gp_rho, contrast_s,
                      kappa_parent = kappa_parent,
                      sigma_kernel = sigma_kernel)
  
  res_local_lasso <- c(list(method="local_lasso"),    local_lasso(prep))
  res_local_scad  <- c(list(method="local_scad"),     local_scad(prep))
  res_lasso       <- c(list(method="lasso"),          lasso(prep))
  res_scad        <- c(list(method="scad"),           scad(prep))
  res_adalasso    <- c(list(method="adaptive_lasso"), adaptive_lasso(prep, gamma_adapt=1.0, eps_w=1e-6, init="ridge"))
  
  out <- bind_rows_fill(res_local_lasso, res_local_scad, res_lasso, res_scad, res_adalasso)
  
  for (cn in setdiff(names(out), "method")) out[[cn]] <- as.numeric(out[[cn]])
  out
}

## ---------- Heuristic for quadrature resolution ----------
choose_qtiles <- function(n_expected, k_per_point = 5, min_dummy = 1000, max_side = 128) {
  D <- max(k_per_point * n_expected, min_dummy)
  side <- ceiling(sqrt(D))
  side <- min(side, max_side)
  c(side, side)
}
## ===================== Example run =====================
res_table_thomas <- run_all_methods(
  seed        = 1,
  pn          = 500,
  p0          = 10,
  qtiles      = c(16,16),
  J           = 2, j0 = 0,
  mu_target   = 500,
  gp_sigma2   = 1, gp_rho = 0.25,
  contrast_s  = 0.7,  kappa_parent = 60,
  sigma_kernel = 0.10
)
print(res_table_thomas)

set.seed(1)

pn <- 1000
p0 <- 10
mu_grid <- c(100, 200, 500, 900)

## ---- Regime A: Moderate clustering (more parents, wider kernel) ----
results_mod <- list(); rid <- 1L
for (mu in mu_grid) {
  qt <- c(16,16)
  cat(sprintf("\n=== MODERATE: mu_target = %d | qtiles = c(%d,%d) ===\n", mu, qt[1], qt[2]))
  tab <- run_all_methods(
    seed        = 123,
    pn          = pn,
    p0          = p0,
    qtiles      = qt,
    J           = 2, j0 = 0,
    mu_target   = mu,
    gp_sigma2   = 1, gp_rho = 0.25,
    contrast_s  = 0.7,
    kappa_parent = 80,     # more parents
    sigma_kernel = 0.12    # wider clusters
  )
  tab$mu_target   <- mu
  tab$qtiles_side <- qt[1]
  tab$clustering  <- "moderate"
  results_mod[[rid]] <- tab; rid <- rid + 1L
}
results_thomas_moderate <- do.call(rbind, results_mod)
results_thomas_moderate <- results_thomas_moderate[, c("clustering","mu_target","qtiles_side","method",
                                                       "beta_rmspe","TPR_global","FPR_global","runtime_fit","local_TPR")]

row.names(results_thomas_moderate) <- NULL
print(results_thomas_moderate)

## ---- Regime B: High clustering (fewer parents, tighter kernel) ----
results_high <- list(); rid <- 1L
for (mu in mu_grid) {
  qt <- c(16,16)
  cat(sprintf("\n=== HIGH: mu_target = %d | qtiles = c(%d,%d) ===\n", mu, qt[1], qt[2]))
  tab <- run_all_methods(
    seed        = 456,
    pn          = pn,
    p0          = p0,
    qtiles      = qt,
    J           = 2, j0 = 0,
    mu_target   = mu,
    gp_sigma2   = 1, gp_rho = 0.25,
    contrast_s  = 0.7,
    kappa_parent = 30,     # fewer parents
    sigma_kernel = 0.06    # tighter clusters
  )
  tab$mu_target   <- mu
  tab$qtiles_side <- qt[1]
  tab$clustering  <- "high"
  results_high[[rid]] <- tab; rid <- rid + 1L
}
results_thomas_high <- do.call(rbind, results_high)
results_thomas_high <- results_thomas_high[, c("clustering","mu_target","qtiles_side","method",
                                               "beta_rmspe","TPR_global","FPR_global","runtime_fit","local_TPR")]

row.names(results_thomas_high) <- NULL
print(results_thomas_high)

## Save
save(results_thomas_moderate, file = "results_thomas_moderate_simulation.RData")
save(results_thomas_high, file = "results_thomas_high_simulation.RData")

print(results_thomas_moderate)
print(results_thomas_high)

# ====================== Visualization helpers =========================
# Needs: spatstat.geom, spatstat.random, grDevices
viz_thomas <- function(prep,
                       main_prefix = "Inhomogeneous Thomas",
                       add_contours = TRUE,
                       point_cex = 0.35) {
  stopifnot(all(c("grid","gx","gy","lambda_true_grid",
                  "parents","sigma_kernel","kappa_parent",
                  "U","yi") %in% names(prep)))
  
  library(spatstat.geom)
  
  xs <- sort(unique(prep$grid[,1]))
  ys <- sort(unique(prep$grid[,2]))
  gx <- prep$gx; gy <- prep$gy
  
  # --- 1) Final intensity image lambda(u) ---
  Lambda_mat <- matrix(prep$lambda_true_grid, nrow = gx, ncol = gy, byrow = FALSE)
  Lambda_im  <- im(xcol = xs, yrow = ys, mat = Lambda_mat)
  
  # --- 2) Cluster field S(u) ---
  Sdens_im <- density.ppp(prep$parents, sigma = prep$sigma_kernel,
                          at = "pixels", edge = TRUE, dimyx = c(gx, gy))
  S_im <- eval.im(Sdens_im / prep$kappa_parent)
  
  # --- 3) Baseline lambda_base(u) ≈ lambda(u)/S(u) ---
  eps <- 1e-12
  Base_im <- eval.im(Lambda_im / pmax(S_im, eps))
  
  events_df  <- prep$U[prep$yi == 1L, , drop = FALSE]
  parents_pp <- prep$parents
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1,3), mar = c(3.5,3.5,2.5,5), oma = c(0,0,0,0))
  
  # Helper: now contour knows what to contour (the im object)
  add_layers <- function(imobj) {
    if (add_contours) {
      contour(imobj, add = TRUE, drawlabels = FALSE, nlevels = 7, lwd = 0.6)
    }
    points(parents_pp, pch = 3, cex = 0.6, col = rgb(0,0,0,0.7))
    if (nrow(events_df)) points(events_df$x, events_df$y, pch = 16, cex = point_cex, col = rgb(0,0,0,0.5))
    legend("topright", inset = 0.02, cex = 0.8, bty = "n",
           legend = c("parents", "events"),
           pch = c(3,16), col = c("black","black"))
  }
  
  plot(Base_im,
       main = sprintf("%s — Baseline λ_base(u)", main_prefix),
       ribargs = list(las = 1))
  add_layers(Base_im)
  
  plot(S_im,
       main = sprintf("%s — Cluster field S(u; C)", main_prefix),
       ribargs = list(las = 1))
  add_layers(S_im)
  
  plot(Lambda_im,
       main = sprintf("%s — Final intensity λ(u)", main_prefix),
       ribargs = list(las = 1))
  add_layers(Lambda_im)
  
  invisible(list(lambda_im = Lambda_im, S_im = S_im, base_im = Base_im))
}

# ====================== Example usage =========================
# 1) Build a scenario (your existing function)
#    e.g., "moderate" clustering:
prep <- prep_thomas(seed = 123, pn = 1000, p0 = 10, qtiles = c(16,16),
                    J = 2, j0 = 0, mu_target = 500,
                    gp_sigma2 = 1, gp_rho = 0.25, contrast_s = 0.7,
                    kappa_parent = 80, sigma_kernel = 0.12)

# 2) Visualize
viz <- viz_thomas(prep, main_prefix = "Thomas (moderate)")

# --------------------- Two-regime panel (optional) ---------------------
viz_compare_two <- function(prepA, prepB,
                            titles = c("Regime A", "Regime B"),
                            which_field = c("lambda","base","cluster")) {
  library(spatstat.geom)
  which_field <- match.arg(which_field)
  to_im <- function(prep) {
    xs <- sort(unique(prep$grid[,1])); ys <- sort(unique(prep$grid[,2]))
    gx <- prep$gx; gy <- prep$gy
    Lambda_im <- im(xcol = xs, yrow = ys,
                    mat = matrix(prep$lambda_true_grid, nrow = gx, ncol = gy, byrow = FALSE))
    Sdens_im  <- density.ppp(prep$parents, sigma = prep$sigma_kernel,
                             at = "pixels", edge = TRUE, dimyx = c(gx, gy))
    S_im      <- eval.im(Sdens_im / prep$kappa_parent)
    Base_im   <- eval.im(Lambda_im / pmax(S_im, 1e-12))
    list(lambda = Lambda_im, base = Base_im, cluster = S_im)
  }
  A <- to_im(prepA); B <- to_im(prepB)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1,2))
  plot(A[[which_field]], main = paste0(titles[1], " — ", which_field))
  points(prepA$parents, pch = 3, cex = 0.6)
  plot(B[[which_field]], main = paste0(titles[2], " — ", which_field))
  points(prepB$parents, pch = 3, cex = 0.6)
  invisible(NULL)
}

# Example:
prep_mod  <- prep_thomas(123, 1000, 10, c(16,16), 2, 0, 500, 1, 0.25, 0.7, kappa_parent=80, sigma_kernel=0.12)
prep_high <- prep_thomas(123, 1000, 10, c(16,16), 2, 0, 500, 1, 0.25, 0.7, kappa_parent=60, sigma_kernel=0.10)
viz_compare_two(prep_mod, prep_high, titles = c("Moderate", "High"), which_field = "cluster")
