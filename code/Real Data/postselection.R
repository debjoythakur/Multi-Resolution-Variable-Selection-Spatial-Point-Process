# ============================================================
# COMMENTED SUMMARY: Post-selection unpenalized BT refit +
#                    intensity comparison + block-wise maps
# ============================================================

# ------------------------------------------------------------
# GOAL (big picture)
# ------------------------------------------------------------
# 1) You already did penalized selection (Local LASSO/SCAD with BT–Haar,
#    and Global LASSO/SCAD/Adaptive LASSO) using BIC to pick lambda.
# 2) Now you do a *post-selection refit*:
#    - Take the selected variables (and selected Haar atoms for local methods),
#    - Refit an *unpenalized* Poisson BT likelihood model (MLE) using BFGS,
#    - Produce:
#       (i) MLE coefficients,
#      (ii) variable importance summaries,
#     (iii) fitted intensity lambda_hat(u) on quadrature grid,
#      (iv) plots comparing empirical intensity vs model intensities,
#       (v) block-wise (j=1 and j=2) maps showing which covariates drive which tiles.

# ============================================================
# 0) Helper functions (standardize, BT objective/gradient, refit, plotting)
# ============================================================

# standardize_cols(X):
# - Column-wise z-score standardization: (X - mean)/sd
# - Stores mu and sd as attributes for possible inverse transforms later.

# bt_loglik_and_grad(Z_design, yi, wi):
# - Defines negative Poisson log-likelihood for BT quadrature:
#     ell(w) = sum_i [ yi_i * eta_i - wi_i * exp(eta_i) ], eta = Z_design %*% w
# - Returns:
#     fn(w) = -ell(w)     (for minimization),
#     gr(w) = -∂ell/∂w    (analytic gradient for faster/stable optim)

# bt_refit(Z_design, yi, wi):
# - Runs optim(method="BFGS") starting from zeros to get unpenalized MLE.
# - Warns if convergence code != 0.
# - Returns coefficient vector w_hat (named by columns of Z_design).

# compute_intensity_and_plot(Z_design, w_hat, U, title_tag):
# - Computes eta_hat = Z_design %*% w_hat and lambda_hat = exp(eta_hat).
# - Builds a point plot over quadrature locations U with color = lambda_hat.
# - Returns lambda_hat vector (aligned with U).

# ============================================================
# 1) Build design matrices for post-selection refit
# ============================================================

# build_Z_design_local(sel_df, Z, B, Xq_reduced, method_label):
# - For LOCAL methods, selection happens at (covariate, Haar atom) level.
# - From sel_df (selection_by_resolution) it extracts unique pairs
#   (covariate, atom_index) [later version keeps full meta too].
# - Converts each pair to a column in Z via:
#       col_idx_Z = (covariate_index-1)*R + atom_index
#   where R = #Haar atoms columns in B.
# - Standardizes selected columns, then adds an intercept.
# - Returns:
#     Z_design: [Intercept + standardized selected (covariate×atom) columns]
#     parent_covariate or meta_cols: info to aggregate coefficients back to covariates.

# build_Z_design_global(sel_names, Xq_reduced, method_label):
# - For GLOBAL methods, selection is only at covariate level.
# - Keeps selected covariate columns from Xq_reduced.
# - Standardizes columns, adds intercept.
# - Returns Z_design + parent_covariate vector.

# ============================================================
# 2) Variable-importance summaries after refit
# ============================================================

# variable_importance(w_hat, parent_covariate, method_label):
# - Removes intercept, takes |beta|.
# - Aggregates by covariate name (summing abs coefficients).
# - Prints sorted importance (larger = more influence).

# variable_importance_local(w_hat, meta_cols, method_label) [later extension]:
# - Builds a full table with:
#     covariate, j, type (H/V/D), kx, ky, tile bounds + beta + |beta|
# - Prints:
#   (a) atom-level ranked table (fine-grain local importance),
#   (b) aggregated summary by (covariate, j) to show which scales matter.

# ============================================================
# 3) Run post-selection refit for EACH method
# ============================================================

# A) Localized LASSO (BT–Haar):
# - Build local Z_design using selected (covariate, atom_index) from lli_sel_cov.
# - Refit unpenalized BT-MLE with bt_refit().
# - Compute variable importance (global by covariate, and optionally by atom/level).
# - Compute lambda_hat_lli and plot the fitted intensity.

# B) Localized SCAD (BT–Haar):
# - Same pipeline as Localized LASSO, using lls_sel_cov -> w_hat_lls -> lambda_hat_lls.

# C) Global LASSO:
# - Build Z_design from selected covariates lasso_sel_cov (from global model).
# - Refit unpenalized BT-MLE -> w_hat_gl_lasso -> lambda_hat_gl_lasso -> plot.

# D) Global SCAD:
# - Same pipeline: scad_sel_cov -> w_hat_gl_scad -> lambda_hat_gl_scad -> plot.

# E) Global Adaptive LASSO:
# - Same pipeline: adlasso_sel_cov -> w_hat_gl_adalasso -> lambda_hat_gl_adalasso -> plot.

# ============================================================
# 4) Empirical intensity baseline (kernel density of event points)
# ============================================================

# sigma_emp <- bw.ppl(Xpp):
# - Automatic bandwidth for point process kernel smoothing.

# emp_field <- density(Xpp, sigma = sigma_emp, at="pixels", dimyx=c(300,300)):
# - Computes a kernel-smoothed empirical intensity field on a fine grid.

# lambda_emp_U <- emp_field[U]:
# - Extract empirical intensity at your BT quadrature locations U
#   so it aligns with lambda_hat_<method> vectors.

# compare_lambda_df:
# - A single dataframe with:
#     (x,y) + lambda_emp + lambda_hat_lli + lambda_hat_lls + lambda_hat_global_* ...

# ============================================================
# 5) Side-by-side intensity comparison panels (2×3)
# ============================================================

# make_intensity_point_plot(df, lambda_col, title_tag):
# - Generic plot: points at (x,y), color by chosen lambda column.

# Build panels:
# - p_emp, p_lli, p_lls, p_gl_lasso, p_gl_scad, p_gl_adalasso

# Arrange with patchwork:
# - Top row: empirical vs local LASSO vs local SCAD
# - Bottom row: global LASSO vs global SCAD vs global adaptive LASSO

# ============================================================
# 6) Make “block-wise” local interpretation maps on the roads
# ============================================================

# Step 6.1: Attach refit coefficients to local meta
# - meta_beta_lli = meta_lli + beta_lli and absbeta (same for SCAD as meta_beta_lls)
# - This gives covariate + (j,type,kx,ky,tile) + estimated beta.

# Step 6.2: Normalize road geometry to [0,1]^2 for easy overlay of Haar blocks
# - st_geometry(roads_events_01) = (geom - min) / (range)  (componentwise)
# - So Haar tiles (kx,ky) match simple unit-square partitions.

# Step 6.3: Build label text for each Haar tile
# - make_block_labels(meta_beta, j_level, kx_vals, ky_vals, top_per_block):
#   For each (kx,ky) tile at level j:
#     - subset coefficients in that tile
#     - split by type: H(psi_x), V(psi_y), D(psi_xy)
#     - within each type, aggregate beta by covariate and take top few (by |beta|)
#     - create a multi-line label like:
#         H(kx,ky): cov1 (0.123) / cov2 (-0.045)
#         V(kx,ky): ...
#         D(kx,ky): ...

# Step 6.4: Plot j=1 (2×2) and j=2 (4×4) labeled block maps
# - blocks_j1 / blocks_j2 from make_block_labels()
# - ggplot(roads_events_01_crisp or roads_events_01) +
#     geom_sf(roads) +
#     geom_rect(block boundaries) +
#     geom_label(block labels)
# - Repeat for LLI (local lasso) and LLS (local scad) separately.

# Step 6.5 (visual clarity): “crisp” road plotting fix
# - st_make_valid() for geometry
# - st_simplify(dTolerance=5) to reduce vertex clutter
# - slice_sample() to cap number of road segments (speed + readability)
# - normalize only after simplification/sampling

# ============================================================
# 7) Optional: Add simple numbered grids (2×2 or 4×4) on original roads plot
# ============================================================

# add_grid_and_labels(p, roads_sf, n=2 or 4):
# - Computes bounding box of roads_sf.
# - Draws vertical/horizontal lines splitting bbox into n×n grid.
# - Places a numbered bubble in the top-left of each cell for reference.

# p_base:
# - Your base road plot: grey vs black by has_incident, thin vs thick linewidth.

# p_4, p_16:
# - Adds 2×2 or 4×4 grids with bright label bubbles for easy “cell indexing”.

# ============================================================
# FINAL PRODUCTS from this whole script
# ============================================================
# - Unpenalized BT refit coefficients for each method: w_hat_*
# - Variable importance summaries: vi_* (and local atom/level tables if used)
# - Fitted intensities on quadrature points: lambda_hat_*
# - Empirical baseline intensity: lambda_emp_U
# - 6-panel intensity comparison figure (empirical + 5 methods)
# - j=1 and j=2 local interpretation maps (tiles labeled with top covariates + betas)
# - Optional numbered grid overlays (2×2 and 4×4) for quick referencing.
# ============================================================

## ==========================================================
## Post-selection *unpenalized* BT refit for all methods
##   Assumes objects are already in memory:
##   - Z, B, Xq_reduced, yi, wi, U
##   - lli_sel_cov, lls_sel_cov, lasso_sel_cov, scad_sel_cov, adlasso_sel_cov
## ==========================================================

library(dplyr)
library(ggplot2)

## ---------- 0) Helpers: standardization, BT loglik, plotting ----------

standardize_cols <- function(X) {
  X <- as.matrix(X)
  mu <- colMeans(X)
  sd <- apply(X, 2, sd)
  sd[sd == 0] <- 1
  Xs <- sweep(sweep(X, 2, mu, "-"), 2, sd, "/")
  attr(Xs, "mu") <- mu
  attr(Xs, "sd") <- sd
  Xs
}

bt_loglik_and_grad <- function(Z_design, yi, wi) {
  # returns functions (fn, gr) to be used by optim
  fn <- function(w) {
    eta <- as.numeric(Z_design %*% w)
    # negative log-likelihood for minimization
    -sum(yi * eta - wi * exp(eta))
  }
  gr <- function(w) {
    eta <- as.numeric(Z_design %*% w)
    resid <- yi - wi * exp(eta)
    # gradient of *negative* loglik
    -drop(t(Z_design) %*% resid)
  }
  list(fn = fn, gr = gr)
}

bt_refit <- function(Z_design, yi, wi, method_label = "") {
  M <- nrow(Z_design)
  stopifnot(length(yi) == M, length(wi) == M)
  
  w0 <- rep(0, ncol(Z_design))
  lg <- bt_loglik_and_grad(Z_design, yi, wi)
  
  opt <- optim(
    par     = w0,
    fn      = lg$fn,
    gr      = lg$gr,
    method  = "BFGS",
    control = list(maxit = 200, reltol = 1e-8)
  )
  
  if (opt$convergence != 0) {
    warning("BT optimization for ", method_label,
            " did not fully converge (code = ", opt$convergence, ")")
  }
  
  w_hat <- opt$par
  names(w_hat) <- colnames(Z_design)
  w_hat
}

compute_intensity_and_plot <- function(Z_design, w_hat, U, title_tag) {
  eta_hat    <- as.numeric(Z_design %*% w_hat)
  lambda_hat <- exp(eta_hat)
  
  df_int <- data.frame(
    x      = U$x,
    y      = U$y,
    lambda = lambda_hat
  )
  
  p <- ggplot(df_int, aes(x = x, y = y, color = lambda)) +
    geom_point(size = 0.5) +
    coord_equal() +
    scale_color_viridis_c() +
    theme_minimal() +
    labs(
      title = paste("BT intensity –", title_tag, "(unpenalized refit)"),
      color = expression(hat(lambda)(u))
    )
  
  print(p)
  
  lambda_hat
}
## ---------- 1) Local methods: shared builder for Z_design ----------

build_Z_design_local <- function(sel_df, Z, B, Xq_reduced, method_label) {
  if (is.null(sel_df) || nrow(sel_df) == 0) {
    message("No localized selections for ", method_label)
    return(NULL)
  }
  
  R         <- ncol(B)
  cov_names <- colnames(Xq_reduced)
  
  # unique (covariate, atom_index) pairs
  pairs <- sel_df %>%
    distinct(covariate, atom_index)
  
  p_idx <- match(pairs$covariate, cov_names)
  if (any(is.na(p_idx))) {
    stop("Some covariate names in ", method_label,
         " not found in Xq_reduced: ",
         paste(unique(pairs$covariate[is.na(p_idx)]), collapse = ", "))
  }
  
  col_idx_Z <- (p_idx - 1L) * R + pairs$atom_index
  
  Z_sel <- Z[, col_idx_Z, drop = FALSE]
  colnames(Z_sel) <- paste0(pairs$covariate, "_a", pairs$atom_index)
  
  # standardize and add intercept
  Z_sel_std <- standardize_cols(Z_sel)
  Z_design  <- cbind("(Intercept)" = 1, Z_sel_std)
  colnames(Z_design)[1] <- "(Intercept)"
  
  list(
    Z_design        = Z_design,
    parent_covariate = pairs$covariate  # one per *non-intercept* column
  )
}

## ---------- 2) Global methods: shared builder for X-design ----------

build_Z_design_global <- function(sel_names, Xq_reduced, method_label) {
  if (is.null(sel_names) || length(sel_names) == 0) {
    message("No global selections for ", method_label)
    return(NULL)
  }
  
  keep <- intersect(sel_names, colnames(Xq_reduced))
  if (length(keep) == 0) {
    message("Selected names for ", method_label,
            " not found in Xq_reduced.")
    return(NULL)
  }
  
  X_sel     <- Xq_reduced[, keep, drop = FALSE]
  X_sel_std <- standardize_cols(X_sel)
  Z_design  <- cbind("(Intercept)" = 1, X_sel_std)
  colnames(Z_design)[1] <- "(Intercept)"
  
  list(
    Z_design        = Z_design,
    parent_covariate = keep  # one per *non-intercept* column
  )
}

## ---------- 3) Variable importance (aggregate over atoms if needed) ----------

variable_importance <- function(w_hat, parent_covariate, method_label) {
  beta_hat <- w_hat[names(w_hat) != "(Intercept)"]
  
  # make sure lengths match
  stopifnot(length(beta_hat) == length(parent_covariate))
  
  vi <- tapply(abs(beta_hat), parent_covariate, sum)
  vi <- sort(vi, decreasing = TRUE)
  
  cat("\n===== ", method_label, " – variable importance (BT refit) =====\n", sep = "")
  print(vi)
  
  vi
}

## ==========================================================
## A. Localized LASSO (BT–Haar)
## ==========================================================

loc_lasso_obj <- build_Z_design_local(
  sel_df       = lli_sel_cov,
  Z            = Z,
  B            = B,
  Xq_reduced   = Xq_reduced,
  method_label = "Localized LASSO (BT–Haar)"
)

if (!is.null(loc_lasso_obj)) {
  Z_design_lli <- loc_lasso_obj$Z_design
  parent_cov_lli <- loc_lasso_obj$parent_covariate
  
  w_hat_lli <- bt_refit(Z_design_lli, yi, wi,
                        method_label = "Localized LASSO (BT–Haar)")
  vi_lli <- variable_importance(w_hat_lli, parent_cov_lli,
                                method_label = "Localized LASSO (BT–Haar)")
  
  lambda_hat_lli <- compute_intensity_and_plot(
    Z_design = Z_design_lli,
    w_hat    = w_hat_lli,
    U        = U,
    title_tag = "Localized LASSO (BT–Haar)"
  )
}

## ==========================================================
## B. Localized SCAD (BT–Haar)
## ==========================================================

loc_scad_obj <- build_Z_design_local(
  sel_df       = lls_sel_cov,
  Z            = Z,
  B            = B,
  Xq_reduced   = Xq_reduced,
  method_label = "Localized SCAD (BT–Haar)"
)

if (!is.null(loc_scad_obj)) {
  Z_design_lls   <- loc_scad_obj$Z_design
  parent_cov_lls <- loc_scad_obj$parent_covariate
  
  w_hat_lls <- bt_refit(Z_design_lls, yi, wi,
                        method_label = "Localized SCAD (BT–Haar)")
  vi_lls <- variable_importance(w_hat_lls, parent_cov_lls,
                                method_label = "Localized SCAD (BT–Haar)")
  
  lambda_hat_lls <- compute_intensity_and_plot(
    Z_design = Z_design_lls,
    w_hat    = w_hat_lls,
    U        = U,
    title_tag = "Localized SCAD (BT–Haar)"
  )
}

## ==========================================================
## C. Global LASSO
## ==========================================================

gl_lasso_obj <- build_Z_design_global(
  sel_names    = lasso_sel_cov,
  Xq_reduced   = Xq_reduced,
  method_label = "Global LASSO"
)

if (!is.null(gl_lasso_obj)) {
  Z_design_gl_lasso   <- gl_lasso_obj$Z_design
  parent_cov_gl_lasso <- gl_lasso_obj$parent_covariate
  
  w_hat_gl_lasso <- bt_refit(Z_design_gl_lasso, yi, wi,
                             method_label = "Global LASSO")
  vi_gl_lasso <- variable_importance(w_hat_gl_lasso, parent_cov_gl_lasso,
                                     method_label = "Global LASSO")
  
  lambda_hat_gl_lasso <- compute_intensity_and_plot(
    Z_design = Z_design_gl_lasso,
    w_hat    = w_hat_gl_lasso,
    U        = U,
    title_tag = "Global LASSO"
  )
}

## ==========================================================
## D. Global SCAD
## ==========================================================

gl_scad_obj <- build_Z_design_global(
  sel_names    = scad_sel_cov,
  Xq_reduced   = Xq_reduced,
  method_label = "Global SCAD"
)

if (!is.null(gl_scad_obj)) {
  Z_design_gl_scad   <- gl_scad_obj$Z_design
  parent_cov_gl_scad <- gl_scad_obj$parent_covariate
  
  w_hat_gl_scad <- bt_refit(Z_design_gl_scad, yi, wi,
                            method_label = "Global SCAD")
  vi_gl_scad <- variable_importance(w_hat_gl_scad, parent_cov_gl_scad,
                                    method_label = "Global SCAD")
  
  lambda_hat_gl_scad <- compute_intensity_and_plot(
    Z_design = Z_design_gl_scad,
    w_hat    = w_hat_gl_scad,
    U        = U,
    title_tag = "Global SCAD"
  )
}

## ==========================================================
## E. Global Adaptive LASSO
## ==========================================================

gl_adalasso_obj <- build_Z_design_global(
  sel_names    = adlasso_sel_cov,
  Xq_reduced   = Xq_reduced,
  method_label = "Global adaptive LASSO"
)

if (!is.null(gl_adalasso_obj)) {
  Z_design_gl_adalasso   <- gl_adalasso_obj$Z_design
  parent_cov_gl_adalasso <- gl_adalasso_obj$parent_covariate
  
  w_hat_gl_adalasso <- bt_refit(Z_design_gl_adalasso, yi, wi,
                                method_label = "Global adaptive LASSO")
  vi_gl_adalasso <- variable_importance(w_hat_gl_adalasso, parent_cov_gl_adalasso,
                                        method_label = "Global adaptive LASSO")
  
  lambda_hat_gl_adalasso <- compute_intensity_and_plot(
    Z_design = Z_design_gl_adalasso,
    w_hat    = w_hat_gl_adalasso,
    U        = U,
    title_tag = "Global adaptive LASSO"
  )
}


# Kernel bandwidth (automatic or manual)
sigma_emp <- bw.ppl(Xpp)   # recommended
# sigma_emp <- 0.005       # you can set manually too

# Create density object over a fine pixel grid
emp_field <- density(Xpp, sigma = sigma_emp, at = "pixels",
                     dimyx = c(300, 300))   # dense grid for accuracy

# Now extract λ_emp(U) for your BT quadrature grid
lambda_emp_U <- emp_field[U]

compare_lambda_df <- data.frame(
  x = U[,1],
  y = U[,2],
  lambda_emp   = lambda_emp_U,
  lambda_hat_lli = lambda_hat_lli,
  lambda_hat_lls = lambda_hat_lls,
  lambda_hat_gl_lasso = lambda_hat_gl_lasso,
  lambda_hat_gl_scad = lambda_hat_gl_scad,
  lambda_hat_gl_adalasso = lambda_hat_gl_adalasso
)


library(ggplot2)
library(patchwork)
## ----------------------------------------------------------
## 1) Helper: one BT-style scatter intensity plot
## ----------------------------------------------------------
make_intensity_point_plot <- function(df, lambda_col, title_tag) {
  ggplot(df, aes(x = x, y = y, color = .data[[lambda_col]])) +
    geom_point(size = 0.5) +
    coord_equal() +
    scale_color_gradient(low = "grey95", high = "grey10") +
    theme_bw() +
    labs(
      title = title_tag,
      color = expression(hat(lambda)(u)),
      x = "x", y = "y"
    )
}

## ----------------------------------------------------------
## 2) Build the six panels
## ----------------------------------------------------------
p_emp <- make_intensity_point_plot(
  compare_lambda_df,
  "lambda_emp",
  "Empirical kernel intensity"
)

p_lli <- make_intensity_point_plot(
  compare_lambda_df,
  "lambda_hat_lli",
  "Local Lasso Based Intensity"
)

p_lls <- make_intensity_point_plot(
  compare_lambda_df,
  "lambda_hat_lls",
  "Local SCAD Based Intensity"
)

p_gl_lasso <- make_intensity_point_plot(
  compare_lambda_df,
  "lambda_hat_gl_lasso",
  "Lasso Based Intensity"
)

p_gl_scad <- make_intensity_point_plot(
  compare_lambda_df,
  "lambda_hat_gl_scad",
  "Scad Based Based Intensity"
)

p_gl_adalasso <- make_intensity_point_plot(
  compare_lambda_df,
  "lambda_hat_gl_adalasso",
  "Adaptive Lasso Based Intensity"
)

## ----------------------------------------------------------
## 3) Arrange in 2 rows × 3 columns
## ----------------------------------------------------------
p_all <- (p_emp | p_lli | p_lls) /
  (p_gl_lasso | p_gl_scad | p_gl_adalasso)

print(p_all)


## ==========================================================
## Done: you now have, for each method:
##  - w_hat_<method>     : BT-MLE coefficients (unpenalized)
##  - vi_<method>        : variable importance by covariate
##  - lambda_hat_<method>: intensity at quadrature points
##  - scatter plot of intensity over quadrature points
## ==========================================================
build_Z_design_local <- function(sel_df, Z, B, Xq_reduced, method_label) {
  if (is.null(sel_df) || nrow(sel_df) == 0) {
    message("No localized selections for ", method_label)
    return(NULL)
  }
  
  R         <- ncol(B)
  cov_names <- colnames(Xq_reduced)
  
  # Use ALL meta columns from sel_df (covariate, atom_index, j, type, kx, ky, ...)
  pairs <- sel_df %>%
    distinct(covariate, atom_index, j, type, kx, ky,
             x_lo01, x_hi01, y_lo01, y_hi01,
             x_lo, x_hi, y_lo, y_hi)
  
  p_idx <- match(pairs$covariate, cov_names)
  if (any(is.na(p_idx))) {
    stop("Some covariate names in ", method_label,
         " not found in Xq_reduced: ",
         paste(unique(pairs$covariate[is.na(p_idx)]), collapse = ", "))
  }
  
  col_idx_Z <- (p_idx - 1L) * R + pairs$atom_index
  
  Z_sel <- Z[, col_idx_Z, drop = FALSE]
  colnames(Z_sel) <- paste0(pairs$covariate, "_a", pairs$atom_index)
  
  # standardize and add intercept
  Z_sel_std <- standardize_cols(Z_sel)
  Z_design  <- cbind("(Intercept)" = 1, Z_sel_std)
  colnames(Z_design)[1] <- "(Intercept)"
  
  list(
    Z_design  = Z_design,
    meta_cols = pairs  # full local meta for each *non-intercept* column
  )
}


variable_importance_local <- function(w_hat, meta_cols, method_label) {
  beta_hat <- w_hat[names(w_hat) != "(Intercept)"]
  
  if (length(beta_hat) != nrow(meta_cols)) {
    stop("Length of beta_hat does not match number of meta rows in ", method_label)
  }
  
  meta_beta <- meta_cols %>%
    mutate(
      beta    = as.numeric(beta_hat),
      absbeta = abs(beta_hat)
    )
  
  # (a) Full local table: each (covariate, j, type, kx, ky, tile)
  cat("\n===== ", method_label,
      " – local variable importance by atom (BT refit) =====\n", sep = "")
  print(
    meta_beta %>%
      arrange(desc(absbeta))
  )
  
  # (b) Aggregated at “local level”: covariate × Haar level j
  vi_level <- meta_beta %>%
    group_by(covariate, j) %>%
    summarise(
      sum_absbeta = sum(absbeta),
      n_atoms     = n(),
      .groups     = "drop"
    ) %>%
    arrange(covariate, j)
  
  cat("\n===== ", method_label,
      " – aggregated importance by (covariate, level j) =====\n", sep = "")
  print(vi_level)
  
  invisible(list(atom_level = meta_beta, level_summary = vi_level))
}

## ==========================================================
## A. Localized LASSO (BT–Haar)
## ==========================================================

loc_lasso_obj <- build_Z_design_local(
  sel_df       = lli_sel_cov,
  Z            = Z,
  B            = B,
  Xq_reduced   = Xq_reduced,
  method_label = "Localized LASSO (BT–Haar)"
)

if (!is.null(loc_lasso_obj)) {
  Z_design_lli <- loc_lasso_obj$Z_design
  meta_lli     <- loc_lasso_obj$meta_cols   # <--- local meta
  
  w_hat_lli <- bt_refit(
    Z_design_lli, yi, wi,
    method_label = "Localized LASSO (BT–Haar)"
  )
  
  vi_lli_local <- variable_importance_local(
    w_hat_lli, meta_lli,
    method_label = "Localized LASSO (BT–Haar)"
  )
  
  lambda_hat_lli <- compute_intensity_and_plot(
    Z_design = Z_design_lli,
    w_hat    = w_hat_lli,
    U        = U,
    title_tag = "Localized LASSO (BT–Haar)"
  )
}

## ==========================================================
## B. Localized SCAD (BT–Haar)
## ==========================================================

loc_scad_obj <- build_Z_design_local(
  sel_df       = lls_sel_cov,
  Z            = Z,
  B            = B,
  Xq_reduced   = Xq_reduced,
  method_label = "Localized SCAD (BT–Haar)"
)

if (!is.null(loc_scad_obj)) {
  Z_design_lls <- loc_scad_obj$Z_design
  meta_lls     <- loc_scad_obj$meta_cols
  
  w_hat_lls <- bt_refit(
    Z_design_lls, yi, wi,
    method_label = "Localized SCAD (BT–Haar)"
  )
  
  vi_lls_local <- variable_importance_local(
    w_hat_lls, meta_lls,
    method_label = "Localized SCAD (BT–Haar)"
  )
  
  lambda_hat_lls <- compute_intensity_and_plot(
    Z_design = Z_design_lls,
    w_hat    = w_hat_lls,
    U        = U,
    title_tag = "Localized SCAD (BT–Haar)"
  )
}
library(sf)
library(dplyr)
library(ggplot2)

## ==========================================================
## 0) Build meta_beta_lli (local LASSO meta + coeffs)
##    (assumes you already ran build_Z_design_local and bt_refit)
## ==========================================================
meta_lli <- loc_lasso_obj$meta_cols
beta_lli <- w_hat_lli[names(w_hat_lli) != "(Intercept)"]

stopifnot(nrow(meta_lli) == length(beta_lli))

meta_beta_lli <- meta_lli %>%
  mutate(
    beta    = as.numeric(beta_lli),
    absbeta = abs(beta)
  )

## ==========================================================
## 1) NORMALIZE GEOMETRY TO [0,1]^2
## ==========================================================
bb   <- st_bbox(roads_events)
xr   <- c(bb["xmin"], bb["xmax"])
yr   <- c(bb["ymin"], bb["ymax"])
dx   <- xr[2] - xr[1]
dy   <- yr[2] - yr[1]

roads_events_01 <- roads_events
st_geometry(roads_events_01) <-
  (st_geometry(roads_events_01) - c(xr[1], yr[1])) / c(dx, dy)

## Helper to build labels for a given j and block grid
make_block_labels <- function(meta_beta, j_level, kx_vals, ky_vals,
                              top_per_block = 2) {
  blocks <- expand.grid(kx = kx_vals, ky = ky_vals)
  
  # normalized block coordinates in [0,1]
  xs01 <- seq(0, 1, length.out = length(kx_vals) + 1)
  ys01 <- seq(0, 1, length.out = length(ky_vals) + 1)
  
  blocks <- blocks %>%
    mutate(
      xmin = xs01[kx + 1],
      xmax = xs01[kx + 2],
      ymin = ys01[ky + 1],
      ymax = ys01[ky + 2],
      cx   = (xmin + xmax) / 2,
      cy   = (ymin + ymax) / 2
    )
  
  labels <- character(nrow(blocks))
  
  for (i in seq_len(nrow(blocks))) {
    kx_i <- blocks$kx[i]
    ky_i <- blocks$ky[i]
    
    sel <- subset(meta_beta, j == j_level & kx == kx_i & ky == ky_i)
    
    make_line <- function(type_code, letter) {
      s <- sel[sel$type == type_code, , drop = FALSE]
      if (nrow(s) == 0) {
        # no j-level in the label
        return(sprintf("%s(%d,%d): none", letter, kx_i, ky_i))
      }
      # aggregate by covariate within this tile
      agg <- aggregate(beta ~ covariate, data = s, FUN = sum)
      agg <- agg[order(-abs(agg$beta)), , drop = FALSE]
      agg <- head(agg, top_per_block)
      
      agg$beta_str <- sprintf("%.3f", agg$beta)
      cov_lines <- paste0("  ", agg$covariate, " (", agg$beta_str, ")",
                          collapse = "\n")
      
      # again only (kx, ky)
      sprintf("%s(%d,%d):\n%s", letter, kx_i, ky_i, cov_lines)
    }
    
    
    H_line <- make_line("psi_x",  "H")
    V_line <- make_line("psi_y",  "V")
    D_line <- make_line("psi_xy", "D")
    
    labels[i] <- paste(H_line, V_line, D_line, sep = "\n")
  }
  
  blocks$label <- labels
  blocks
}

## ==========================================================
## 2) j = 1  → 2×2 blocks in [0,1]^2
## ==========================================================
blocks_j1 <- make_block_labels(
  meta_beta = meta_beta_lli,
  j_level   = 1,
  kx_vals   = 0:1,
  ky_vals   = 0:1,
  top_per_block = 2
)

ggplot(roads_events_01) +
  geom_sf(aes(color = factor(has_incident)), linewidth = 0.3, na.rm = TRUE) +
  scale_color_manual(
    values = c("0" = "grey80", "1" = "red"),
    labels = c("No", "Yes"),
    name   = "Has incident"
  ) +
  theme_minimal() +
  labs(title = "Roads with Any Incident (accident or crime)")


# ------------------------------------------------------------
# Single-snippet fix: make roads_events_01 crisp (no haze)
# by (i) keeping UTM geometry, (ii) simplifying + sampling,
# and (iii) then normalizing to [0,1] only for plotting.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
})

# 0) Start from the *clear* object with valid CRS
roads <- roads_events %>% st_make_valid()

# 1) OPTIONAL BUT HIGHLY RECOMMENDED: simplify + sample (speed + clarity)
#    dTolerance is in meters (UTM). Try 2, 5, 10 depending on desired detail.
roads_plot <- roads %>%
  st_simplify(dTolerance = 5, preserveTopology = TRUE)

set.seed(1)
roads_plot <- roads_plot %>% slice_sample(n = min(nrow(roads_plot), 50000))

# 2) Normalize to unit square [0,1] only for plotting
bb <- st_bbox(roads_plot)
xr <- as.numeric(bb["xmax"] - bb["xmin"])
yr <- as.numeric(bb["ymax"] - bb["ymin"])

roads_events_01_crisp <- roads_plot %>%
  mutate(
    geom01 = (st_geometry(.) - c(bb["xmin"], bb["ymin"])) / c(xr, yr)
  ) %>%
  st_set_geometry("geom01") %>%
  st_set_crs(NA)

# 3) Plot (ggplot) + blocks overlay (assumes blocks_j1 already in [0,1] coords)
p_blocks_j1 <- ggplot(roads_events_01_crisp) +
  geom_sf(aes(color = factor(has_incident)),
          linewidth = 0.35, alpha = 0.9, na.rm = TRUE) +
  scale_color_manual(
    values = c("0" = "black", "1" = "red"),
    labels = c("No", "Yes"),
    name   = "Has incident"
  ) +
  coord_sf(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    legend.position  = "right",
    plot.title       = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Selected variables and estimated regression coefficient at j=1 by LLI",
    x = NULL, y = NULL
  ) +
  geom_rect(
    data = blocks_j1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "black", linewidth = 0.8
  ) +
  geom_label(
    data = blocks_j1,
    aes(x = cx, y = cy, label = label),
    inherit.aes = FALSE,
    size = 2.8,
    fill = "blue", colour = "white",
    label.size = 0
  )

print(p_blocks_j1)


## ==========================================================
## 3) j = 2  → 4×4 blocks in [0,1]^2
## ==========================================================
blocks_j2 <- make_block_labels(
  meta_beta = meta_beta_lli,
  j_level   = 2,
  kx_vals   = 0:3,
  ky_vals   = 0:3,
  top_per_block = 2
)

p_blocks_j2 <- ggplot(roads_events_01) +
  geom_sf(aes(color = factor(has_incident)),
          linewidth = 0.35, alpha = 0.9, na.rm = TRUE) +
  scale_color_manual(
    values = c("0" = "grey80", "1" = "red"),
    labels = c("No", "Yes"),
    name   = "Has incident"
  ) +
  coord_sf(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    legend.position  = "right",
    plot.title       = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Selected variables and estimated regression coefficeint at j=2 by LLI",
    x = NULL, y = NULL
  ) +
  geom_rect(
    data = blocks_j2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill   = NA,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_label(
    data = blocks_j2,
    aes(x = cx, y = cy, label = label),
    inherit.aes = FALSE,
    size = 2.2,
    fill = "blue",
    colour = "white",
    label.size = 0
  )

print(p_blocks_j2)




library(sf)
library(dplyr)
library(ggplot2)

## ==========================================================
## 0) meta_beta_lls (local SCAD meta + coefficients)
##    assumes loc_scad_obj & w_hat_lls already exist
## ==========================================================
meta_lls <- loc_scad_obj$meta_cols
beta_lls <- w_hat_lls[names(w_hat_lls) != "(Intercept)"]

stopifnot(nrow(meta_lls) == length(beta_lls))

meta_beta_lls <- meta_lls %>%
  mutate(
    beta    = as.numeric(beta_lls),
    absbeta = abs(beta)
  )

## ==========================================================
## 1) Reuse same normalization roads_events_01
##    (if not in memory, repeat the normalization)
## ==========================================================
# (Assuming roads_events_01 already created in part A)

## ==========================================================
## 2) j = 1 blocks → 2×2 grid (same as above)
## ==========================================================
j_level <- 1
xs <- seq(0, 1, length.out = 3)
ys <- seq(0, 1, length.out = 3)

blocks_j1_scad <- expand.grid(kx = 0:1, ky = 0:1)
blocks_j1_scad <- blocks_j1_scad %>%
  mutate(
    xmin = xs[kx + 1],
    xmax = xs[kx + 2],
    ymin = ys[ky + 1],
    ymax = ys[ky + 2],
    cx   = (xmin + xmax) / 2,
    cy   = (ymin + ymax) / 2
  )

## ==========================================================
## 3) Labels per j = 1 tile (SCAD)
## ==========================================================
top_per_block <- 2
labels_j1_lls <- character(nrow(blocks_j1_scad))

for (i in seq_len(nrow(blocks_j1_scad))) {
  kx_i <- blocks_j1_scad$kx[i]
  ky_i <- blocks_j1_scad$ky[i]
  
  sel <- subset(meta_beta_lls, j == j_level & kx == kx_i & ky == ky_i)
  
  make_line <- function(type_code, letter) {
    s <- sel[sel$type == type_code, , drop = FALSE]
    if (nrow(s) == 0) {
      # no j-level in the label
      return(sprintf("%s(%d,%d): none", letter, kx_i, ky_i))
    }
    # aggregate by covariate within this tile
    agg <- aggregate(beta ~ covariate, data = s, FUN = sum)
    agg <- agg[order(-abs(agg$beta)), , drop = FALSE]
    agg <- head(agg, top_per_block)
    
    agg$beta_str <- sprintf("%.3f", agg$beta)
    cov_lines <- paste0("  ", agg$covariate, " (", agg$beta_str, ")",
                        collapse = "\n")
    
    # again only (kx, ky)
    sprintf("%s(%d,%d):\n%s", letter, kx_i, ky_i, cov_lines)
  }
  
  
  H_line <- make_line("psi_x",  "H")
  V_line <- make_line("psi_y",  "V")
  D_line <- make_line("psi_xy", "D")
  
  labels_j1_lls[i] <- paste(H_line, V_line, D_line, sep = "\n")
}

blocks_j1_scad$label_lls <- labels_j1_lls

## ==========================================================
## 4) Plot: j = 1 blocks for Localized SCAD
## ==========================================================
p_blocks_j1_lls <- ggplot(roads_events_01) +
  geom_sf(aes(color = factor(has_incident)),
          linewidth = 0.35,
          alpha = 0.9,
          na.rm = TRUE) +
  scale_color_manual(
    values = c("0" = "grey80", "1" = "red"),
    labels = c("No", "Yes"),
    name   = "Has incident"
  ) +
  coord_sf(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    legend.position  = "right",
    plot.title       = element_text(hjust = 0.5)
  ) +
  labs(
    title =  "Selected variables and estimated regression coefficeint at j=1 by LLS",
    x = NULL, y = NULL
  ) +
  geom_rect(
    data = blocks_j1_scad,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill   = NA,
    colour = "black",
    linewidth = 0.7
  ) +
  geom_label(
    data = blocks_j1_scad,
    aes(x = cx, y = cy, label = label_lls),
    inherit.aes = FALSE,
    size = 2.4,
    fill = "blue",
    colour = "white",
    label.size = 0
  )

print(p_blocks_j1_lls)


## ==========================================================
## 0) Build meta_beta_lls (local SCAD meta + coefficients)
##     (assumes loc_scad_obj and w_hat_lls already exist)
## ==========================================================
meta_lls <- loc_scad_obj$meta_cols                 # from build_Z_design_local()
beta_lls <- w_hat_lls[names(w_hat_lls) != "(Intercept)"]

stopifnot(nrow(meta_lls) == length(beta_lls))

meta_beta_lls <- meta_lls %>%
  mutate(
    beta    = as.numeric(beta_lls),
    absbeta = abs(beta)
  )

## ==========================================================
## 1) Normalize roads to [0,1]^2 so region fills all blocks
## ==========================================================
bb   <- st_bbox(roads_events)
xr   <- c(bb["xmin"], bb["xmax"])
yr   <- c(bb["ymin"], bb["ymax"])
dx   <- xr[2] - xr[1]
dy   <- yr[2] - yr[1]

roads_events_01 <- roads_events
st_geometry(roads_events_01) <-
  (st_geometry(roads_events_01) - c(xr[1], yr[1])) / c(dx, dy)

## ==========================================================
## 2) j = 2 blocks in [0,1]^2  → 4×4 grid
## ==========================================================
xs <- seq(0, 1, length.out = 5)
ys <- seq(0, 1, length.out = 5)

blocks <- expand.grid(kx = 0:3, ky = 0:3)
blocks <- blocks %>%
  mutate(
    xmin = xs[kx + 1],
    xmax = xs[kx + 2],
    ymin = ys[ky + 1],
    ymax = ys[ky + 2],
    cx   = (xmin + xmax) / 2,
    cy   = (ymin + ymax) / 2
  )

## ==========================================================
## 3) Build label for each j = 2 tile (SCAD):
##    H/V/D with covariates + coefficients
## ==========================================================
top_per_block <- 2   # covariates per H/V/D line
labels <- character(nrow(blocks))

for (i in seq_len(nrow(blocks))) {
  kx_i <- blocks$kx[i]
  ky_i <- blocks$ky[i]
  
  # all local atoms in this (j=2, kx, ky) tile
  sel <- subset(meta_beta_lls, j == 2 & kx == kx_i & ky == ky_i)
  
  make_line <- function(type_code, letter) {
    s <- sel[sel$type == type_code, , drop = FALSE]
    if (nrow(s) == 0) {
      # no j-level in the label
      return(sprintf("%s(%d,%d): none", letter, kx_i, ky_i))
    }
    # aggregate by covariate within this tile
    agg <- aggregate(beta ~ covariate, data = s, FUN = sum)
    agg <- agg[order(-abs(agg$beta)), , drop = FALSE]
    agg <- head(agg, top_per_block)
    
    agg$beta_str <- sprintf("%.3f", agg$beta)
    cov_lines <- paste0("  ", agg$covariate, " (", agg$beta_str, ")",
                        collapse = "\n")
    
    # again only (kx, ky)
    sprintf("%s(%d,%d):\n%s", letter, kx_i, ky_i, cov_lines)
  }
  
  
  H_line <- make_line("psi_x",  "H")   # horizontal
  V_line <- make_line("psi_y",  "V")   # vertical
  D_line <- make_line("psi_xy", "D")   # diagonal
  
  labels[i] <- paste(H_line, V_line, D_line, sep = "\n")
}

blocks$label <- labels

## ==========================================================
## 4) Plot: normalized roads + 4×4 j=2 blocks + blue labels
##    for Localized SCAD (BT–Haar)
## ==========================================================

p_blocks_j2_scad <- ggplot(roads_events_01) +
  geom_sf(aes(color = factor(has_incident)),
          linewidth = 0.35,
          alpha = 0.9,
          na.rm = TRUE) +
  scale_color_manual(
    values = c("0" = "grey80", "1" = "red"),
    labels = c("No", "Yes"),
    name   = "Has incident"
  ) +
  coord_sf(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    legend.position  = "right",
    plot.title       = element_text(hjust = 0.5)
  ) +
  labs(
    title =  "Selected variables and estimated regression coefficeint at j=2 by LLS",
    x = NULL, y = NULL
  ) +
  geom_rect(
    data = blocks,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill   = NA,
    colour = "black",
    linewidth = 0.7
  ) +
  geom_label(
    data = blocks,
    aes(x = cx, y = cy, label = label),
    inherit.aes = FALSE,
    size = 2.0,
    fill = "blue",
    colour = "white",
    label.size = 0
  )

print(p_blocks_j2_scad)



add_grid_and_labels <- function(p, roads_sf, n = 2,
                                line_col = "black", line_lwd = 0.4,
                                bubble_fill = "gold",
                                bubble_alpha = 0.9,
                                label_size = 4,
                                pad_frac = 0.08) {
  
  bb <- st_bbox(roads_sf)
  
  xs <- seq(bb["xmin"], bb["xmax"], length.out = n + 1)
  ys <- seq(bb["ymin"], bb["ymax"], length.out = n + 1)
  
  # grid lines
  vlines <- data.frame(x = xs, xend = xs, y = bb["ymin"], yend = bb["ymax"])
  hlines <- data.frame(x = bb["xmin"], xend = bb["xmax"], y = ys, yend = ys)
  
  # cell label positions: TOP-LEFT of each cell (with padding)
  cells <- expand.grid(ix = 1:n, iy = 1:n)
  
  dx <- xs[2] - xs[1]
  dy <- ys[2] - ys[1]
  
  cells$x <- xs[cells$ix] + pad_frac * dx          # a bit right from left edge
  cells$y <- ys[cells$iy + 1] - pad_frac * dy      # a bit down from top edge
  
  # numbering: top-left starts at 1, goes row-wise left->right, top->bottom
  cells$label <- with(cells, (n - iy) * n + ix)
  
  p +
    geom_segment(data = vlines,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, color = line_col, linewidth = line_lwd) +
    geom_segment(data = hlines,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, color = line_col, linewidth = line_lwd) +
    # circle bubble + number (high visibility)
    geom_label(
      data = cells,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      label.size = 0,                 # no border (clean)
      label.r = unit(0.35, "lines"),   # rounded -> circle-like bubble
      fill = bubble_fill,
      alpha = bubble_alpha,
      color = "black",
      size = label_size,
      fontface = "bold"
    )
}

# ---- your base plot (same as before) ----
p_base <- ggplot(roads_events) +
  geom_sf(aes(color = factor(has_incident),
              linewidth = factor(has_incident)),
          na.rm = TRUE) +
  scale_color_manual(values = c("0" = "grey80", "1" = "black"),
                     labels = c("No", "Yes"),
                     name   = "Has incident") +
  scale_linewidth_manual(values = c("0" = 0.3, "1" = 1.2),
                         guide  = "none") +
  theme_minimal() +
  labs(title = "Roads with Any Incident (accident or crime)")

# ---- make 2x2 and 4x4 with visible numbered bubbles ----
p_4  <- add_grid_and_labels(p_base, roads_events, n = 2, bubble_fill = "deepskyblue2") +
  labs(subtitle = "2×2 grid")
p_16 <- add_grid_and_labels(p_base, roads_events, n = 4, bubble_fill = "gold") +
  labs(subtitle = "4×4 grid")

print(p_4)
print(p_16)