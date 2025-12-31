# ============================================================
# COMMENTED SUMMARY (What this whole block does)
# ============================================================
# A) Model-selection helpers (BIC along regularization paths; NO standardization)
#   1) select_by_bic_glmnet():
#      - Fits Poisson LASSO path via glmnet (alpha=1) with standardize=FALSE.
#      - Evaluates BIC for every lambda using your bic_poisson_offset().
#      - Returns beta/lambda achieving minimum BIC + full (lambda, BIC) path.
#   2) select_by_bic_ncvreg():
#      - Fits Poisson SCAD path via ncvreg.
#      - Evaluates BIC for every lambda; returns best beta/lambda + full path.
#
# B) Localized coefficient mapping from vector -> (Haar atom × covariate) matrix
#   3) nonzero_map_local(Gamma, meta, cov_names):
#      - Gamma is R×P where rows correspond to Haar atoms (meta) and columns to covariates.
#      - Finds all nonzero entries (above tol) and returns:
#          (i) selected_covariates (unique covariate names)
#         (ii) selection_by_resolution: detailed table of which atom (j,type,kx,ky, tile bounds)
#              is active for which covariate (interpretable “where/at what scale” selection).
#
# C) Runners: localized vs global variable selection (all using BIC tuning)
#   4) run_local_lasso_bic():
#      - Builds localized design Z = build_localized_Z(Xq, B) of size M × (R·P).
#      - Fits Poisson LASSO (glmnet) with offset off; selects lambda by minimum BIC.
#      - Reshapes coefficient vector into Gamma (R×P) and maps nonzeros via nonzero_map_local().
#
#   5) run_local_scad_bic():
#      - Same as above but uses Poisson SCAD (ncvreg) instead of LASSO.
#
#   6) run_global_lasso_bic() / run_global_scad_bic():
#      - Fits “global” models using ONLY Xq (no Haar localization).
#      - Selects lambda by BIC; returns selected covariates (no resolution table).
#
#   7) run_global_adaptive_lasso_bic():
#      - Computes adaptive weights from an initial ridge (alpha=0, standardize=FALSE),
#        w_j = 1/(|b0_j|+eps)^gamma.
#      - Re-fits Poisson adaptive LASSO via glmnet with penalty.factor = w, lambda chosen by BIC.
#
# D) Run models on a curated covariate subset
#   8) important_covariates:
#      - A hand-curated set of interpretable features: road class, controls/intersections,
#        transit, lighting, parking, buildings/land-use, POIs, rail/waterway, and crime_total_neigh.
#   9) Xq_reduced:
#      - Restricts Xq to columns in important_covariates that actually exist in Xq.
#  10) Fit all 5 models on Xq_reduced:
#      - Localized LASSO (BT–Haar), Localized SCAD (BT–Haar),
#        Global LASSO, Global SCAD, Global Adaptive LASSO.
#  11) Export selections:
#      - lli_sel_cov, lls_sel_cov: detailed (atom-level) selections for localized methods.
#      - lasso_sel_cov, scad_sel_cov, adlasso_sel_cov: selected covariate names for global methods.
#
# E) Haar basis visualization: 64 atoms (Jmax=3) on a single page (8×8 grid)
#  12) enumerate_haar2d(Jmax=3):
#      - Enumerates all 2D Haar atoms up to level Jmax-1:
#          1 scaling atom (LL) + 3 * (sum_{j=0}^{Jmax-1} 4^j) detail atoms = 64 total for Jmax=3.
#      - Stores each atom’s (j, type, kx, ky) and its support tile [x_lo,x_hi)×[y_lo,y_hi) in [0,1]^2.
#  13) draw_atom_panel(...):
#      - Draws one atom: tile outline + +/- sign pattern (H/V/D) or “LL”.
#      - Adds compact label (and optional tile interval text).
#  14) render_8x8(...):
#      - Writes a PDF (“haar64_atoms_8x8.pdf”) and lays out the 64 atoms in an 8×8 grid,
#        with a title and optional support-interval labels if SHOW_TILE=TRUE.
#
# Final outputs produced by this block:
#   - Variable-selection results for 5 models (localized + global) with BIC-tuned lambda.
#   - Interpretable localized selection tables (by Haar resolution/tile) for BT–Haar methods.
#   - A one-page PDF visualization of all 64 2D Haar atoms for Jmax=3.
# ============================================================

## ========================= Robust path selectors (NO STANDARDIZATION) =========================
select_by_bic_glmnet <- function(X, y, offset,
                                 penalty.factor = NULL,
                                 family = "poisson",
                                 n_eff = length(y)) {
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

select_by_bic_ncvreg <- function(X, y, offset, n_eff = length(y)) {
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

## ========================= Mapping & runners =========================
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

run_local_lasso_bic <- function(U, Xq, B, meta, yi, off, cov_names = colnames(Xq),
                                n_eff = length(yi)) {
  Z  <- build_localized_Z(Xq, B)
  cz <- .prep_design(Z)     # no drops, no impute
  Zc <- cz$X
  sel <- select_by_bic_glmnet(Zc, yi, off, penalty.factor = NULL, family = "poisson", n_eff = n_eff)
  R <- ncol(B); P <- ncol(Xq); RP <- R * P
  gamma_hat_full <- numeric(RP); gamma_hat_full[cz$keep] <- sel$beta
  Gamma <- matrix(gamma_hat_full, nrow = R, ncol = P, byrow = FALSE)
  map <- nonzero_map_local(Gamma, meta, cov_names)
  list(
    selected_covariates     = map$selected_covariates,
    selection_by_resolution = map$selection_by_resolution,
    lambda = sel$lambda, bic_min = sel$bic_min, bic_path = sel$bic_path
  )
}

run_local_scad_bic <- function(U, Xq, B, meta, yi, off, cov_names = colnames(Xq),
                               n_eff = length(yi)) {
  Z  <- build_localized_Z(Xq, B)
  cz <- .prep_design(Z)     # no drops, no impute
  Zc <- cz$X
  sel <- select_by_bic_ncvreg(Zc, yi, off, n_eff = n_eff)
  R <- ncol(B); P <- ncol(Xq); RP <- R * P
  gamma_hat_full <- numeric(RP); gamma_hat_full[cz$keep] <- sel$beta
  Gamma <- matrix(gamma_hat_full, nrow = R, ncol = P, byrow = FALSE)
  map <- nonzero_map_local(Gamma, meta, cov_names)
  list(
    selected_covariates     = map$selected_covariates,
    selection_by_resolution = map$selection_by_resolution,
    lambda = sel$lambda, bic_min = sel$bic_min, bic_path = sel$bic_path
  )
}

run_global_adaptive_lasso_bic <- function(
    Xq, yi, off,
    cov_names    = colnames(Xq),
    gamma_adapt  = 1.0,
    eps_w        = 1e-6,
    n_eff        = length(yi)
) {
  cz <- .prep_design(Xq)  # no drops, no impute
  Xc <- cz$X
  # initial ridge (no standardization)
  cv0 <- glmnet::cv.glmnet(
    x = Xc, y = yi, family = "poisson", alpha = 0,
    standardize = FALSE, offset = off       # <-- NO STANDARDIZATION
  )
  b0 <- as.numeric(stats::coef(cv0, s = "lambda.min"))[-1]
  w  <- 1 / (abs(b0) + eps_w)^gamma_adapt
  # adaptive lasso (no standardization)
  sel <- select_by_bic_glmnet(
    X = Xc, y = yi, offset = off,
    penalty.factor = w,
    family = "poisson",
    n_eff = n_eff
  )
  beta_full <- numeric(ncol(Xq)); beta_full[cz$keep] <- sel$beta
  sel_idx <- which(beta_full != 0)
  list(
    selected_covariates = cov_names[sel_idx],
    coefficients = beta_full,
    lambda   = sel$lambda,
    bic_min  = sel$bic_min,
    bic_path = sel$bic_path
  )
}

run_global_lasso_bic <- function(Xq, yi, off, cov_names = colnames(Xq),
                                 n_eff = length(yi)) {
  cg <- .prep_design(Xq)  # no drops, no impute
  Xc <- cg$X
  sel <- select_by_bic_glmnet(Xc, yi, off, penalty.factor = NULL, family = "poisson", n_eff = n_eff)
  beta_full <- numeric(ncol(Xq)); beta_full[cg$keep] <- sel$beta
  nz <- which(abs(beta_full) > 0)
  list(
    selected_covariates     = cov_names[nz],
    selection_by_resolution = data.frame(),
    lambda = sel$lambda, bic_min = sel$bic_min, bic_path = sel$bic_path
  )
}

run_global_scad_bic <- function(Xq, yi, off, cov_names = colnames(Xq),
                                n_eff = length(yi)) {
  cg <- .prep_design(Xq)  # no drops, no impute
  Xc <- cg$X
  sel <- select_by_bic_ncvreg(Xc, yi, off, n_eff = n_eff)
  beta_full <- numeric(ncol(Xq)); beta_full[cg$keep] <- sel$beta
  nz <- which(abs(beta_full) > 0)
  list(
    selected_covariates     = cov_names[nz],
    selection_by_resolution = data.frame(),
    lambda = sel$lambda, bic_min = sel$bic_min, bic_path = sel$bic_path
  )
}

## ========================= Run all models =========================
# Localized models
important_covariates <- c(
  # Road class
  "road_fclass_motorway","road_fclass_motorway_link",
  "road_fclass_trunk","road_fclass_trunk_link",
  "road_fclass_primary","road_fclass_primary_link",
  "road_fclass_secondary","road_fclass_secondary_link",
  "road_fclass_tertiary","road_fclass_tertiary_link",
  "road_fclass_residential","road_fclass_service",
  
  # Intersections / control
  "traffic_motorway_junction","traffic_traffic_signals","traffic_crossing",
  "traffic_stop","traffic_mini_roundabout","traffic_turning_circle",
  
  # Transit
  "trans_bus_stop","trans_tram_stop","trans_railway_station","trans_railway_halt",
  
  # Lighting
  "traffic_street_lamp",
  
  # Parking / access
  "traffic_parking","traffic_parking_multistorey","traffic_parking_bicycle",
  "trafficch_parking","trafficch_parking_multistorey","trafficch_parking_underground",
  
  # Land use / buildings
  "bldg_apartments","bldg_residential","bldg_house","bldg_detached","bldg_terrace",
  "bldg_commercial","bldg_office","bldg_retail","bldg_parking",
  "bldg_school","bldg_college","bldg_university","bldg_hospital",
  "bldg_abandoned","bldg_demolished","bldg_ruins",
  
  # POIs
  "pois_bar","pois_pub","pois_nightclub","pois_beverages",
  "pois_convenience","pois_supermarket","pois_atm","pois_bank",
  "pois_school","pois_kindergarten","pois_university","pois_college",
  "pois_hospital","pois_clinic","pois_courthouse",
  "pois_park","pois_playground","pois_stadium",
  
  # Rail / waterway context
  "rail_rail","rail_tram","rail_light_rail",
  "waterway_river","waterway_canal",
  
  # Neighborhood baseline
  "crime_total_neigh"
)

keep_now <- intersect(important_covariates, colnames(Xq))
Xq_reduced <- as.matrix(Xq[, keep_now, drop = FALSE])

# then rerun your models with Xq_reduced:
out_loc_lasso_imp <- run_local_lasso_bic(U, Xq_reduced, B, meta, yi, off, cov_names = colnames(Xq_reduced))
out_loc_scad_imp  <- run_local_scad_bic (U, Xq_reduced, B, meta, yi, off, cov_names = colnames(Xq_reduced))
out_gl_lasso_imp  <- run_global_lasso_bic(Xq_reduced, yi, off, cov_names = colnames(Xq_reduced))
out_gl_scad_imp   <- run_global_scad_bic (Xq_reduced, yi, off, cov_names = colnames(Xq_reduced))
out_global_adalasso_imp <- run_global_adaptive_lasso_bic(
  Xq = Xq_reduced, yi = yi, off = off,
  cov_names = colnames(Xq_reduced),
  gamma_adapt = 1.0, eps_w = 1e-6
)


## ========================= Expose selections =========================
lli_sel_cov      <- out_loc_lasso_imp$selection_by_resolution
lls_sel_cov      <- out_loc_scad_imp$selection_by_resolution
adlasso_sel_cov  <- out_global_adalasso_imp$selected_covariates
lasso_sel_cov    <- out_gl_lasso_imp$selected_covariates
scad_sel_cov     <- out_gl_scad_imp$selected_covariates


# ==========================================================
# 64 Haar atoms (Jmax = 3) on ONE page in an 8×8 grid
# ----------------------------------------------------------
# Output: haar64_atoms_8x8.pdf
# Each cell shows: support tile + Haar +/- pattern
# Label: LL or H/V/D(kx,ky)  j=...
# Optionally show tile intervals (set SHOW_TILE = TRUE)
# ==========================================================

library(grid)

# ---------- Enumerate full 2D Haar basis up to Jmax ----------
enumerate_haar2d <- function(Jmax = 3) {
  stopifnot(Jmax >= 1)
  rows <- list()
  # Scaling (global)
  rows[[length(rows)+1]] <- data.frame(
    name="LL", type="phi", j=0, kx=0, ky=0,
    x_lo=0, x_hi=1, y_lo=0, y_hi=1, stringsAsFactors = FALSE
  )
  # Details j = 0..Jmax-1
  for (j in 0:(Jmax-1)) {
    K <- 2^j
    for (kx in 0:(K-1)) for (ky in 0:(K-1)) {
      x_lo <- kx/K; x_hi <- (kx+1)/K
      y_lo <- ky/K; y_hi <- (ky+1)/K
      rows[[length(rows)+1]] <- data.frame(name="H", type="psi_x", j=j,kx=kx,ky=ky,
                                           x_lo=x_lo,x_hi=x_hi,y_lo=y_lo,y_hi=y_hi)
      rows[[length(rows)+1]] <- data.frame(name="V", type="psi_y", j=j,kx=kx,ky=ky,
                                           x_lo=x_lo,x_hi=x_hi,y_lo=y_lo,y_hi=y_hi)
      rows[[length(rows)+1]] <- data.frame(name="D", type="psi_xy", j=j,kx=kx,ky=ky,
                                           x_lo=x_lo,x_hi=x_hi,y_lo=y_lo,y_hi=y_hi)
    }
  }
  df <- do.call(rbind, rows)
  # Order: LL, then by j (0,1,2), then ky, kx, and H,V,D
  ord <- with(df, order(type != "phi", j, ky, kx, match(name, c("H","V","D"))))
  df[ord, , drop=FALSE]
}

# ---------- Draw a single atom panel ----------
draw_atom_panel <- function(xlo, xhi, ylo, yhi, type,
                            label_main, label_tile = NULL,
                            base_cex = 0.7) {
  # frame
  grid.rect(gp=gpar(fill=NA, col="grey40", lwd=0.6))
  # support tile
  grid.rect(x=unit(xlo,"npc"), y=unit(ylo,"npc"),
            width=unit(xhi-xlo,"npc"), height=unit(yhi-ylo,"npc"),
            just=c("left","bottom"),
            gp=gpar(fill="grey85", alpha=0.9, col="black", lwd=0.7))
  # sign pattern
  xm <- (xlo+xhi)/2; ym <- (ylo+yhi)/2
  if (type == "psi_x") {
    grid.lines(x=unit(c(xm,xm),"npc"), y=unit(c(ylo,yhi),"npc"),
               gp=gpar(col="grey20", lwd=0.8))
    grid.text("-", x=unit((xlo+xm)/2,"npc"), y=unit(ym,"npc"), gp=gpar(cex=base_cex))
    grid.text("+", x=unit((xm+xhi)/2,"npc"), y=unit(ym,"npc"), gp=gpar(cex=base_cex))
  } else if (type == "psi_y") {
    grid.lines(x=unit(c(xlo,xhi),"npc"), y=unit(c(ym,ym),"npc"),
               gp=gpar(col="grey20", lwd=0.8))
    grid.text("-", x=unit(xm,"npc"), y=unit((ylo+ym)/2,"npc"), gp=gpar(cex=base_cex))
    grid.text("+", x=unit(xm,"npc"), y=unit((ym+yhi)/2,"npc"), gp=gpar(cex=base_cex))
  } else if (type == "psi_xy") {
    grid.lines(x=unit(c(xm,xm),"npc"), y=unit(c(ylo,yhi),"npc"),
               gp=gpar(col="grey20", lwd=0.8))
    grid.lines(x=unit(c(xlo,xhi),"npc"), y=unit(c(ym,ym),"npc"),
               gp=gpar(col="grey20", lwd=0.8))
    grid.text("-", x=unit((xlo+xm)/2,"npc"), y=unit((ylo+ym)/2,"npc"), gp=gpar(cex=base_cex))
    grid.text("+", x=unit((xm+xhi)/2,"npc"), y=unit((ylo+ym)/2,"npc"), gp=gpar(cex=base_cex))
    grid.text("+", x=unit((xlo+xm)/2,"npc"), y=unit((ym+yhi)/2,"npc"), gp=gpar(cex=base_cex))
    grid.text("-", x=unit((xm+xhi)/2,"npc"), y=unit((ym+yhi)/2,"npc"), gp=gpar(cex=base_cex))
  } else {
    grid.text("LL", x=0.5, y=0.5, gp=gpar(cex=base_cex+0.1))
  }
  # label inside (top-left) to save space
  grid.text(label_main, x=unit(0.02,"npc"), y=unit(0.98,"npc"),
            just=c("left","top"), gp=gpar(cex=base_cex*0.9))
  if (!is.null(label_tile) && nzchar(label_tile))
    grid.text(label_tile, x=unit(0.02,"npc"), y=unit(0.80,"npc"),
              just=c("left","top"), gp=gpar(cex=base_cex*0.75))
}

# ---------- One-page 8×8 renderer ----------
render_8x8 <- function(df, file="haar64_atoms_8x8.pdf",
                       SHOW_TILE = FALSE,  # set TRUE to show intervals
                       page_width_in = 11, page_height_in = 11) {
  
  stopifnot(nrow(df) == 64)
  pdf(file, width = page_width_in, height = page_height_in)
  on.exit(dev.off(), add=TRUE)
  
  grid.newpage()
  # Title
  grid.text("Full 2D Haar Basis (Jmax = 3) — 64 Atoms (8×8)",
            y=unit(1, "npc") - unit(0.02, "npc"),
            gp=gpar(fontface="bold", cex=1.1))
  
  # Margins
  top_gap <- 0.06; bottom_gap <- 0.03; left_gap <- 0.03; right_gap <- 0.02
  pushViewport(viewport(x=unit(left_gap,"npc"), y=unit(bottom_gap,"npc"),
                        width=unit(1-left_gap-right_gap,"npc"),
                        height=unit(1-top_gap-bottom_gap,"npc"),
                        just=c("left","bottom")))
  
  lay <- grid.layout(nrow=8, ncol=8,
                     widths=unit(rep(1,8),"null"),
                     heights=unit(rep(1,8),"null"))
  pushViewport(viewport(layout=lay))
  
  idx <- 1
  for (r in 1:8) {
    for (c in 1:8) {
      pushViewport(viewport(layout.pos.row=r, layout.pos.col=c))
      # Compose compact label
      main_lab <- if (df$name[idx] == "LL") {
        "LL"
      } else {
        sprintf("%s(%d,%d)  j=%d", df$name[idx], df$kx[idx], df$ky[idx], df$j[idx])
      }
      tile_lab <- if (SHOW_TILE) {
        sprintf("[%.2f,%.2f)×[%.2f,%.2f)", df$x_lo[idx], df$x_hi[idx],
                df$y_lo[idx], df$y_hi[idx])
      } else ""
      draw_atom_panel(df$x_lo[idx], df$x_hi[idx], df$y_lo[idx], df$y_hi[idx],
                      df$type[idx], main_lab, tile_lab, base_cex = 0.70)
      popViewport()
      idx <- idx + 1
    }
  }
  
  popViewport(); popViewport()
  message("Saved: ", normalizePath(file))
}

# ---------- Build and render ----------
atoms <- enumerate_haar2d(Jmax = 3)  # 64 atoms
# Single-page 8×8
render_8x8(atoms, file = "haar64_atoms_8x8.pdf",
           SHOW_TILE = FALSE,           # switch to TRUE if you want bounds in each cell
           page_width_in = 12, page_height_in = 12)

