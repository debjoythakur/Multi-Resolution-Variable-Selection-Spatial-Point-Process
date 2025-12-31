# ============================================================
# SUMMARY / COMMENTARY (What this entire code block does)
# ============================================================
# This script has TWO main parts:
#
# PART A) Incident visualization on the road network (ggplot2)
#   - Plot road segments colored/thickened by whether the road has any incident
#     (accident or crime) attached to it (has_incident = 1).
#
# PART B) Haar basis visualization + Localized BT–Haar Poisson LASSO (spatial point-process style)
#   - Draw a grayscale visualization of the 2D Haar basis (J=2) and export to PDF.
#   - Build a point-pattern quadrature scheme over the spatial domain.
#   - Build a 2D Haar basis evaluated on quadrature points.
#   - Convert road-level covariates into smoothed spatial “images” on the SAME quadrature grid.
#   - Form a localized design matrix Z = (Haar atoms) × (covariate images) and fit a Poisson
#     model with LASSO-type penalties (framework prepared here; fitting occurs later).
#
# ============================================================
# PART A: Plot roads with incident highlighting (ggplot2 + sf)
# ============================================================
# - roads_events is an sf LINE layer with per-road attributes:
#     has_incident, acc_count_road, crime_count_road, killed_total, injured_total, etc.
# - The map draws:
#     has_incident=0 → grey thin lines
#     has_incident=1 → black thick lines
# - linewidth is used visually but its legend is hidden for cleaner plot.
#
# Output:
#   A simple map to verify snapping/aggregation worked and to visually inspect hotspots.
#
# ============================================================
# PART B1: Haar basis figure export (grid graphics)
# ============================================================
# 1) draw_haar_atom_grey():
#    - Uses grid to draw one Haar atom (phi, psi_x, psi_y, psi_xy) using grayscale tiles.
#    - Shows +/- regions via darker/lighter greys and optionally outlines the support.
#
# 2) render_haar2d_J2_greyscale():
#    - Creates a PDF panel figure for the full 2D Haar basis at J=2:
#        Top: j=0 → 1 scaling (LL) + 3 global details (H, V, D)
#        Bottom: j=1 → 12 localized details (H/V/D in 4 quadrants)
#    - Adds a controlled vertical GAP and separator line between j=0 and j=1 blocks.
#
# Output:
#   "haar2d_J2_greyscale_gap.pdf" (a publication-quality Haar basis illustration).
#
# ============================================================
# PART B2: Prepare modeling data (subset to roads with events)
# ============================================================
# - model_df: road-level modeling dataframe already assembled upstream.
# - temp_event = crime_count_road + acc_count_road
# - model_df_temp = rows with temp_event > 0 (i.e., incident-positive roads only)
#
# Then:
# - s_n: spatial coordinates (columns 2 and 3 in model_df_temp; must be planar x/y)
# - X  : covariate matrix (drops selected columns by index; keeps all others)
#
# Notes:
# - The code later enforces: NO IMPUTATION / NO STANDARDIZATION / NO DROPPING.
# - Any non-finite value in X triggers a hard stop (by design).
#
# ============================================================
# PART B3: Build a quadrature scheme (spatstat) for Poisson intensity modeling
# ============================================================
# Purpose:
#   A Poisson point-process likelihood is approximated using Berman–Turner quadrature.
#
# Steps:
# - Create an observation window W covering the point locations (with small padding).
# - Create a ppp object Xpp from the coordinates s_n.
# - Construct a grid-based quadrature scheme Q via quadscheme(..., ntile = c(16,16)).
# - Extract, in a consistent internal ordering:
#     U  : dataframe of quadrature point coordinates (data + dummy points)
#     yi : 1 if quadrature point is a data point, 0 if dummy
#     wi : quadrature weights (repaired to be positive finite if needed)
# - Build offset off = log(wi) which is used in the BT log-likelihood.
#
# Key objects produced:
#   W, Xpp, Q, U (x,y), yi (0/1), wi (>0), off = log(wi).
#
# ============================================================
# PART B4: Construct 2D Haar basis on quadrature points
# ============================================================
# Steps:
# - normalize_to_unit_square(U): rescales coordinates to [0,1]^2 for Haar evaluation.
# - build_haar2d_meta(..., J=3):
#     * constructs scaling atoms + wavelet atoms (psi_x, psi_y, psi_xy) up to level J-1
#     * returns:
#         B    : (M x R) matrix, Haar basis evaluated at all quadrature points (M = nrow(U))
#         meta : dataframe describing each atom (j, type, kx, ky, support box, etc.)
#
# Key objects:
#   B (Haar basis matrix), meta (atom descriptors for interpretability/plotting).
#
# ============================================================
# PART B5: Convert road covariates into spatial images on the same grid
# ============================================================
# Purpose:
#   We need covariate values at quadrature points U, not only at observed road points.
#   The script creates smooth spatial surfaces (images) for each covariate.
#
# Steps (for each covariate column Xdf[[j]]):
# - Assign covariate values as marks on Xpp (the observed points).
# - Smooth the marked point pattern into a pixel image using Smooth(..., sigma=0.1, dimyx=qtiles).
# - Evaluate the image at each quadrature point (nearest pixel) → gives a value at U.
#
# Result:
# - Xq: (M x P) matrix of covariate values aligned with quadrature points.
#
# Important:
# - No missing/non-finite marks are allowed; the code stops if they appear.
#
# ============================================================
# PART B6: Build localized design Z = B ⊙ Xq (tensor-style interaction)
# ============================================================
# build_localized_Z(Xq, B):
# - For each covariate p, multiply that covariate surface Xq[,p] by every Haar atom column in B.
# - Stack blocks horizontally:
#     Z = [ B*Xq[,1] | B*Xq[,2] | ... | B*Xq[,P] ]
#
# Dimensions:
#   If B is (M x R) and Xq is (M x P), then Z is (M x (P*R)).
#
# Interpretation:
#   Each coefficient corresponds to a covariate effect localized to a Haar tile/scale.
#
# ============================================================
# PART B7: Likelihood and BIC helpers (Poisson with offset)
# ============================================================
# - loglik_poisson_offset(): computes sum( y*eta - exp(eta) ) optionally weighted.
# - bic_poisson_offset(): computes -2*logLik + (#nonzero coefficients)*log(n_eff).
#
# These are utilities for model selection/penalty tuning after fitting.
#
# ============================================================
# PART B8: Design preparation stub (explicitly “do nothing”)
# ============================================================
# .prep_design():
# - Intentionally does not impute, standardize, drop NA-only columns, or drop constants.
# - Returns X unchanged (as matrix) and bookkeeping fields.
#
# Why:
#   You requested a pipeline with NO automatic cleaning; errors should surface upstream.
#
# ============================================================
# Final deliverables from this block (ready for fitting step)
# ============================================================
# - A road map plot (incident vs no-incident)
# - A Haar basis PDF illustration
# - Quadrature objects: U, yi, wi, off
# - Haar basis: B, meta
# - Quadrature-aligned covariates: Xq
# - Localized design matrix: Z
# - Likelihood/BIC helper functions
#
# Next step (not shown here):
# - Fit penalized Poisson model using glmnet / ncvreg on:
#     response: yi
#     design  : Z (and possibly intercept)
#     offset  : off
#     weights : wi (optional; depends on your chosen BT formulation)
# ============================================================


library(ggplot2)

ggplot(roads_events) +
  geom_sf(
    aes(color = factor(has_incident),
        linewidth = factor(has_incident)),
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c("0" = "grey80", "1" = "black"),
    labels = c("No", "Yes"),
    name   = "Has incident"
  ) +
  scale_linewidth_manual(
    values = c("0" = 0.3, "1" = 1.2),   # thin vs bold
    guide  = "none"                     # hide linewidth legend
  ) +
  theme_minimal() +
  labs(title = "Roads with Any Incident (accident or crime)")

library(grid)

# ---- draw a single atom in grayscale (no inside "H/V/D" label) ----
draw_haar_atom_grey <- function(xlo, xhi, ylo, yhi, type,
                                show_bounds = FALSE,
                                base_cex = 0.85) {
  
  # panel frame
  grid.rect(gp=gpar(fill=NA, col="grey35", lwd=0.7))
  
  # optional support tile outline
  if (show_bounds) {
    grid.rect(x=unit(xlo,"npc"), y=unit(ylo,"npc"),
              width=unit(xhi-xlo,"npc"), height=unit(yhi-ylo,"npc"),
              just=c("left","bottom"),
              gp=gpar(fill="grey92", col="grey20", lwd=0.7))
  }
  
  xm <- (xlo + xhi)/2
  ym <- (ylo + yhi)/2
  
  subrect <- function(ax0, ax1, ay0, ay1, fill_grey) {
    grid.rect(x=unit(ax0,"npc"), y=unit(ay0,"npc"),
              width=unit(ax1-ax0,"npc"), height=unit(ay1-ay0,"npc"),
              just=c("left","bottom"),
              gp=gpar(fill=fill_grey, col=NA))
  }
  
  # grayscale for signs
  g_plus  <- "grey25"  # dark
  g_minus <- "grey82"  # light
  
  if (type == "phi") {
    subrect(xlo, xhi, ylo, yhi, "grey55")
    grid.text("LL", x=0.5, y=0.5, gp=gpar(cex=base_cex+0.1, col="white"))
  }
  
  if (type == "psi_x") {   # H
    subrect(xlo, xm, ylo, yhi, g_plus)
    subrect(xm, xhi, ylo, yhi, g_minus)
    grid.lines(x=unit(c(xm,xm),"npc"), y=unit(c(ylo,yhi),"npc"),
               gp=gpar(col="grey30", lwd=0.8))
    grid.text("+", x=unit((xlo+xm)/2,"npc"), y=unit((ylo+yhi)/2,"npc"),
              gp=gpar(cex=base_cex, col="white"))
    grid.text("-", x=unit((xm+xhi)/2,"npc"), y=unit((ylo+yhi)/2,"npc"),
              gp=gpar(cex=base_cex, col="grey10"))
  }
  
  if (type == "psi_y") {   # V
    subrect(xlo, xhi, ylo, ym, g_plus)
    subrect(xlo, xhi, ym, yhi, g_minus)
    grid.lines(x=unit(c(xlo,xhi),"npc"), y=unit(c(ym,ym),"npc"),
               gp=gpar(col="grey30", lwd=0.8))
    grid.text("+", x=unit((xlo+xhi)/2,"npc"), y=unit((ylo+ym)/2,"npc"),
              gp=gpar(cex=base_cex, col="white"))
    grid.text("-", x=unit((xlo+xhi)/2,"npc"), y=unit((ym+yhi)/2,"npc"),
              gp=gpar(cex=base_cex, col="grey10"))
  }
  
  if (type == "psi_xy") {  # D
    subrect(xlo, xm, ylo, ym, g_plus)
    subrect(xm, xhi, ylo, ym, g_minus)
    subrect(xlo, xm, ym, yhi, g_minus)
    subrect(xm, xhi, ym, yhi, g_plus)
    grid.lines(x=unit(c(xm,xm),"npc"), y=unit(c(ylo,yhi),"npc"),
               gp=gpar(col="grey30", lwd=0.8))
    grid.lines(x=unit(c(xlo,xhi),"npc"), y=unit(c(ym,ym),"npc"),
               gp=gpar(col="grey30", lwd=0.8))
    grid.text("+", x=unit((xlo+xm)/2,"npc"), y=unit((ylo+ym)/2,"npc"),
              gp=gpar(cex=base_cex, col="white"))
    grid.text("-", x=unit((xm+xhi)/2,"npc"), y=unit((ylo+ym)/2,"npc"),
              gp=gpar(cex=base_cex, col="grey10"))
    grid.text("-", x=unit((xlo+xm)/2,"npc"), y=unit((ym+yhi)/2,"npc"),
              gp=gpar(cex=base_cex, col="grey10"))
    grid.text("+", x=unit((xm+xhi)/2,"npc"), y=unit((ym+yhi)/2,"npc"),
              gp=gpar(cex=base_cex, col="white"))
  }
}


# --- (keep your draw_haar_atom_grey() as-is) ---

render_haar2d_J2_greyscale <- function(file="haar2d_J2_greyscale.pdf",
                                       width_in=11, height_in=9,
                                       show_bounds=FALSE,
                                       gap_mm = 6,      # <-- space between j=0 and j=1
                                       line_mm = 1.2) { # <-- thickness of separator row
  
  pdf(file, width=width_in, height=height_in)
  on.exit(dev.off(), add=TRUE)
  
  grid.newpage()
  
  # Titles
  grid.text("Full 2D Haar Basis for J = 2",
            x=0.02, y=0.98, just=c("left","top"),
            gp=gpar(fontface="bold", cex=1.25))
  grid.text("Top: j = 0 (1 scaling + 3 details).   Bottom: j = 1 (12 localized details).",
            x=0.02, y=0.945, just=c("left","top"),
            gp=gpar(cex=0.95, col="grey20"))
  
  # Outer area
  pushViewport(viewport(x=0.04, y=0.06, width=0.92, height=0.84,
                        just=c("left","bottom")))
  
  # >>> KEY CHANGE: 4-row layout (top, gap, line, bottom)
  lay <- grid.layout(
    nrow = 4, ncol = 1,
    heights = unit.c(
      unit(0.40, "npc"),    # j=0 block height
      unit(gap_mm, "mm"),   # GAP
      unit(line_mm, "mm"),  # separator line row
      unit(1, "null")       # j=1 block (takes remaining)
    )
  )
  pushViewport(viewport(layout=lay))
  
  # ---------------- Row 1: TOP (j=0) ----------------
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
  grid.text("j = 0 (Scaling + Global Details)",
            x=0, y=1.02, just=c("left","bottom"),
            gp=gpar(fontface="bold", cex=1.05))
  
  lay_top <- grid.layout(nrow=1, ncol=4, widths=unit(rep(1,4),"null"))
  pushViewport(viewport(layout=lay_top))
  
  top_names <- c("LL","H","V","D")
  top_types <- c("phi","psi_x","psi_y","psi_xy")
  strip_h <- 0.13
  
  for (cc in 1:4) {
    pushViewport(viewport(layout.pos.row=1, layout.pos.col=cc))
    grid.rect(x=0.5, y=1-strip_h/2, width=1, height=strip_h,
              gp=gpar(fill="grey85", col="grey50", lwd=0.7))
    grid.text(top_names[cc], x=0.5, y=1-strip_h/2,
              gp=gpar(cex=0.95, col="grey10"))
    
    pushViewport(viewport(y=0, height=1-strip_h, just=c("center","bottom")))
    draw_haar_atom_grey(0,1,0,1, top_types[cc],
                        show_bounds=show_bounds, base_cex=0.9)
    popViewport()
    popViewport()
  }
  popViewport(); popViewport()
  
  # ---------------- Row 2: GAP ----------------
  # (Nothing drawn; this row is pure whitespace)
  
  # ---------------- Row 3: SEPARATOR LINE ----------------
  pushViewport(viewport(layout.pos.row=3, layout.pos.col=1))
  grid.lines(x=unit(c(0,1),"npc"), y=unit(c(0.5,0.5),"npc"),
             gp=gpar(col="grey40", lwd=1))
  popViewport()
  
  # ---------------- Row 4: BOTTOM (j=1) ----------------
  pushViewport(viewport(layout.pos.row=4, layout.pos.col=1))
  grid.text("j = 1 (Detail Atoms in Quadrants)",
            x=0, y=1.02, just=c("left","bottom"),
            gp=gpar(fontface="bold", cex=1.05))
  
  lay_bot <- grid.layout(
    nrow=3, ncol=5,
    widths=unit.c(unit(0.10,"npc"), unit(rep(0.225,4),"npc")),
    heights=unit(rep(1,3),"null")
  )
  pushViewport(viewport(layout=lay_bot))
  
  row_lab <- c("H","V","D")
  row_typ <- c("psi_x","psi_y","psi_xy")
  quads <- data.frame(kx=c(0,0,1,1), ky=c(0,1,0,1))
  strip_h2 <- 0.14
  
  for (r in 1:3) {
    pushViewport(viewport(layout.pos.row=r, layout.pos.col=1))
    grid.text(row_lab[r], x=0.5, y=0.5, gp=gpar(fontface="bold", cex=1.05))
    popViewport()
    
    for (c in 1:4) {
      kx <- quads$kx[c]; ky <- quads$ky[c]
      K <- 2
      xlo <- kx/K; xhi <- (kx+1)/K
      ylo <- ky/K; yhi <- (ky+1)/K
      
      pushViewport(viewport(layout.pos.row=r, layout.pos.col=c+1))
      
      grid.rect(x=0.5, y=1-strip_h2/2, width=1, height=strip_h2,
                gp=gpar(fill="grey85", col="grey50", lwd=0.7))
      grid.text(sprintf("(k1=%d, k2=%d)", kx, ky),
                x=0.5, y=1-strip_h2/2,
                gp=gpar(cex=0.85, col="grey10"))
      
      pushViewport(viewport(y=0, height=1-strip_h2, just=c("center","bottom")))
      draw_haar_atom_grey(xlo,xhi,ylo,yhi, row_typ[r],
                          show_bounds=show_bounds, base_cex=0.85)
      popViewport()
      
      popViewport()
    }
  }
  
  popViewport(); popViewport()
  popViewport(); popViewport()
  
  message("Saved: ", normalizePath(file))
}

render_haar2d_J2_greyscale(
  file = "haar2d_J2_greyscale_gap.pdf",
  gap_mm = 10,
  line_mm = 1.5
)

temp_event <- model_df$crime_count_road + model_df$acc_count_road
model_df$temp_event <- temp_event
temp_kill <- model_df$killed_total
temp_injured <- model_df$injured_total


model_df_temp <- model_df[model_df$temp_event > 0, ]

names(model_df_temp)
s_n <- model_df_temp[,c(2,3)]; X<- model_df_temp[,-c(1,2,3,5,10,276)]
## ========================= Localized BT–Haar Poisson LASSO (with model_df_temp) =========================
## ========================= Localized BT–Haar Poisson LASSO (with model_df_temp) =========================
## ========================= Localized BT–Haar Poisson LASSO (NO IMPUTATION / NO STANDARDIZATION / NO DROPPING) =========================
suppressPackageStartupMessages({
  library(spatstat)       # metapkg (geom/core/explore)
  library(glmnet)         # LASSO / ridge / adaptive LASSO
  library(ncvreg)         # SCAD
  library(dplyr)
})

## ------------------------- Inputs from model_df_temp -------------------------
stopifnot(exists("model_df_temp"))

# EXACT selection per your request
s_n <- as.matrix(model_df_temp[, c(2, 3)])                      # coordinates (x,y)
X   <- model_df_temp[, -c(1, 2, 3, 5, 10, 276), drop = FALSE]   # covariates

# Coerce to numeric ONLY (no imputation), keep names
X <- as.data.frame(X)
if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
for (j in seq_along(X)) {
  if (!is.numeric(X[[j]])) X[[j]] <- suppressWarnings(as.numeric(X[[j]]))
}

# Hard checks: we do not impute; require finite values
if (any(!is.finite(as.matrix(X)))) {
  bad_cols <- names(which(colSums(!is.finite(as.matrix(X))) > 0))
  stop(sprintf(
    "Non-finite values detected in X for columns: %s.\nThis script does not impute by design. Please clean those values upstream.",
    paste(bad_cols, collapse = ", ")
  ))
}
stopifnot(nrow(s_n) == nrow(X))
Xdf <- X

## ========================= WINDOW & PPP =========================
make_padded_window <- function(xr, yr) {
  dx <- diff(xr); dy <- diff(yr)
  if (dx <= 0) xr <- xr + c(-1e-3, 1e-3)
  if (dy <= 0) yr <- yr + c(-1e-3, 1e-3)
  spatstat.geom::owin(xrange = xr, yrange = yr)
}
W   <- make_padded_window(range(s_n[,1]), range(s_n[,2]))
Xpp <- spatstat.geom::ppp(s_n[,1], s_n[,2], window = W, checkdup = FALSE)

## ========================= QUADRATURE (single source of truth) =========================
qtiles <- c(16, 16)  # grid resolution for BT quadrature
Q  <- quadscheme(Xpp, method = "grid", ntile = qtiles)

# 1) data/dummy indicator in internal order
isD_int <- as.logical(is.data(Q))

# 2) quadrature weights in internal order
wi <- try(as.numeric(weights(Q)), silent = TRUE)
if (inherits(wi, "try-error") || length(wi) == 0L) {
  wi <- if (!is.null(Q$w)) as.numeric(Q$w) else numeric(0)
}
if (length(wi) == 0L) {
  if (exists("w.quad", where = asNamespace("spatstat.model"), inherits = FALSE)) {
    wi <- as.numeric(getFromNamespace("w.quad", "spatstat.model")(Q))
  }
}

# 3) coordinates of quadrature points in *internal order*
Upp_try <- try(spatstat.geom::as.ppp(Q), silent = TRUE)
if (!inherits(Upp_try, "try-error")) {
  U  <- data.frame(x = Upp_try$x, y = Upp_try$y)
  yi <- as.integer(isD_int)
} else {
  # fallback: (data, dummy) order + reorder weights
  Xd  <- if (!is.null(Q$data))  Q$data  else spatstat.model::data.quad(Q)
  Xdm <- if (!is.null(Q$dummy)) Q$dummy else spatstat.model::dummy.quad(Q)
  U   <- data.frame(x = c(Xd$x, Xdm$x), y = c(Xd$y, Xdm$y))
  yi  <- c(rep(1L, spatstat.geom::npoints(Xd)),
           rep(0L, spatstat.geom::npoints(Xdm)))
  if (length(wi) > 0L) {
    idx_data  <- which(isD_int)
    idx_dummy <- which(!isD_int)
    wi <- c(wi[idx_data], wi[idx_dummy])
  }
}

# 4) validate/repair weights (no imputation into X; but BT needs positive weights)
M <- nrow(U)
Wow <- spatstat.geom::as.owin(W)
cell_area <- (diff(Wow$xrange)/qtiles[1]) * (diff(Wow$yrange)/qtiles[2])
if (length(wi) != M) wi <- rep(NA_real_, M)
bad <- !is.finite(wi) | wi <= 0
if (all(bad)) {
  wi <- ifelse(yi == 1L, 0, cell_area)
} else {
  wi[bad] <- cell_area
}
off <- log(pmax(wi, 1e-12))

stopifnot(length(yi) == M, length(wi) == M, nrow(U) == M, all(is.finite(off)))

## ========================= HAAR BASIS (on quadrature points U) =========================
normalize_to_unit_square <- function(s) {
  s <- as.matrix(s)
  r1 <- range(s[,1]); r2 <- range(s[,2])
  cbind((s[,1]-r1[1])/max(r1[2]-r1[1], 1e-12),
        (s[,2]-r2[1])/max(r2[2]-r2[1], 1e-12))
}
haar_phi_1d <- function(j,k,x){ a <- k/2^j; b <- (k+1)/2^j; as.numeric(x>=a & x<b)*2^(j/2) }
haar_psi_1d <- function(j,k,x){
  a <- k/2^j; m <- (k+0.5)/2^j; b <- (k+1)/2^j
  (as.numeric(x>=a & x<m) - as.numeric(x>=m & x<b)) * 2^(j/2)
}
build_haar2d_meta <- function(s01, j0=0, J=2, original_xy = NULL){
  x <- s01[,1]; y <- s01[,2]
  Phi <- list(); meta <- list(); cnt <- 0L
  n_k <- max(2^j0, 1)
  # scaling (phi)
  for (kx in 0:(n_k-1)) for (ky in 0:(n_k-1)) {
    cnt <- cnt+1L
    Phi[[cnt]] <- haar_phi_1d(j0,kx,x)*haar_phi_1d(j0,ky,y)
    meta[[cnt]] <- list(j=j0, type="phi", kx=kx, ky=ky)
  }
  # wavelets
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
  meta_df$atom_index <- seq_len(nrow(meta_df))
  meta_df$x_lo01 <- meta_df$kx / (2^meta_df$j)
  meta_df$x_hi01 <- (meta_df$kx+1) / (2^meta_df$j)
  meta_df$y_lo01 <- meta_df$ky / (2^meta_df$j)
  meta_df$y_hi01 <- (meta_df$ky+1) / (2^meta_df$j)
  if (!is.null(original_xy)) {
    xr <- range(original_xy[,1]); yr <- range(original_xy[,2])
    mapx <- function(u01) xr[1] + u01 * (xr[2]-xr[1])
    mapy <- function(u01) yr[1] + u01 * (yr[2]-yr[1])
    meta_df$x_lo <- mapx(meta_df$x_lo01); meta_df$x_hi <- mapx(meta_df$x_hi01)
    meta_df$y_lo <- mapy(meta_df$y_lo01); meta_df$y_hi <- mapy(meta_df$y_hi01)
  } else {
    meta_df$x_lo <- NA_real_; meta_df$x_hi <- NA_real_
    meta_df$y_lo <- NA_real_; meta_df$y_hi <- NA_real_
  }
  meta_df$desc <- sprintf("j=%d, %s, tile(kx=%d,ky=%d), x∈[%.3f,%.3f), y∈[%.3f,%.3f)",
                          meta_df$j, meta_df$type, meta_df$kx, meta_df$ky,
                          meta_df$x_lo01, meta_df$x_hi01, meta_df$y_lo01, meta_df$y_hi01)
  list(B=B, meta=meta_df)
}

s01      <- normalize_to_unit_square(U)
haar_out <- build_haar2d_meta(s01, j0 = 0, J = 3, original_xy = as.matrix(U))
B        <- haar_out$B
meta     <- haar_out$meta
stopifnot(nrow(B) == M, nrow(meta) == ncol(B))

## ========================= Covariate images on the SAME quadrature grid =========================
# No imputation: require finite marks (spatstat smoothing expects numeric marks)
if (any(!is.finite(as.matrix(Xdf)))) {
  bad_cols <- names(which(colSums(!is.finite(as.matrix(Xdf))) > 0))
  stop(sprintf(
    "Non-finite values detected in X (prior to image step) for columns: %s.\nThis script does not impute by design. Please clean those values upstream.",
    paste(bad_cols, collapse = ", ")
  ))
}

eval_cov_im_at <- function(im, x, y) {
  np  <- spatstat.geom::nearest.pixel(x, y, im)
  col <- pmin(pmax(np$col, 1L), ncol(im$v))
  row <- pmin(pmax(np$row, 1L), nrow(im$v))
  im$v[cbind(row, col)]
}

cov_images <- vector("list", ncol(Xdf)); names(cov_images) <- colnames(Xdf)
for (j in seq_len(ncol(Xdf))) {
  v <- Xdf[[j]]
  if (!is.numeric(v)) v <- suppressWarnings(as.numeric(v))
  if (any(!is.finite(v))) {
    stop(sprintf("Column '%s' has non-finite marks; please clean upstream (no imputation performed here).", names(Xdf)[j]))
  }
  Xm <- Xpp; spatstat.geom::marks(Xm) <- v
  imj <- spatstat.explore::Smooth(Xm, sigma = 0.1, at = "pixels",
                                  normalise = TRUE, dimyx = qtiles)
  if (!all(is.finite(imj$v))) {
    stop(sprintf("Smoothed image for '%s' produced non-finite values; adjust sigma/dimyx or clean input.", names(Xdf)[j]))
  }
  cov_images[[j]] <- imj
}
cov_list <- lapply(cov_images, function(im) eval_cov_im_at(im, U$x, U$y))
Xq <- as.matrix(as.data.frame(cov_list))
colnames(Xq) <- colnames(Xdf)

## ========================= Localized design Z and alignment =========================
build_localized_Z <- function(X, B) {
  X <- as.matrix(X); B <- as.matrix(B)
  stopifnot(nrow(X) == nrow(B))
  n <- nrow(X); P <- ncol(X); R <- ncol(B)
  Z <- matrix(0.0, n, P*R)
  for (p in seq_len(P)) {
    cols <- ((p-1)*R + 1):(p*R)
    Z[, cols] <- B * X[, p]
  }
  Z
}
Z <- build_localized_Z(Xq, B)
stopifnot(nrow(Z) == M, nrow(Xq) == M)
## ========================= Likelihood + BIC helpers =========================
loglik_poisson_offset <- function(y, eta, weights = NULL) {
  mu <- exp(eta)
  if (is.null(weights)) sum(y * eta - mu) else sum(weights * (y * eta - mu))
}
bic_poisson_offset <- function(y, X, beta, offset, n_eff = length(y), weights = NULL) {
  eta <- as.numeric(X %*% beta + offset)
  ll  <- loglik_poisson_offset(y, eta, weights)
  s   <- sum(abs(beta) > 0)
  -2*ll + s * log(n_eff)
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

