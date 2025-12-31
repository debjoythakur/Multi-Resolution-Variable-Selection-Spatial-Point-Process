# ============================================================
# OVERVIEW (What this script does)
# ============================================================
# Goal:
#   Build a road-segment level modeling dataset for St. Louis where each road segment is a unit of analysis.
#   For each road, we (i) attach incident counts (accidents + crimes), and (ii) construct neighborhood-based
#   binary covariates (presence of buildings/POIs/rail/traffic/etc. within a buffer), plus optional road attributes.
#
# Main Steps:
#   1) Read layers from a GeoPackage (roads, buildings, POIs, railways, transport, waterways, traffic, accidents, crimes)
#      and transform everything to a metric CRS (UTM 15N) so distances and buffers are in meters.
#
#   2) Clean accident timestamps:
#      - Parse Accident_D as Date
#      - Convert Accident_T (e.g., 645 -> "06:45") to a HH:MM time
#      - Build a POSIXct datetime (America/Chicago)
#      - Subset accidents from 2025-01-01 onward
#
#   3) Compute per-road geometry properties:
#      - Road segment length (meters) from st_length()
#      - Drop invalid/zero-length segments
#
#   4) Snap point events to roads (spatial linking):
#      - For each accident/crime point, find nearest road segment (st_nearest_feature)
#      - Keep the match only if distance <= snap_thr_m (30m); otherwise set road_id = NA
#
#   5) Aggregate incidents to road-level:
#      - Per road: accident count, crime count
#      - Robustly detect/convert killed/injured fields (column names may vary across exports)
#      - Create totals (killed_total, injured_total) and has_incident indicator
#
#   6) Create road “centers” and neighborhood buffers:
#      - Use st_point_on_surface() to get a representative point per road
#      - Store lon/lat (WGS84) for mapping or downstream models
#      - Build buffers (buffer_m) around road centers for “neighborhood” covariates
#
#   7) Build neighborhood covariates (presence dummies):
#      - For each layer (buildings, POIs, rail, traffic, etc.), create binary indicators:
#          prefix_category = 1 if any feature of that category intersects the road buffer
#      - Also compute neighborhood crime counts within the buffer (crime_total_neigh)
#
#   8) Add road’s own attributes (if available):
#      - e.g., fclass, maxspeed, bridge, tunnel
#      - Convert bridge/tunnel flags to 0/1 and one-hot encode road fclass
#
#   9) Assemble model_df:
#      - Join lon/lat + core incident variables + road attributes + neighborhood dummies
#      - Fill remaining NA dummy columns with 0 (keep key numeric columns intact)
#
#  10) Visualization utilities:
#      - Plot each covariate layer within a zoomed bounding box (full legend, large fonts)
#      - Identify “incident-rich” roads (many nonzero covariates) and draw local road network context
#      - Overlay predictor “marks” near the road centroid to show which covariates are present
#
# Notes:
#   - sf_use_s2(FALSE) is used to avoid s2 geometry issues for planar distance/buffer workflows.
#   - buffer_m controls neighborhood features; buf_m controls visualization extent around a chosen road.
# ============================================================

rm(list = ls(all = TRUE))
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(data.table); library(units); library(ggplot2)
})

## Data Collection and Data Cleaning

sf_use_s2(FALSE)

library(dplyr)
library(lubridate)
library(stringr)
library(sf)


## --------- 0) Paths & CRS ----------
gpkg_path <- "/Users/debjoythakur/Downloads/LMPPLE/missouri_stl.gpkg"
crs_metric <- 26915        # NAD83 / UTM zone 15N
snap_thr_m <- set_units(100, m)  # distance threshold for mapping events to roads
buffer_m   <- 0

## --------- 1) Read layers ----------
st_layers(gpkg_path)

roads <- st_read(gpkg_path, layer = "clipped-road") %>%
  st_transform(crs_metric) %>%
  mutate(road_id = row_number())

buildings   <- st_read(gpkg_path, layer = "clipped-buildings") %>% st_transform(crs_metric)
pois        <- st_read(gpkg_path, layer = "clipped-pois") %>% st_transform(crs_metric)
railways    <- st_read(gpkg_path, layer = "clipped-rail") %>% st_transform(crs_metric)
transport   <- st_read(gpkg_path, layer = "clipped-transport") %>% st_transform(crs_metric)
waterways   <- st_read(gpkg_path, layer = "clipped-waterways") %>% st_transform(crs_metric)
traffic     <- st_read(gpkg_path, layer = "clipped-traffic") %>% st_transform(crs_metric)
traffic_ch  <- st_read(gpkg_path, layer = "clipped-traffic-ch") %>% st_transform(crs_metric)
accident    <- st_read(gpkg_path, layer = "clipped-accident") %>% st_transform(crs_metric)
crime_data  <- st_read(gpkg_path, layer = "clipped-crime") %>% st_transform(crs_metric)


# 1) Parse date and time robustly
accident <- accident %>%
  mutate(
    # Coerce Accident_D to Date (handles character/factor)
    Accident_D = as.Date(Accident_D),
    
    # Coerce time like 645 -> "06:45", 1045 -> "10:45"
    Accident_T_chr = sprintf("%04d", as.integer(Accident_T)),
    Accident_HH = substr(Accident_T_chr, 1, 2),
    Accident_MM = substr(Accident_T_chr, 3, 4),
    
    # Build a POSIXct datetime in local (St. Louis is America/Chicago)
    acc_datetime = ymd_hm(
      paste(Accident_D, paste0(Accident_HH, ":", Accident_MM)),
      tz = "America/Chicago"
    )
  )

# 2) Subset from 2025-01-01 onward (choose date or datetime filter)
start_dt <- ymd("2025-01-01", tz = "America/Chicago")

# If you care only about the date (ignore time within the day):
accident <- accident %>% filter(Accident_D >= as.Date(start_dt))

## --------- 2) Per-road length ----------
roads$len_m <- as.numeric(st_length(roads))
roads <- roads[is.finite(roads$len_m) & roads$len_m > 0, ]

## --------- 3) Map accidents & crimes to nearest road within 30 m ----------
## Accidents
idx_acc <- st_nearest_feature(accident, roads)
d_acc   <- st_distance(accident, roads[idx_acc, ], by_element = TRUE)

accident$road_id <- roads$road_id[idx_acc]
accident$road_id[d_acc > snap_thr_m] <- NA

## Crimes
idx_cr <- st_nearest_feature(crime_data, roads)
d_cr   <- st_distance(crime_data, roads[idx_cr, ], by_element = TRUE)

crime_data$road_id <- roads$road_id[idx_cr]
crime_data$road_id[d_cr > snap_thr_m] <- NA

## --------- 4) Build per-road counts & totals (robust to column names) ----------
to_num <- function(x) suppressWarnings(as.numeric(x))

## columns that *might* exist for killed/injured
kill_candidates_acc  <- c("Number_Kil","Number_Killed","Killed","Fatalities")
inj_candidates_acc   <- c("Number_Inj","Number_Injured","Injured")
kill_candidates_crim <- c("Number_Kil","Homicide_Victims","Victims_Killed","Killed","Fatalities")
inj_candidates_crim  <- c("Number_Inj","Victims_Injured","Injured")

## Replace the helper with this safer version
to_num_scalar <- function(el){
  # Convert a single cell to a numeric scalar (or NA)
  if (length(el) == 0) return(NA_real_)
  if (is.list(el)) el <- el[[1]]
  # If it's a character like "3" or "3.0", coerce; else NA
  val <- suppressWarnings(as.numeric(el))
  if (length(val) == 0) return(NA_real_)
  val
}

pick_col_or_zero <- function(df, candidates){
  cols <- intersect(candidates, names(df))
  if (length(cols) == 0) return(rep(0, nrow(df)))
  
  mats <- lapply(cols, function(cl){
    x <- df[[cl]]
    
    # Handle list-columns, character, factor, numeric robustly
    if (is.list(x)) {
      x <- vapply(x, to_num_scalar, numeric(1))
    } else {
      # factor/character/numeric -> numeric
      x <- suppressWarnings(as.numeric(x))
    }
    
    x[is.na(x)] <- 0
    x
  })
  
  rowSums(do.call(cbind, mats), na.rm = TRUE)
}


accident <- accident %>%
  mutate(
    killed_acc = pick_col_or_zero(., kill_candidates_acc),
    inj_acc    = pick_col_or_zero(., inj_candidates_acc)
  )

crime_data <- crime_data %>%
  mutate(
    killed_crime = pick_col_or_zero(., kill_candidates_crim),
    inj_crime    = pick_col_or_zero(., inj_candidates_crim)
  )

## Summaries by road
acc_by_road <- accident %>%
  st_drop_geometry() %>%
  filter(!is.na(road_id)) %>%
  group_by(road_id) %>%
  summarise(
    acc_count_road = n(),
    killed_acc_sum = sum(killed_acc, na.rm = TRUE),
    inj_acc_sum    = sum(inj_acc,    na.rm = TRUE),
    .groups = "drop"
  )

crime_by_road <- crime_data %>%
  st_drop_geometry() %>%
  filter(!is.na(road_id)) %>%
  group_by(road_id) %>%
  summarise(
    crime_count_road = n(),
    killed_crime_sum = sum(killed_crime, na.rm = TRUE),
    inj_crime_sum    = sum(inj_crime,    na.rm = TRUE),
    .groups = "drop"
  )

## Join to roads; create totals & has_incident
roads_events <- roads %>%
  left_join(acc_by_road,  by = "road_id") %>%
  left_join(crime_by_road, by = "road_id") %>%
  mutate(
    across(c(acc_count_road, crime_count_road,
             killed_acc_sum, inj_acc_sum,
             killed_crime_sum, inj_crime_sum),
           ~ replace_na(., 0L)),
    killed_total  = killed_acc_sum + killed_crime_sum,
    injured_total = inj_acc_sum    + inj_crime_sum,
    has_incident  = as.integer((acc_count_road + crime_count_road) > 0)
  )

## --------- 5) Road centers (for lon/lat) & 100 m buffers ----------
mid_sfc <- st_point_on_surface(st_geometry(roads_events))
roads_centers <- st_sf(st_drop_geometry(roads_events), geometry = mid_sfc, crs = st_crs(roads_events))

cent_ll <- st_transform(roads_centers, 4326)
coords  <- st_coordinates(cent_ll)
roads_centers$lon <- coords[,1]
roads_centers$lat <- coords[,2]

buf100 <- st_buffer(roads_centers, buffer_m)

## --------- 6) Helpers for neighborhood features ----------
sanitize_cat <- function(x){
  x <- tolower(trimws(as.character(x)))
  x[is.na(x) | x %in% c("", "na")] <- "na"
  x <- str_replace_all(x, "[^a-z0-9]+", "_")
  x <- str_replace_all(x, "_+", "_")
  str_remove_all(x, "^_|_$")
}

presence_block <- function(buffers, layer, field, prefix){
  if (!inherits(layer, "sf")) stop("Layer is not sf: ", deparse(substitute(layer)))
  layer <- st_transform(layer, st_crs(buffers)) %>%
    mutate(.cat = sanitize_cat(.data[[field]])) %>%
    select(.cat)
  hits <- st_join(buffers %>% select(road_id), layer %>% select(.cat),
                  join = st_intersects, left = FALSE) %>%
    st_drop_geometry() %>%
    distinct(road_id, .cat)
  if (nrow(hits) == 0L) return(buffers %>% st_drop_geometry() %>% select(road_id))
  hits %>%
    mutate(val = 1L) %>%
    pivot_wider(names_from = .cat, values_from = val, values_fill = 0L,
                names_prefix = paste0(prefix, "_")) %>%
    group_by(road_id) %>%
    summarise(across(everything(), ~ as.integer(any(. == 1L))), .groups = "drop")
}

## Optional: total count of crimes within 100 m (neighborhood; different from “on the road”)
crime_total_block <- function(buffers, crime_sf){
  crime_sf <- st_transform(crime_sf, st_crs(buffers))
  idx <- st_intersects(buffers, crime_sf, sparse = TRUE)
  tibble::tibble(
    road_id = buffers$road_id,
    crime_total_neigh = lengths(idx)
  )
}

## --------- 7) Neighborhood dummies (presence within 100 m) ----------
bldg_dummy   <- presence_block(buf100, buildings,  "type",    "bldg")
pois_dummy   <- presence_block(buf100, pois,       "fclass",  "pois")
rail_dummy   <- presence_block(buf100, railways,   "fclass",  "rail")
trans_dummy  <- presence_block(buf100, transport,  "fclass",  "trans")
water_dummy  <- presence_block(buf100, waterways,  "fclass",  "waterway")
trafp_dummy  <- presence_block(buf100, traffic,    "fclass",  "traffic")
trafch_dummy <- presence_block(buf100, traffic_ch, "fclass",  "trafficch")
crime_neigh  <- crime_total_block(buf100, crime_data)

## --------- 8) Road own attributes (if present) ----------
road_attr_cols <- intersect(c("fclass","maxspeed","bridge","tunnel"), names(roads_events))
road_attr <- roads_events %>%
  st_drop_geometry() %>%
  select(any_of(c("road_id", road_attr_cols)))

road_fclass_dmy <- NULL
if ("fclass" %in% names(road_attr)) {
  road_fclass_dmy <- road_attr %>%
    mutate(fclass = sanitize_cat(fclass)) %>%
    transmute(road_id, col = paste0("road_fclass_", fclass), val = 1L) %>%
    distinct() %>%
    pivot_wider(names_from = col, values_from = val, values_fill = 0L)
}

road_attr_final <- road_attr %>%
  { if ("bridge" %in% names(.)) mutate(., bridge = if_else(bridge %in% c("T","true","True","1","Y","y"), 1L, 0L)) else . } %>%
  { if ("tunnel" %in% names(.)) mutate(., tunnel = if_else(tunnel %in% c("T","true","True","1","Y","y"), 1L, 0L)) else . } %>%
  left_join(road_fclass_dmy %||% tibble(road_id = .$road_id), by = "road_id")

## --------- 9) Assemble FINAL model_df (no risk_per_100m) ----------
centers_ll <- roads_centers %>% st_drop_geometry() %>%
  select(road_id, lon, lat)

core_covars <- roads_events %>%
  st_drop_geometry() %>%
  select(road_id, len_m, has_incident,
         crime_count_road, acc_count_road,
         killed_total, injured_total)

model_df <- list(
  centers_ll,
  core_covars,
  road_attr_final,
  bldg_dummy, pois_dummy, rail_dummy, trans_dummy, water_dummy, trafp_dummy, trafch_dummy,
  crime_neigh
) %>% reduce(~ left_join(.x, .y, by = "road_id"))

## Fill remaining NAs in dummies with 0 (keep lon/lat/len_m/maxspeed intact)
zero_cols <- setdiff(names(model_df),
                     c("road_id","lon","lat","len_m","has_incident",
                       "crime_count_road","acc_count_road","killed_total","injured_total","maxspeed"))
for (cl in zero_cols) {
  if (is.numeric(model_df[[cl]]) || is.integer(model_df[[cl]])) {
    model_df[[cl]][is.na(model_df[[cl]])] <- 0
  }
}

## --------- 10) Inspect ----------
print(dim(model_df))
print(head(model_df[, c("road_id","lon","lat","len_m",
                        "has_incident","crime_count_road","acc_count_road",
                        "killed_total","injured_total")]))

# --- 1) Diagnose
print(inherits(roads_events, "sf"))        # should be TRUE
print(st_crs(roads_events))                 # should be non-NULL
print(is.null(st_geometry(roads_events)))   # should be FALSE

# --- 2) Repair if geometry got dropped along the way
if (!inherits(roads_events, "sf") || is.null(st_geometry(roads_events))) {
  # Rebuild by joining the attributes you computed back onto the original roads geometry
  attrs <- roads_events %>% 
    dplyr::select(-any_of(attr(roads, "sf_column"))) %>% 
    sf::st_drop_geometry()                  # attributes only (has_incident, counts, etc.)
  
  roads_events <- roads %>%                 # roads still has geometry
    dplyr::select(road_id, geometry = !!attr(roads, "sf_column")) %>%
    dplyr::left_join(attrs, by = "road_id") %>%
    sf::st_as_sf(crs = sf::st_crs(roads))   # ensure sf tagged with same CRS
}


# ============================================================
# ONE-SNIPPET (UPDATED + bigger fonts): Zoomed base-R plots for ALL covariates
# - Focus bbox: lon/lat bbox (edit below)
# - Each layer plotted SEPARATELY (one plot per covariate layer)
# - fclass/type auto-detected per layer and shown via different pch
# - NO "top k" restriction: legend shows ALL types/classes
# - Legends placed OUTSIDE on the RIGHT (no clipping if device wide enough)
# - Slightly larger font sizes for titles + legend text
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

sf_use_s2(FALSE)

# ----------------------------
# 0) BBOX (lon/lat) and CRS handling
# ----------------------------
bb_ll <- st_as_sfc(st_bbox(c(xmin=-90.30, xmax=-90.28, ymin=38.78, ymax=38.80), crs=4326))

ref_crs <- st_crs(roads)
if (is.na(ref_crs)) ref_crs <- st_crs(roads_events)
bb_ref <- st_transform(bb_ll, ref_crs)

# Use st_crop for a bbox (simpler + usually faster than st_intersects)
clip_sf <- function(x, bb_poly) {
  if (!is.na(st_crs(x)) && !is.na(st_crs(bb_poly)) && st_crs(x) != st_crs(bb_poly)) {
    x <- st_transform(x, st_crs(bb_poly))
  }
  x <- st_make_valid(x)
  st_crop(x, st_bbox(bb_poly))
}

# ----------------------------
# 1) Auto-detect the "class" column for each layer
# ----------------------------
pick_class_col <- function(x, candidates) {
  candidates <- intersect(candidates, names(x))
  if (length(candidates) == 0) return(NA_character_)
  scores <- sapply(candidates, function(cc) {
    v <- as.character(x[[cc]])
    v <- v[!is.na(v) & v != "" & v != "<NA>"]
    length(unique(v))
  })
  candidates[which.max(scores)]
}

# ----------------------------
# 2) Plotting helper (FULL legend, bigger fonts, legend outside)
# ----------------------------
plot_layer_bw_all <- function(x, base_roads, bb_poly, main,
                              class_candidates = c("fclass","type","building","amenity","class","subclass","ref","name"),
                              max_points = 2500,
                              legend_ncol = 2,
                              title_cex  = 1.6,   # bigger title
                              legend_cex = 0.95,  # bigger legend text
                              pt_cex     = 0.85,  # bigger symbols
                              right_mar  = 18) {  # extra right margin for legend
  
  bb <- st_bbox(bb_poly)
  
  # margins: big right margin to fit legend text
  op <- par(no.readonly = TRUE)
  par(mar = c(2, 2, 4, right_mar), xpd = NA, cex.main = title_cex)
  
  on.exit(par(op), add = TRUE)
  
  # empty plot handler
  if (nrow(x) == 0) {
    plot(NA, xlim=c(bb["xmin"], bb["xmax"]), ylim=c(bb["ymin"], bb["ymax"]),
         asp=1, axes=FALSE, xlab="", ylab="", main=main)
    box()
    if (!is.null(base_roads) && nrow(base_roads) > 0) {
      plot(st_geometry(base_roads), add=TRUE, col="grey85", lwd=0.7)
    }
    legend("topright", inset = c(-0.30, 0), bty="n", cex=legend_cex,
           title="(none)", legend="(no features in bbox)")
    return(invisible(NULL))
  }
  
  # match CRS for base roads
  if (!is.null(base_roads) && nrow(base_roads) > 0 &&
      !is.na(st_crs(base_roads)) && !is.na(st_crs(x)) && st_crs(base_roads) != st_crs(x)) {
    base_roads <- st_transform(base_roads, st_crs(x))
  }
  
  # detect class column
  class_col <- pick_class_col(x, class_candidates)
  cls <- if (!is.na(class_col) && class_col %in% names(x)) as.character(x[[class_col]]) else rep("unknown", nrow(x))
  cls[is.na(cls) | cls=="" | cls=="<NA>"] <- "unknown"
  
  # FULL set of classes (no restriction)
  ucls <- sort(unique(cls))
  
  # pch mapping (repeats if too many classes)
  pch_base <- 0:25
  pch_map  <- setNames(pch_base[(seq_along(ucls)-1) %% length(pch_base) + 1], ucls)
  if (length(ucls) > length(pch_base)) {
    message("NOTE: ", length(ucls), " classes but only ", length(pch_base),
            " base-R pch symbols. Symbols WILL REPEAT in legend.")
  }
  
  # geometry mode
  gt <- unique(as.character(st_geometry_type(x)))
  mode <- if (any(grepl("POINT", gt))) "point" else if (any(grepl("LINE", gt))) "line" else "poly"
  
  # plot canvas
  plot(NA, xlim=c(bb["xmin"], bb["xmax"]), ylim=c(bb["ymin"], bb["ymax"]),
       asp=1, axes=FALSE, xlab="", ylab="", main=main)
  box()
  
  # base roads
  if (!is.null(base_roads) && nrow(base_roads) > 0) {
    plot(st_geometry(base_roads), add=TRUE, col="grey88", lwd=0.7)
  }
  
  # draw geometry lightly
  if (mode == "line") {
    plot(st_geometry(x), add=TRUE, col="grey60", lwd=0.9)
  } else if (mode == "poly") {
    plot(st_geometry(x), add=TRUE, border="grey60", col=NA, lwd=0.7)
  }
  
  # coordinates for marks
  if (mode == "point") {
    xy <- st_coordinates(st_geometry(x))
  } else {
    xy <- st_coordinates(st_centroid(st_geometry(x)))
  }
  
  # downsample for readability (optional)
  if (nrow(xy) > max_points) {
    set.seed(1)
    keep <- sample(seq_len(nrow(xy)), max_points)
    xy  <- xy[keep, , drop=FALSE]
    cls <- cls[keep]
  }
  
  # marks
  for (k in seq_along(ucls)) {
    sel <- cls == ucls[k]
    if (any(sel)) points(xy[sel,1], xy[sel,2], pch=pch_map[[ucls[k]]], cex=pt_cex)
  }
  
  # legend OUTSIDE to the right
  legend("topright",
         inset = c(-0.32, 0),
         bty = "n",
         cex = legend_cex,
         title = {
           # Capitalize the class name in legend title (e.g., "type" -> "Type", "fclass" -> "Fclass")
           tt <- ifelse(is.na(class_col), "Class", class_col)
           paste0(toupper(substr(tt,1,1)), substring(tt,2))
         },
         legend = ucls,
         pch = unname(pch_map),
         ncol = legend_ncol)
}

# ----------------------------
# 3) Clip ALL layers to bbox
# ----------------------------
roads_z      <- clip_sf(roads,      bb_ref)
buildings_z  <- clip_sf(buildings,  bb_ref)
pois_z       <- clip_sf(pois,       bb_ref)
railways_z   <- clip_sf(railways,   bb_ref)
traffic_z    <- clip_sf(traffic,    bb_ref)
traffic_ch_z <- clip_sf(traffic_ch, bb_ref)
transport_z  <- clip_sf(transport,  bb_ref)
waterways_z  <- clip_sf(waterways,  bb_ref)

# ----------------------------
# 4) Make SEPARATE plots (each call creates its own figure)
# ----------------------------
plot_layer_bw_all(buildings_z,  base_roads=roads_z, bb_poly=bb_ref, main="Buildings",
                  class_candidates=c("type","building","fclass","class","subclass"),
                  legend_ncol=2, title_cex=1.7, legend_cex=0.95, pt_cex=0.90) 

plot_layer_bw_all(pois_z,       base_roads=roads_z, bb_poly=bb_ref, main="POIs",
                  class_candidates=c("fclass","amenity","type","class","subclass","name"),
                  legend_ncol=1, title_cex=1.7, legend_cex=0.95, pt_cex=0.90)

plot_layer_bw_all(railways_z,   base_roads=roads_z, bb_poly=bb_ref, main="Railways",
                  class_candidates=c("fclass","type","class","subclass"),
                  legend_ncol=2, title_cex=1.7, legend_cex=0.95, pt_cex=0.90)

plot_layer_bw_all(waterways_z,  base_roads=roads_z, bb_poly=bb_ref, main="Waterways",
                  class_candidates=c("fclass","type","class","subclass","name"),
                  legend_ncol=2, title_cex=1.7, legend_cex=0.95, pt_cex=0.90)

plot_layer_bw_all(traffic_z,    base_roads=roads_z, bb_poly=bb_ref, main="Traffic",
                  class_candidates=c("fclass","type","class","subclass","name"),
                  legend_ncol=1, title_cex=1.7, legend_cex=0.95, pt_cex=0.90)

plot_layer_bw_all(traffic_ch_z, base_roads=roads_z, bb_poly=bb_ref, main="Traffic_ch",
                  class_candidates=c("fclass","type","class","subclass","name"),
                  legend_ncol=1, title_cex=1.7, legend_cex=0.95, pt_cex=0.90)

plot_layer_bw_all(transport_z,  base_roads=roads_z, bb_poly=bb_ref, main="Transport",
                  class_candidates=c("fclass","type","class","subclass","name"),
                  legend_ncol=1, title_cex=1.7, legend_cex=0.95, pt_cex=0.90)


# --- 3) Make sure has_incident is there and well-typed
if (!"has_incident" %in% names(roads_events)) {
  stop("Column 'has_incident' is missing. Re-run the step that creates it.")
}
roads_events$has_incident <- as.integer(roads_events$has_incident > 0)


sf_use_s2(FALSE)

# -----------------------------
# Settings
# -----------------------------
buf_m     <- 100     # radius used ONLY for extent / local network selection (no circle drawn)
start_nz  <- 10      # start threshold; will auto-relax if needed
max_marks <- 12      # max predictor symbols drawn per road

exclude_cols <- c(
  "road_id", "lon", "lat", "len_m", "maxspeed",
  "has_incident",
  "crime_count_road", "acc_count_road",
  "killed_total", "injured_total"
)

# Use your important_covariates if available, otherwise count all non-excluded cols
cov_pool <- if (exists("important_covariates")) {
  intersect(important_covariates, names(model_df))
} else {
  setdiff(names(model_df), exclude_cols)
}
count_cols <- setdiff(cov_pool, exclude_cols)

# -----------------------------
# 0) Helper: pick a road-name column (robust to different OSM exports)
# -----------------------------
pick_name_col <- function(x) {
  cand <- c("name", "road_name", "street_name", "fullname", "ref", "label", "osm_name")
  cand <- intersect(cand, names(x))
  if (length(cand) == 0) return(NA_character_)
  # prefer columns with more non-missing values
  nn <- sapply(cand, function(cc) sum(!is.na(x[[cc]]) & as.character(x[[cc]]) != ""))
  cand[which.max(nn)]
}

name_col <- pick_name_col(roads_events)

get_road_label <- function(road_id) {
  rs <- roads_events %>% filter(road_id == !!road_id)
  if (nrow(rs) == 0) return(as.character(road_id))
  if (!is.na(name_col) && name_col %in% names(rs)) {
    nm <- as.character(rs[[name_col]][1])
    nm <- trimws(nm)
    if (!is.na(nm) && nm != "") return(nm)
  }
  # fallback
  paste0("road_id=", road_id)
}

# -----------------------------
# 1) Compute nonzero_count for all roads
# -----------------------------
model_df2 <- model_df %>%
  mutate(
    nonzero_count = rowSums(across(all_of(count_cols), ~ !is.na(.) & (. != 0)))
  )

# -----------------------------
# 2) Find 5 incident roads, relaxing threshold if necessary
# -----------------------------
pick_5_incident_rich <- function(df, start_thr = 10) {
  for (thr in seq(start_thr, 1, by = -1)) {
    cand <- df %>%
      filter(has_incident == 1, nonzero_count >= thr) %>%
      arrange(desc(nonzero_count), desc(acc_count_road), desc(crime_count_road))
    if (nrow(cand) >= 5) {
      message("Using threshold nonzero_count >= ", thr, " (found ", nrow(cand), " roads).")
      return(list(threshold = thr, roads5 = cand %>% slice_head(n = 5)))
    }
  }
  stop("No incident roads found at any threshold. Check has_incident coding in model_df.")
}

pick <- pick_5_incident_rich(model_df2, start_thr = start_nz)
thr_used <- pick$threshold
roads5_df <- pick$roads5

# -----------------------------
# 3) Plot helper (one road)
# -----------------------------
plot_one_road_with_marks <- function(road_id) {
  
  road_sf <- roads_events %>% filter(road_id == !!road_id) %>% st_make_valid()
  if (nrow(road_sf) != 1) stop("road_id not found (or duplicated) in roads_events: ", road_id)
  
  # buffer used only for extent + local-road selection (circle NOT drawn)
  buf <- st_union(st_buffer(road_sf, buf_m))
  bb  <- st_bbox(buf)
  
  # local road network within buffer
  roads_local <- roads_events %>%
    st_make_valid() %>%
    st_transform(st_crs(road_sf)) %>%
    { .[ as.vector(st_intersects(., buf, sparse = FALSE)), ] }
  
  df1 <- roads5_df %>% filter(road_id == !!road_id)
  
  present_covs <- count_cols[ which(df1[1, count_cols, drop=TRUE] != 0) ]
  present_covs <- head(present_covs, max_marks)
  
  road_label <- get_road_label(road_id)
  
  main_txt <- paste0(
    road_label,
    " | nz=", df1$nonzero_count,
    " | acc=", df1$acc_count_road,
    " | crime=", df1$crime_count_road,
    " | thr>=", thr_used
  )
  
  plot(NA,
       xlim = c(bb["xmin"], bb["xmax"]),
       ylim = c(bb["ymin"], bb["ymax"]),
       asp = 1, axes = FALSE, xlab = "", ylab = "", main = main_txt)
  box()
  
  # (REMOVED) buffer circle drawing:
  # plot(st_geometry(buf), add = TRUE, border = "grey70", lwd = 1.0, lty = 2)
  
  if (nrow(roads_local) > 0) plot(st_geometry(roads_local), add = TRUE, col = "grey80", lwd = 0.7)
  plot(st_geometry(road_sf), add = TRUE, col = "black", lwd = 3.0)
  
  # predictor marks near centroid (jittered)
  cen <- st_centroid(st_geometry(road_sf))
  xy  <- st_coordinates(cen)
  
  set.seed(1)
  J <- 22
  jit <- function() c(xy[1] + runif(1, -J, J), xy[2] + runif(1, -J, J))
  
  pch_pool <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
  pch_use  <- pch_pool[seq_len(min(length(present_covs), length(pch_pool)))]
  
  for (k in seq_along(present_covs)) {
    p <- jit()
    points(p[1], p[2], pch = pch_use[k], cex = 1.15)
  }
  
  legend("bottomleft", bty = "n",
         legend = c("Local road network", "Selected road"),
         col = c("grey80", "black"),
         lwd = c(0.7, 3.0), lty = c(1,1), cex = 0.85)
  
  if (length(present_covs) > 0) {
    legend("right", inset = c(-0.28, 0), xpd = TRUE, bty = "n",
           title = paste0("Predictors shown (", length(present_covs), ")"),
           legend = present_covs, pch = pch_use, cex = 0.65)
  }
}

# -----------------------------
# 4) Make 5 plots (3x2 layout)
# -----------------------------
op <- par(mfrow = c(3, 2), mar = c(1,1,3,10))

for (rid in roads5_df$road_id) {
  plot_one_road_with_marks(rid)
}
plot.new()  # empty 6th slot

par(op)


# ---- you already have: roads_events, model_df, count_cols, buf_m ----
# If you don't have model_df2:
model_df2 <- model_df %>%
  mutate(nonzero_count = rowSums(across(all_of(count_cols), ~ !is.na(.) & (. != 0))))

# --- road-name helpers (same as before) ---
pick_name_col <- function(x) {
  cand <- c("name","road_name","street_name","fullname","ref","label","osm_name")
  cand <- intersect(cand, names(x))
  if (length(cand) == 0) return(NA_character_)
  nn <- sapply(cand, function(cc) sum(!is.na(x[[cc]]) & as.character(x[[cc]]) != ""))
  cand[which.max(nn)]
}
name_col <- pick_name_col(roads_events)

get_road_label <- function(road_id) {
  rs <- roads_events %>% filter(road_id == !!road_id)
  if (nrow(rs) == 0) return(as.character(road_id))
  if (!is.na(name_col) && name_col %in% names(rs)) {
    nm <- trimws(as.character(rs[[name_col]][1]))
    if (!is.na(nm) && nm != "") return(nm)
  }
  paste0("road_id=", road_id)
}

# --- SAFE "present covariate" detector (prevents NAs introduced by coercion) ---
is_present_val <- function(v) {
  if (is.list(v)) v <- v[[1]]
  v <- v[1]
  
  if (is.logical(v)) return(!is.na(v) && v)
  if (is.numeric(v) || is.integer(v)) return(!is.na(v) && v != 0)
  
  # character / factor
  vc <- as.character(v)
  if (is.na(vc) || vc == "") return(FALSE)
  
  suppressWarnings(num <- as.numeric(vc))
  if (!is.na(num)) return(num != 0)
  
  vc %in% c("T","TRUE","true","Y","YES","Yes","yes")
}

# -----------------------------
# Plot helper (one road)
# -----------------------------
max_marks <- 12  # if you want, increase
plot_one_road_with_marks <- function(road_id) {
  
  road_sf <- roads_events %>% filter(road_id == !!road_id) %>% st_make_valid()
  if (nrow(road_sf) != 1) stop("road_id not found (or duplicated) in roads_events: ", road_id)
  
  buf <- st_union(st_buffer(road_sf, buf_m))     # used ONLY for extent/local network
  bb  <- st_bbox(buf)
  
  roads_local <- roads_events %>%
    st_make_valid() %>%
    st_transform(st_crs(road_sf)) %>%
    { .[ as.vector(st_intersects(., buf, sparse = FALSE)), ] }
  
  df1 <- model_df2 %>% filter(road_id == !!road_id)
  if (nrow(df1) != 1) stop("road_id not found (or duplicated) in model_df2: ", road_id)
  
  vals <- df1[1, count_cols, drop = FALSE]
  present_covs <- count_cols[ sapply(vals, is_present_val) ]
  present_covs <- head(present_covs, max_marks)
  
  main_txt <- paste0(
    get_road_label(road_id),
    " | accident count=", df1$acc_count_road,
    " | crime count =", df1$crime_count_road
  )
  
  plot(NA,
       xlim = c(bb["xmin"], bb["xmax"]),
       ylim = c(bb["ymin"], bb["ymax"]),
       asp = 1, axes = FALSE, xlab = "", ylab = "", main = main_txt)
  box()
  
  # NO buffer circle drawn
  
  if (nrow(roads_local) > 0) plot(st_geometry(roads_local), add = TRUE, col = "grey80", lwd = 0.7)
  plot(st_geometry(road_sf), add = TRUE, col = "black", lwd = 3.0)
  
  # predictor marks near centroid (jittered)
  cen <- st_centroid(st_geometry(road_sf))
  xy  <- st_coordinates(cen)
  
  set.seed(1)
  J <- 22
  jit <- function() c(xy[1] + runif(1, -J, J), xy[2] + runif(1, -J, J))
  
  pch_pool <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
  pch_use  <- pch_pool[seq_len(min(length(present_covs), length(pch_pool)))]
  
  for (k in seq_along(present_covs)) {
    p <- jit()
    points(p[1], p[2], pch = pch_use[k], cex = 1.15)
  }
  
  legend("bottomleft", bty = "n",
         legend = c("Local road network", "Selected road"),
         col = c("grey80", "black"),
         lwd = c(0.7, 3.0), lty = c(1,1), cex = 0.9)
  
  if (length(present_covs) > 0) {
    legend("topright", inset = c(-0.28, 0), xpd = TRUE, bty = "n",
           title = paste0("Predictors shown"),
           legend = present_covs, pch = pch_use, cex = 0.75)
  }
}

# -----------------------------
# ONLY keep these two roads in a 2x1 layout
# -----------------------------
keep_ids <- c(58100, 279705)

op <- par(mfrow = c(2, 1), mar = c(1,1,3,18))
plot_one_road_with_marks(keep_ids[1])
plot_one_road_with_marks(keep_ids[2])
par(op)