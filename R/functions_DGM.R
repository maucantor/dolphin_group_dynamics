# dolphin_group_metrics.R
# title: "Dolphin-group-dynamics-metrics"
# author: "Mauricio Cantor and João Valle-Pereira"
# date: "2025-11-26, 27", 2025-12-02"
#
# Modular pipeline to compute group-level metrics from drone detections.
#
# Required columns in the CSV output from dolphin-tracker:
#  - FrameIndex   : integer frame number
#  - ObjectID     : unique ID per tracked individual (may be numeric or character)
#  - CenterX_m    : x coordinate in meters
#  - CenterY_m    : y coordinate in meters
#  - MovingAvgAngle_deg : heading angle in degrees (optional but recommended)
# Optional:
#  - BreathCol : name of a binary column marking surfacing/breathing (1 / 0). Presence detection will be used as a proxy.
#
# Defaults:
#  fps = 30, window_sec = 5, Nmin = 2, interp_gap = 1 (frames)
#


# ---------------------------
#  Load packages and Helper functions
# ---------------------------

library(tidyverse)
library(zoo)       
library(igraph)    
library(viridis)   
library(patchwork) 




# safe_chull_area:
# Compute convex hull area for a set of 2D points.
# Returns 0 for two points, NA if no points.
safe_chull_area <- function(x, y) {
  if (length(x) == 0) return(NA_real_)
  if (length(x) == 1) return(0)
  if (length(x) == 2) {
    # area of degenerate polygon = 0
    return(0)
  }
  idx <- chull(x, y)
  # close polygon
  xs <- c(x[idx], x[idx[1]])
  ys <- c(y[idx], y[idx[1]])
  area <- abs(sum(xs[-1] * ys[-length(ys)] - xs[-length(xs)] * ys[-1])) / 2
  return(area)
}

# interpolate_tracks:
# Fill short detection gaps per ObjectID using linear interpolation up to maxgap frames.
# - df: tibble with FrameIndex, ObjectID, CenterX_m, CenterY_m, optional angle col
# - frame_col, id_col, x_col, y_col, angle_col (string names)
# - maxgap: maximum number of consecutive missing frames to interpolate
interpolate_tracks <- function(df,
                               frame_col = "FrameIndex",
                               id_col = "ObjectID",
                               x_col = "CenterX_m",
                               y_col = "CenterY_m",
                               angle_col = "MovingAvgAngle_deg",
                               maxgap = 1) {
  # ensure frame col numeric and sorted
  df <- df %>% arrange(.data[[id_col]], .data[[frame_col]])
  all_ids <- unique(df[[id_col]])
  out_list <- vector("list", length(all_ids))
  i <- 1
  for (id in all_ids) {
    sub <- df %>% filter(.data[[id_col]] == id)
    frmin <- min(sub[[frame_col]], na.rm = TRUE)
    frmax <- max(sub[[frame_col]], na.rm = TRUE)
    full_idx <- tibble(!!frame_col := seq(frmin, frmax))
    sub2 <- full_idx %>%
      left_join(sub, by = frame_col)
    # ensure id column exists
    sub2[[id_col]] <- id
    # interpolate x,y,angle with maxgap
    sub2[[x_col]] <- na.approx(sub2[[x_col]], maxgap = maxgap, na.rm = FALSE)
    sub2[[y_col]] <- na.approx(sub2[[y_col]], maxgap = maxgap, na.rm = FALSE)
    if (angle_col %in% names(sub2)) {
      sub2[[angle_col]] <- na.approx(sub2[[angle_col]], maxgap = maxgap, na.rm = FALSE)
    }
    out_list[[i]] <- sub2
    i <- i + 1
  }
  out <- bind_rows(out_list)
  # place columns consistently
  return(out)
}

# frames_presence_matrix:
# Build binary presence matrix (rows = ObjectID, columns = FrameIndex) using x_col as presence indicator
frames_presence_matrix <- function(df_interp, id_col = "ObjectID", frame_col = "FrameIndex", x_col = "CenterX_m") {
  pres <- df_interp %>%
    mutate(pres = ifelse(!is.na(.data[[x_col]]), 1, 0)) %>%
    select(all_of(c(id_col, frame_col, "pres"))) %>%
    pivot_wider(names_from = all_of(frame_col), values_from = pres, values_fill = 0)
  rownames_pres <- pres[[id_col]]
  M <- pres %>% select(-all_of(id_col)) %>% as.matrix()
  rownames(M) <- rownames_pres
  return(M)
}


# ---------------------------
# Per-frame metric functions
# ---------------------------

# 1. MPD_f: Mean Pairwise Distance per frame.
# Formula: MPD = mean_{i<j} d_{ij}
# Rationale: direct measure of group spatial compactness. Lower = tighter group.
# Biological meaning: smaller MPD indicates cohesive group behavior.
compute_mpd_frame <- function(points_df) {
  # points_df: tibble with CenterX_m, CenterY_m for that frame
  pts <- points_df %>% select(CenterX_m, CenterY_m) %>% drop_na()
  n <- nrow(pts)
  if (n < 2) return(NA_real_)
  dvec <- as.vector(dist(as.matrix(pts)))
  mpd <- mean(dvec)
  return(mpd)
}

# 2. MNN_f: Median Nearest Neighbor distance per frame (robust).
# Formula: MNN = median_i ( min_{j != i} d_{ij} )
# Rationale: robust local spacing metric, less sensitive to far-out outliers.
# Biological meaning: typical spacing to nearest neighbor; informs local social spacing.
compute_mnn_frame <- function(points_df) {
  pts <- points_df %>% select(CenterX_m, CenterY_m) %>% drop_na()
  n <- nrow(pts)
  if (n < 2) return(NA_real_)
  dmat <- as.matrix(dist(as.matrix(pts)))
  diag(dmat) <- NA
  nn <- apply(dmat, 1, min, na.rm = TRUE)
  return(median(nn, na.rm = TRUE))
}

# 3. CHAI_f: Convex Hull Area per individual (area normalized by N)
# Formula: CHAI = Area(convex hull of points) / N
# Rationale: normalizes group footprint by size so area is comparable across N.
# Biological meaning: captures overall spread per animal; lower CHAI implies tighter per-capita spacing.
compute_chai_frame <- function(points_df) {
  pts <- points_df %>% select(CenterX_m, CenterY_m) %>% drop_na()
  n <- nrow(pts)
  if (n == 0) return(NA_real_)
  area <- safe_chull_area(pts$CenterX_m, pts$CenterY_m)
  return(area / max(1, n))
}

# 4. Group Centroid Distance (MDC) (not sure this is useful)
# Formula: MDC = mean_i || p_i - centroid ||
# Rationale: measures central concentration and anisotropy.
# Biological meaning: reveals presence of a core and peripheral individuals.
compute_centroid_distance_frame <- function(points_df) {
  pts <- points_df %>% select(CenterX_m, CenterY_m) %>% drop_na()
  n <- nrow(pts)
  if (n == 0) return(NA_real_)
  cx <- mean(pts$CenterX_m)
  cy <- mean(pts$CenterY_m)
  dists <- sqrt((pts$CenterX_m - cx)^2 + (pts$CenterY_m - cy)^2)
  return(mean(dists))
}

# 5. Standard Distance (SD)
# Formula: SD = sqrt[(Σ(xi - x̄)²)/(n-1) + (Σ(yi - ȳ)²)/(n-1)]
# Rationale: measures group spread from centroid, combining x and y deviations.
# Biological meaning: lower SD = tighter clustering around center; higher SD = more dispersed.
# Reference: Proposed metrics from research framework PDF
compute_standard_distance_frame <- function(points_df) {
  pts <- points_df %>% select(CenterX_m, CenterY_m) %>% drop_na()
  n <- nrow(pts)
  if (n < 2) return(NA_real_)
  cx <- mean(pts$CenterX_m)
  cy <- mean(pts$CenterY_m)
  sd_x <- sum((pts$CenterX_m - cx)^2) / (n - 1)
  sd_y <- sum((pts$CenterY_m - cy)^2) / (n - 1)
  return(sqrt(sd_x + sd_y))
}

# 6. Dispersion Index
# Formula: Dispersion = MaxDistance_m / Mean_Distance_m
# Rationale: ratio of maximum to mean distance, indicates evenness of spacing.
# Biological meaning: higher values = more dispersed/uneven groups.
# Reference: Proposed metrics from research framework PDF
compute_dispersion_index_frame <- function(points_df) {
  pts <- points_df %>% select(CenterX_m, CenterY_m) %>% drop_na()
  n <- nrow(pts)
  if (n < 2) return(NA_real_)
  dmat <- as.matrix(dist(as.matrix(pts)))
  max_dist <- max(dmat)
  mean_dist <- mean(dmat[lower.tri(dmat)])
  if (mean_dist == 0) return(NA_real_)
  return(max_dist / mean_dist)
}

# 7. Minimum Spanning Tree (MST) mean edge length
# Formula: construct MST on pairwise distances; MST_mean = mean(edge lengths)
# Rationale: captures topological connectedness, robust to outliers.
# Biological meaning: low MST mean suggests close functional connectivity among individuals.
compute_mst_mean_frame <- function(points_df) {
  pts <- points_df %>% select(CenterX_m, CenterY_m) %>% drop_na()
  n <- nrow(pts)
  if (n < 2) return(NA_real_)
  dmat <- as.matrix(dist(as.matrix(pts)))
  g <- graph_from_adjacency_matrix(dmat, mode = "undirected", weighted = TRUE, diag = FALSE)
  mst_g <- mst(g, weights = E(g)$weight)
  weights <- E(mst_g)$weight
  if (length(weights) == 0) return(0)
  return(mean(weights))
}

# 8. Heading Similarity / Directional Coordination
# We compute two related measures:
#  - Heading_MRL: mean resultant length (MRL) of unit heading vectors (range 0..1)
#    Formula: R = || sum_i (u_i) || / n, where u_i = unit vector of heading angle
#  - Heading_pairwise_similarity: 1 - (mean absolute angular difference / pi), scaled 0..1
# Rationale: captures alignment across group in heading direction.
# Biological meaning: high values indicate traveling/aligned states; low values indicate milling/divergence.
compute_heading_similarity_frame <- function(points_df, angle_col = "MovingAvgAngle_deg") {
  if (!(angle_col %in% names(points_df))) return(list(MRL = NA_real_, PairwiseSim = NA_real_))
  angs <- points_df[[angle_col]] %>% na.omit() %>% as.numeric()
  n <- length(angs)
  if (n == 0) return(list(MRL = NA_real_, PairwiseSim = NA_real_))
  # convert to radians
  rad <- angs * pi / 180
  # Mean resultant length
  R <- sqrt((sum(cos(rad)))^2 + (sum(sin(rad)))^2) / n
  # pairwise absolute angular differences
  if (n >= 2) {
    radmat <- abs(outer(rad, rad, "-"))
    # wrap angles > pi
    radmat <- pmin(radmat, 2*pi - radmat)
    iu <- which(upper.tri(radmat), arr.ind = TRUE)
    diffs <- radmat[iu]
    pair_mean_diff <- mean(diffs, na.rm = TRUE)
    pairwise_sim <- 1 - (pair_mean_diff / pi) # scale to 0..1
  } else {
    pairwise_sim <- NA_real_
  }
  return(list(MRL = R, PairwiseSim = pairwise_sim))
}

# 9. Breathing / Surfacing Synchrony
# Input: either a dedicated Breath column (binary 0/1), or we use presence (detected => 1) as proxy.
# Per-frame synchrony S(t) = fraction of dyads that share the same binary state (0 or 1).
# Range 0..1, where 1 = all same state (all up or all down).
compute_breathing_synchrony_frame <- function(points_df, breath_col = NULL) {
  # points_df should contain ObjectID and either breath_col or presence proxy from CenterX_m
  if (!is.null(breath_col) && breath_col %in% names(points_df)) {
    states <- points_df[[breath_col]]
  } else {
    # presence proxy: 1 = detected (we assume detected at surface), 0 otherwise
    states <- ifelse(!is.na(points_df$CenterX_m), 1, 0)
  }
  states <- na.omit(states)
  n <- length(states)
  if (n <= 1) return(NA_real_)
  # compute fraction of dyads with equal state
  # Number equal dyads:
  tot_pairs <- choose(n, 2)
  eq_pairs <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (states[i] == states[j]) eq_pairs <- eq_pairs + 1
    }
  }
  return(eq_pairs / tot_pairs)
}

# 10. Synchrony Index (pairwise heading synchrony)
# Formula: Synchrony = Σ cos(θᵢ - θⱼ) / (n(n-1)/2)
# Rationale: calculates mean cosine of angular differences between all pairs.
# Biological meaning: values near 1 indicate high synchrony; can be calculated for different time windows.
# Reference: Proposed metrics from research framework PDF
compute_synchrony_index_frame <- function(points_df, angle_col = "MovingAvgAngle_deg") {
  if (!(angle_col %in% names(points_df))) return(NA_real_)
  angs <- points_df[[angle_col]] %>% na.omit() %>% as.numeric()
  n <- length(angs)
  if (n < 2) return(NA_real_)
  
  rads <- angs * pi / 180
  # Compute pairwise cosine of differences
  sync_sum <- 0
  count <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      sync_sum <- sync_sum + cos(rads[i] - rads[j])
      count <- count + 1
    }
  }
  return(sync_sum / count)
}

# 11. Circular Variance
# Formula: Circular_Variance = 1 - R̄, where R̄ = sqrt[(Σcos(θᵢ))² + (Σsin(θᵢ))²]/n
# Rationale: measures dispersion of heading angles; inverse of mean resultant length.
# Biological meaning: lower variance = higher coordination; 0 = perfect alignment, 1 = random directions.
# Reference: Proposed metrics from research framework PDF
compute_circular_variance_frame <- function(points_df, angle_col = "MovingAvgAngle_deg") {
  if (!(angle_col %in% names(points_df))) return(NA_real_)
  angs <- points_df[[angle_col]] %>% na.omit() %>% as.numeric()
  n <- length(angs)
  if (n == 0) return(NA_real_)
  
  rads <- angs * pi / 180
  R <- sqrt((sum(cos(rads)))^2 + (sum(sin(rads)))^2) / n
  return(1 - R)
}

# 12. Angular Velocity (mean turning rate per frame)
# Formula: Angular_Velocity = |θₜ - θₜ₋₁| / Δt
# Rationale: indicates turning rate and behavioral flexibility; requires ordered data per individual.
# Biological meaning: compare between habitats (open vs. confined); higher = more turning.
# Note: This computes per-individual angular changes within a frame if multiple angles exist.
# Reference: Proposed metrics from research framework PDF
compute_angular_velocity_frame <- function(points_df, angle_col = "MovingAvgAngle_deg", dt = 1) {
  if (!(angle_col %in% names(points_df))) return(NA_real_)
  angs <- points_df[[angle_col]] %>% na.omit() %>% as.numeric()
  n <- length(angs)
  if (n < 2) return(NA_real_)
  
  # Compute angular changes (handling wraparound)
  ang_diffs <- numeric(n-1)
  for (i in 1:(n-1)) {
    diff_raw <- abs(angs[i+1] - angs[i])
    ang_diffs[i] <- ifelse(diff_raw > 180, 360 - diff_raw, diff_raw)
  }
  return(mean(ang_diffs) / dt)
}

# ---------------------------
# Temporal (multi-frame) metric functions
# ---------------------------

# 13. Coordination Stability
# Formula: Standard deviation of polarization (MRL) over a time window
# Rationale: measures temporal stability of coordination patterns.
# Biological meaning: lower SD = more stable coordination patterns over time.
# Note: This is computed over a sequence of frames, not per-frame.
# Reference: Proposed metrics from research framework PDF
compute_coordination_stability <- function(frame_metrics, window_frames = 30) {
  # frame_metrics should be a dataframe with Heading_MRL column
  if (!"Heading_MRL" %in% names(frame_metrics)) return(NA_real_)
  if (nrow(frame_metrics) < window_frames) window_frames <- nrow(frame_metrics)
  if (window_frames < 2) return(NA_real_)
  
  # Use rolling window to compute stability
  mrl_vals <- frame_metrics$Heading_MRL
  n <- length(mrl_vals)
  stabilities <- numeric(n - window_frames + 1)
  
  for (i in 1:(n - window_frames + 1)) {
    window_vals <- mrl_vals[i:(i + window_frames - 1)]
    stabilities[i] <- sd(window_vals, na.rm = TRUE)
  }
  
  return(mean(stabilities, na.rm = TRUE))
}

# 14. Group Turning Events
# Formula: Identify when |mean(θₜ) - mean(θₜ₋₁)| > threshold
# Rationale: detects coordinated direction changes by the group.
# Biological meaning: frequency and magnitude of coordinated turns may reflect habitat complexity or foraging strategy.
# Note: This requires sequential frame data with group mean headings.
# Reference: Proposed metrics from research framework PDF
detect_turning_events <- function(frame_metrics, threshold_deg = 45) {
  # frame_metrics should have FrameIndex and heading data
  # We'll compute this from Heading_MRL direction (mean angle)
  # For simplicity, we'll use a proxy: when Heading_MRL drops significantly between frames
  
  if (nrow(frame_metrics) < 2) return(list(n_events = 0, event_frames = integer(0)))
  
  # We need actual mean heading angle, but MRL doesn't give us direction
  # This is a limitation - ideally we'd store mean heading angle
  # For now, we'll detect turns as sudden drops in MRL (indicating loss of coordination)
  
  mrl_vals <- frame_metrics$Heading_MRL
  mrl_changes <- abs(diff(mrl_vals))
  
  # Detect events where MRL changes more than a threshold (e.g., 0.3)
  mrl_threshold <- 0.3
  event_indices <- which(mrl_changes > mrl_threshold)
  event_frames <- frame_metrics$FrameIndex[event_indices + 1]  # +1 because diff shortens by 1
  
  return(list(
    n_events = length(event_frames),
    event_frames = event_frames,
    mean_mrl_change = mean(mrl_changes[event_indices], na.rm = TRUE)
  ))
}

# 15. Persistence in Direction (Autocorrelation)
# Formula: Calculate autocorrelation of heading angles over time
# Rationale: higher autocorrelation = more persistent movement direction.
# Biological meaning: may indicate different habitat navigation strategies.
# Note: Requires sequential heading data; uses circular correlation.
# Reference: Proposed metrics from research framework PDF
compute_direction_persistence <- function(frame_metrics, lag = 1, angle_col = "Heading_MRL") {
  # This is challenging with circular data
  # We'll use a simplified approach: correlation of MRL values over time
  # True circular autocorrelation would require individual heading angles
  
  if (!angle_col %in% names(frame_metrics)) return(NA_real_)
  if (nrow(frame_metrics) < lag + 2) return(NA_real_)
  
  vals <- frame_metrics[[angle_col]]
  vals <- vals[!is.na(vals)]
  
  if (length(vals) < lag + 2) return(NA_real_)
  
  # Simple lag-1 autocorrelation
  acf_result <- tryCatch(
    acf(vals, lag.max = lag, plot = FALSE)$acf[lag + 1],
    error = function(e) NA_real_
  )
  
  return(acf_result)
}

# ---------------------------
# Core driver functions
# ---------------------------

# compute_metrics_per_frame:
# given an interpolated df (FrameIndex, ObjectID, CenterX_m, CenterY_m, MovingAvgAngle_deg )
# returns a tibble with per-frame metrics
compute_metrics_per_frame <- function(df_interp,
                                      frame_col = "FrameIndex",
                                      id_col = "ObjectID",
                                      x_col = "CenterX_m",
                                      y_col = "CenterY_m",
                                      angle_col = "MovingAvgAngle_deg",
                                      breath_col = NULL,
                                      Nmin = 2) {
  # per-frame grouping
  frames <- sort(unique(df_interp[[frame_col]]))
  out <- vector("list", length(frames))
  k <- 1
  for (f in frames) {
    sub <- df_interp %>% filter(.data[[frame_col]] == f)
    # optionally filter only rows where x,y present
    pts <- sub %>% filter(!is.na(.data[[x_col]]) & !is.na(.data[[y_col]]))
    N <- nrow(pts)
    rec <- tibble(FrameIndex = f, N = N)
    if (N < Nmin) {
      rec <- rec %>% mutate(MPD = NA_real_, MNN = NA_real_, CHAI = NA_real_,
                            CentroidDist = NA_real_, StandardDist = NA_real_, 
                            DispersionIndex = NA_real_, MST_mean = NA_real_,
                            Heading_MRL = NA_real_, Heading_pairSim = NA_real_,
                            SynchronyIndex = NA_real_, CircularVariance = NA_real_,
                            AngularVelocity = NA_real_, BreathSync = NA_real_)
    } else {
      rec <- rec %>% mutate(
        MPD = compute_mpd_frame(pts),
        MNN = compute_mnn_frame(pts),
        CHAI = compute_chai_frame(pts),
        CentroidDist = compute_centroid_distance_frame(pts),
        StandardDist = compute_standard_distance_frame(pts),
        DispersionIndex = compute_dispersion_index_frame(pts),
        MST_mean = compute_mst_mean_frame(pts)
      )
      # heading metrics
      heading_res <- compute_heading_similarity_frame(pts, angle_col = angle_col)
      rec$Heading_MRL <- heading_res$MRL
      rec$Heading_pairSim <- heading_res$PairwiseSim
      rec$SynchronyIndex <- compute_synchrony_index_frame(pts, angle_col = angle_col)
      rec$CircularVariance <- compute_circular_variance_frame(pts, angle_col = angle_col)
      rec$AngularVelocity <- compute_angular_velocity_frame(pts, angle_col = angle_col)
      # breathing/surfacing
      rec$BreathSync <- compute_breathing_synchrony_frame(sub, breath_col = breath_col)
    }
    out[[k]] <- rec
    k <- k + 1
  }
  outdf <- bind_rows(out)
  return(outdf)
}

# sliding_window_summary:
# window_sec (seconds), 
# fps (frames per second of the original video)
# slide_step_frames (how many frames to jump between windows)
# summary_funs: which statistics to compute (defaults median for central tendency)
sliding_window_summary <- function(frame_metrics,
                                   fps = 30,
                                   window_sec = 5,
                                   slide_step_frames = NULL,
                                   stats = c("median", "sd", "mean"),
                                   whole_video = FALSE) {
  
  # OPTION 1: Compute metrics over the entire video
  if (whole_video) {
    sub <- frame_metrics   # the entire dataset
    
    return(tibble(
      WindowStart = min(sub$FrameIndex),
      WindowEnd   = max(sub$FrameIndex),
      N_med = median(sub$N, na.rm = TRUE),
      N_mean = mean(sub$N, na.rm = TRUE),
      N_sd = sd(sub$N, na.rm = TRUE),
      
      MPD_med = median(sub$MPD, na.rm = TRUE),
      MPD_mean = mean(sub$MPD, na.rm = TRUE),
      MPD_sd = sd(sub$MPD, na.rm = TRUE),
      
      MNN_med = median(sub$MNN, na.rm = TRUE),
      MNN_mean = mean(sub$MNN, na.rm = TRUE),
      MNN_sd = sd(sub$MNN, na.rm = TRUE),
      
      CHAI_med = median(sub$CHAI, na.rm = TRUE),
      CHAI_mean = mean(sub$CHAI, na.rm = TRUE),
      CHAI_sd = sd(sub$CHAI, na.rm = TRUE),
      
      CentroidDist_med = median(sub$CentroidDist, na.rm = TRUE),
      CentroidDist_mean = mean(sub$CentroidDist, na.rm = TRUE),
      CentroidDist_sd = sd(sub$CentroidDist, na.rm = TRUE),
      
      StandardDist_med = median(sub$StandardDist, na.rm = TRUE),
      StandardDist_mean = mean(sub$StandardDist, na.rm = TRUE),
      StandardDist_sd = sd(sub$StandardDist, na.rm = TRUE),
      
      DispersionIndex_med = median(sub$DispersionIndex, na.rm = TRUE),
      DispersionIndex_mean = mean(sub$DispersionIndex, na.rm = TRUE),
      DispersionIndex_sd = sd(sub$DispersionIndex, na.rm = TRUE),
      
      MST_med = median(sub$MST_mean, na.rm = TRUE),
      MST_mean = mean(sub$MST_mean, na.rm = TRUE),
      MST_sd = sd(sub$MST_mean, na.rm = TRUE),
      
      Heading_MRL_med = median(sub$Heading_MRL, na.rm = TRUE),
      Heading_MRL_mean = mean(sub$Heading_MRL, na.rm = TRUE),
      Heading_MRL_sd = sd(sub$Heading_MRL, na.rm = TRUE),
      
      SynchronyIndex_med = median(sub$SynchronyIndex, na.rm = TRUE),
      SynchronyIndex_mean = mean(sub$SynchronyIndex, na.rm = TRUE),
      SynchronyIndex_sd = sd(sub$SynchronyIndex, na.rm = TRUE),
      
      CircularVariance_med = median(sub$CircularVariance, na.rm = TRUE),
      CircularVariance_mean = mean(sub$CircularVariance, na.rm = TRUE),
      CircularVariance_sd = sd(sub$CircularVariance, na.rm = TRUE),
      
      AngularVelocity_med = median(sub$AngularVelocity, na.rm = TRUE),
      AngularVelocity_mean = mean(sub$AngularVelocity, na.rm = TRUE),
      AngularVelocity_sd = sd(sub$AngularVelocity, na.rm = TRUE),
      
      BreathSync_med = median(sub$BreathSync, na.rm = TRUE),
      BreathSync_mean = mean(sub$BreathSync, na.rm = TRUE),
      BreathSync_sd = sd(sub$BreathSync, na.rm = TRUE)
    ))
  }
  
  # OPTION 2: Sliding windows
  if (is.null(slide_step_frames))
    slide_step_frames <- max(1, floor((fps * window_sec) / 2))
  
  window_size <- fps * window_sec
  frame_ids <- sort(frame_metrics$FrameIndex)
  starts <- seq(min(frame_ids), max(frame_ids), by = slide_step_frames)
  
  out <- list()
  idx <- 1
  
  for (s in starts) {
    e <- s + window_size - 1
    sub <- frame_metrics %>% filter(FrameIndex >= s & FrameIndex <= e)
    
    if (nrow(sub) == 0) next
    
    rec <- tibble(
      WindowStart = s,
      WindowEnd = e,
      N_med = median(sub$N, na.rm = TRUE),
      N_mean = mean(sub$N, na.rm = TRUE),
      N_sd = sd(sub$N, na.rm = TRUE),
      
      MPD_med = median(sub$MPD, na.rm = TRUE),
      MPD_mean = mean(sub$MPD, na.rm = TRUE),
      MPD_sd = sd(sub$MPD, na.rm = TRUE),
      
      MNN_med = median(sub$MNN, na.rm = TRUE),
      MNN_mean = mean(sub$MNN, na.rm = TRUE),
      MNN_sd = sd(sub$MNN, na.rm = TRUE),
      
      CHAI_med = median(sub$CHAI, na.rm = TRUE),
      CHAI_mean = mean(sub$CHAI, na.rm = TRUE),
      CHAI_sd = sd(sub$CHAI, na.rm = TRUE),
      
      CentroidDist_med = median(sub$CentroidDist, na.rm = TRUE),
      CentroidDist_mean = mean(sub$CentroidDist, na.rm = TRUE),
      CentroidDist_sd = sd(sub$CentroidDist, na.rm = TRUE),
      
      StandardDist_med = median(sub$StandardDist, na.rm = TRUE),
      StandardDist_mean = mean(sub$StandardDist, na.rm = TRUE),
      StandardDist_sd = sd(sub$StandardDist, na.rm = TRUE),
      
      DispersionIndex_med = median(sub$DispersionIndex, na.rm = TRUE),
      DispersionIndex_mean = mean(sub$DispersionIndex, na.rm = TRUE),
      DispersionIndex_sd = sd(sub$DispersionIndex, na.rm = TRUE),
      
      MST_med = median(sub$MST_mean, na.rm = TRUE),
      MST_mean = mean(sub$MST_mean, na.rm = TRUE),
      MST_sd = sd(sub$MST_mean, na.rm = TRUE),
      
      Heading_MRL_med = median(sub$Heading_MRL, na.rm = TRUE),
      Heading_MRL_mean = mean(sub$Heading_MRL, na.rm = TRUE),
      Heading_MRL_sd = sd(sub$Heading_MRL, na.rm = TRUE),
      
      SynchronyIndex_med = median(sub$SynchronyIndex, na.rm = TRUE),
      SynchronyIndex_mean = mean(sub$SynchronyIndex, na.rm = TRUE),
      SynchronyIndex_sd = sd(sub$SynchronyIndex, na.rm = TRUE),
      
      CircularVariance_med = median(sub$CircularVariance, na.rm = TRUE),
      CircularVariance_mean = mean(sub$CircularVariance, na.rm = TRUE),
      CircularVariance_sd = sd(sub$CircularVariance, na.rm = TRUE),
      
      AngularVelocity_med = median(sub$AngularVelocity, na.rm = TRUE),
      AngularVelocity_mean = mean(sub$AngularVelocity, na.rm = TRUE),
      AngularVelocity_sd = sd(sub$AngularVelocity, na.rm = TRUE),
      
      BreathSync_med = median(sub$BreathSync, na.rm = TRUE),
      BreathSync_mean = mean(sub$BreathSync, na.rm = TRUE),
      BreathSync_sd = sd(sub$BreathSync, na.rm = TRUE)
    )
    
    out[[idx]] <- rec
    idx <- idx + 1
  }
  
  return(bind_rows(out))
}



# ---------------------------
# Plotting functions
# ---------------------------

# timeseries plot of one or more metrics
plot_timeseries_metrics <- function(frame_metrics, 
                                    metrics = c("MPD", "MNN", "CentroidDist"), 
                                    colours = NULL) {
  df_long <- frame_metrics %>% select(FrameIndex, all_of(metrics)) %>%
    pivot_longer(cols = -FrameIndex, names_to = "metric", values_to = "value")
  p <- ggplot(df_long, aes(x = FrameIndex, y = value, color = metric)) +
    geom_line() + theme_minimal() +
    labs(x = "FrameIndex", y = "Value", title = paste("Time series:", paste(metrics, collapse = ", ")))
  return(p)
}

# convex hull plot for a specific frame
plot_convex_hull_frame <- function(df_interp, frame_id, x_col = "CenterX_m", y_col = "CenterY_m") {
  sub <- df_interp %>% filter(FrameIndex == frame_id) %>% drop_na(.data[[x_col]], .data[[y_col]])
  if (nrow(sub) < 1) {
    stop("No detections for that frame")
  }
  p <- ggplot(sub, aes_string(x = x_col, y = y_col, label = "ObjectID")) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 3) +
    coord_equal() + theme_minimal() + labs(title = paste("Frame", frame_id))
  if (nrow(sub) >= 3) {
    idx <- chull(sub[[x_col]], sub[[y_col]])
    hull_df <- sub[c(idx, idx[1]), ]
    p <- p + geom_polygon(data = hull_df, aes_string(x = x_col, y = y_col), alpha = 0.2, fill = "steelblue", color = "black")
  }
  return(p)
}

# pairwise presence overlap heatmap
plot_pair_overlap_heatmap <- function(df_interp, id_col = "ObjectID", frame_col = "FrameIndex", x_col = "CenterX_m") {
  M <- frames_presence_matrix(df_interp, id_col = id_col, frame_col = frame_col, x_col = x_col)
  if (nrow(M) < 2) {
    stop("Need at least 2 tracked individuals to plot pairwise overlap.")
  }
  n <- nrow(M)
  overlap <- matrix(NA_real_, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      a <- M[i,]; b <- M[j,]
      denom <- sum((a==1) | (b==1))
      overlap[i, j] <- ifelse(denom == 0, NA_real_, sum((a==1)&(b==1)) / denom)
    }
  }
  rownames(overlap) <- rownames(M)
  colnames(overlap) <- rownames(M)
  mm <- reshape2::melt(overlap, varnames = c("A","B"), value.name = "Overlap")
  p <- ggplot(mm, aes(x = A, y = B, fill = Overlap)) +
    geom_tile() +
    scale_fill_viridis(option = "C", na.value = "grey90") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Pairwise Presence Overlap (Jaccard-like)", x = "Individual", y = "Individual")
  return(p)
}

# ---------------------------
# Master function that runs everything
# ---------------------------
# compute_group_metrics:
# - csv_path: path to input csv
# - fps: frames per second of the video (check with exiftool)
# - sliding = TRUE/FALSE: whether to compute sliding window summaries 
# - window_sec:  length of the sliding window, in seconds: should the metrics be averaged over this window?
# - whole_video = TRUE/FALSE: should the metrics be averaged by the entire video (ie, whole csv)?
# - interp_gap: pass to interpolate_tracks(); should gaps in detection be interpolated, and by how many frame?
# - Nmin: minimum number of detections in the frame to compute the metrics; default is 2 individuals, it will be NA for a single detection
# - angle_col, breath_col: column names if present
# frame_selection: input: a vector with 2 values, min and max frame number. E.g. c(0, 100) to analyse only the first 100 frames of the video (ie, only the 100 rows of the csv file)
# return_plots: FALSE/TRUE

compute_group_metrics <- function(csv_path,
                                  frame_selection = NA,
                                  fps = 30,
                                  window_sec = 5,
                                  whole_video = FALSE,
                                  interp_gap = 1,
                                  Nmin = 2,
                                  angle_col = "MovingAvgAngle_deg",
                                  breath_col = NULL,
                                  sliding = TRUE,
                                  slide_step_frames = NULL,
                                  return_plots = F) {
  # Read
  df_raw <- read_csv(csv_path)
  req_cols <- c("FrameIndex", "ObjectID", "CenterX_m", "CenterY_m")
  missing <- setdiff(req_cols, names(df_raw))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  
  # apply frame restrictions (rows)
  if(!is.na(frame_selection[1])){
    df_raw = df_raw[df_raw$FrameIndex>frame_selection[1] &
                            df_raw$FrameIndex<frame_selection[2], ]
  } 

  # Interpolate short gaps
  df_interp <- interpolate_tracks(df_raw,
                                  frame_col = "FrameIndex",
                                  id_col = "ObjectID",
                                  x_col = "CenterX_m",
                                  y_col = "CenterY_m",
                                  angle_col = angle_col,
                                  maxgap = interp_gap)
  
  # Per-frame metrics
  frame_metrics <- compute_metrics_per_frame(df_interp,
                                             frame_col = "FrameIndex",
                                             id_col = "ObjectID",
                                             x_col = "CenterX_m",
                                             y_col = "CenterY_m",
                                             angle_col = angle_col,
                                             breath_col = breath_col,
                                             Nmin = Nmin)
  
  # sliding window
  window_metrics <- if (sliding) sliding_window_summary(frame_metrics, 
                                                        fps = fps, 
                                                        window_sec = window_sec, 
                                                        slide_step_frames = slide_step_frames,
                                                        whole_video = whole_video) else NULL
  
  # pairwise presence overlap median (synchrony proxy)
  pres_mat <- frames_presence_matrix(df_interp, 
                                     id_col = "ObjectID", 
                                     frame_col = "FrameIndex", 
                                     x_col = "CenterX_m")
  pair_overlaps <- c()
  if (nrow(pres_mat) >= 2) {
    for (i in 1:(nrow(pres_mat)-1)) {
      for (j in (i+1):nrow(pres_mat)) {
        a <- pres_mat[i,]; b <- pres_mat[j,]
        denom <- sum((a==1)|(b==1))
        pair_overlaps <- c(pair_overlaps, ifelse(denom==0, NA_real_, sum((a==1)&(b==1))/denom))
      }
    }
  }
  pair_overlap_median <- ifelse(length(pair_overlaps)>0, median(pair_overlaps, na.rm = TRUE), NA_real_)
  
  # produce plots if requested
  plots <- list()
  if (return_plots) {
    plots$mpd_ts <- plot_timeseries_metrics(frame_metrics, metrics = c("MPD"))
    plots$mnn_ts <- plot_timeseries_metrics(frame_metrics, metrics = c("MNN"))
    plots$head_mrl_ts <- plot_timeseries_metrics(frame_metrics, metrics = c("Heading_MRL"))
    plots$window_mpd <- if (!is.null(window_metrics)) ggplot(window_metrics, aes(x = WindowStart, y = MPD_med)) + geom_line() + theme_minimal() + labs(title = paste0("Windowed median MPD (", window_sec, " s)"), x = "WindowStart", y = "MPD_m") else NULL
    # generate an example convex hull for the median frame (if any)
    if (nrow(frame_metrics) > 0) {
      med_frame <- median(frame_metrics$FrameIndex, na.rm = TRUE)
      plots$convex_hull_example <- tryCatch(plot_convex_hull_frame(df_interp, frame_id = med_frame), error = function(e) NULL)
    }
    # pair overlap heatmap
    plots$pair_overlap_heatmap <- tryCatch(plot_pair_overlap_heatmap(df_interp), error = function(e) NULL)
  }
  
  return(list(
    df_interp = df_interp,
    frame_metrics = frame_metrics,
    window_metrics = window_metrics,
    pair_overlap_median = pair_overlap_median,
    plots = plots
  ))
}
