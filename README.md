# Dolphin Group Dynamics: Automated Metrics from Drone Footage

> **Quantifying collective behavior and social organization in dolphin groups through automated analysis of drone-based video tracking data.**

---

## Overview

This repository provides a modular computational pipeline for extracting **group-level behavioral metrics** from automated dolphin detections in drone footage. The toolkit quantifies multiple dimensions of collective behavior including:

- **Spatial structure** (cohesion, spread, connectivity)
- **Directional coordination** (heading alignment, synchrony)
- **Temporal dynamics** (activity synchronization, state transitions)

These metrics enable:
- Compare social organization across populations and habitats
- Identify behavioral states from aerial patterns
- Quantify environmental influences on group dynamics
- Detect temporal patterns in collective behavior


---

## Metrics Computed

The pipeline computes **7 core metrics** that capture complementary aspects of group organization:

### 1. **Mean Pairwise Distance (MPD)**
- **Formula:** `MPD = mean(all pairwise distances)`
- **Measures:** Overall group compactness
- **Interpretation:** 
  - Low MPD → Tight, cohesive group (coordinated foraging, travel)
  - High MPD → Dispersed group (exploration, low interdependence)

### 2. **Median Nearest Neighbor Distance (MNN)**
- **Formula:** `MNN = median(min distance to any other individual)`
- **Measures:** Local social spacing (robust to outliers)
- **Interpretation:** Typical spacing between closest companions

### 3. **Convex Hull Area per Individual (CHAI)**
- **Formula:** `CHAI = convex_hull_area / N`
- **Measures:** Group footprint normalized by size
- **Interpretation:**
  - Small CHAI → Concentrated formation (vigilance, herding)
  - Large CHAI → Spread formation (searching, prey pursuit)

### 4. **Mean Distance to Centroid (MDC)**
- **Formula:** `MDC = mean(distance from each individual to group center)`
- **Measures:** Radial dispersion from group center
- **Interpretation:** Reveals core-periphery structure and centralization

### 5. **Minimum Spanning Tree Mean Edge Length (MST)**
- **Formula:** `MST = mean(edge lengths in minimum spanning tree)`
- **Measures:** Topological connectivity (robust to outliers)
- **Interpretation:**
  - Short MST → Functionally connected, coordinating
  - Long MST → Fragmented, potential fission

### 6. **Heading Similarity (MRL & Pairwise)**
- **Formula:** 
  - `MRL = ||Σ(unit heading vectors)|| / N`  (range: 0-1)
  - `Pairwise_Sim = 1 - (mean_angular_diff / π)`  (range: 0-1)
- **Measures:** Directional coordination
- **Interpretation:**
  - High alignment → Coordinated travel, herding
  - Low alignment → Milling, exploratory divergence

### 7. **Breathing/Surfacing Synchrony**
- **Formula:** `Synchrony = fraction of dyads sharing same surfacing state`
- **Measures:** Temporal coupling of activity cycles
- **Interpretation:** Proxy for social bonds, shared vigilance, coordinated diving


---


## Documentation

For detailed mathematical formulations and biological justification of each metric, see:

- **`Dolphin-group-metrics.docx`** - Complete documentation with formulas and rationale


---

**Last Updated:** December 2024
**Version:** 1.0.0
**Tested on:** R 4.3.0+ (Mac, Linux, Windows)
