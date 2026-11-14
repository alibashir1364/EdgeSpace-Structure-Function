# Structural Edge Context Governs the Spectral and Event-Based Timescales of Functional Brain Interactions

Official MATLAB implementation for the manuscript:
“Structural edge context governs the spectral and event-based timescales of functional brain interactions”
*Ali Bashirgonbadi, Mohamad Reza Salehi, and Hamid Soltanian-Zadeh*
---
## Overview
While functional interactions occur across the complete pairwise space of cortical regions (Nedge = N_Node(N_Node-1)/2), diffusion MRI tractography reconstructs 
direct white-matter pathways for only a sparse subset of pairs. This repository provides a parameter-free, edge-centric framework that bridges this topological 
mismatch by projecting the human structural connectome into a **complete weighted line graph** (L(KN)).
Every cortical region pair (i, j)—both directly connected (wij > 0) and non-directly connected (wij = 0)—is assigned a principled anatomical descriptor termed 
**Structural Edge Context**:
C_ij^S  = s_i+ s_j-w_ij
Where s_i denotes regional structural connectivity strength. We couple C_ij^S and focal edge weight (w_ij) with functional **Edge Time Series (ETS)** to interrogate 
two complementary dimensions of edge-level dynamics:
1. **Continuous Spectral Timescales:** Parameter-free **Median Spectral Frequency (MSF)** derived from Welch power spectral density (PSD) estimation.
2. **Discrete Co-fluctuation Events:** Supra-threshold **Event Duration** and **Event Amplitude** across the complete interaction space.

## System Requirements & Dependencies

**MATLAB** (tested on R2021a and later)
---
## Required Input Data Format

Before running the primary analysis pipeline, ensure that your preprocessed neuroimaging data (e.g., from the Human Connectome Project, N_subject = 100, TR = 0.72 sec) are formatted as 
`.mat` files containing the following arrays:

1. **Resting-State fMRI Time Series (`roiTS_all`):**
   **Dimensions:** `N_region × N_timepoint × N_subject` (e.g., `200 × 1200 × 100` for Schaefer-200, `68 × 1200 × 100` for Desikan–Killiany DK-68, or `360 × 1200 × 100` for Glasser-360).
   **Description:** Preprocessed region-wise BOLD signals across recording frames and participants.

2. **Subject-Level Structural Connectivity (`SC_all`):**
   **Dimensions:** `N_region × N_region × N_subject` (e.g., `200 × 200 × 100` or `68 × 68 × 100`).
   **Description:** Symmetric, zero-diagonal structural connectivity matrices (e.g., SIFT2-weighted streamline counts).

3. **Group-Level Structural Connectome (`SC_group`):**
   **Dimensions:** `N_region × N_region`
   **Description:** Group-representative structural connectivity matrix (constructed via distance-dependent consensus thresholding or group averaging inside the main script, with zero diagonal and symmetric weights w_ij).
---
## Execution Workflow

The analysis pipeline is modular. **Step 1** executes the core structural line-graph construction and Welch spectral analysis, generating a comprehensive master 
`.mat` file that serves as the input for all downstream statistical, null-model, event-based, and visualization scripts.

### Step 1: Run the Core Edge-Space Structure–Function Pipeline
Open and run **`Main_EdgeSpace_Structure_Function_Analysis.m`**.

* **Configuration:** In `STEP 0 — USER SETTINGS` at the top of the script, specify your target parcellation and input file names:
  ```matlab
  fMRI_file   = 'ROI_TS_68_HCP_100subj.mat';       % or 'ROI_TS_shefer_200_HCP_100subj.mat'
  SC_file     = 'SC_SIFT2_68_100subj.mat';         % or 200-region SC file
  output_file = 'EdgeSpace_Structure_Function_results_68regions.mat'; % or 'EdgeSpace_Structure_Function_results.mat'
  N_nodes     = 68;  % Set to 200 for Schaefer-200, 68 for DK-68, or 360 for Glasser-360
  N_subs      = 100;
  N_time      = 1200;
  TR          = 0.72;
