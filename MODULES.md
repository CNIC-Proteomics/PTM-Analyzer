# Modules

1. [NMpyCompare](#1-nmpycompare)
2. [ReportLimma](#2-reportlimma)
3. [FDRoptimizer](#3-fdroptimizer)
4. [PTMMap](#4-ptmmap)
5. [qTableReport](#5-qtablereport)

### 1. NMpyCompare

This module calculates NM-corrected values from the iSanXoT report. It subtracts the Zpgm2p value of the associated non-modified peptidoform from the Zpgm2p value of each modified peptidoform. Users can specify column names for integration levels, quantitative values, and criteria for identifying non-modified peptidoforms through the configuration file. The NM-corrected values are then appended as new columns in the iSanXoT report.

### 2. ReportLimma

This script performs hypothesis testing for comparisons between different groups across integration levels using the **limma** package. It computes p-values for statistical contrasts, the mean difference between groups, and the LPS value (-log(p-value) adjusted by the sign of the mean difference).

### 3. FDRoptimizer

Implemented in Python, this script applies an iterative algorithm to optimize the spectral count threshold at each integration step. It scans different thresholds to maximize the number of elements with a q-value below a user-defined threshold, using the Benjamini-Hochberg algorithm for multiple testing correction.

### 4. PTMMap

PTMMap is a tool developed with the aim of visualizing, interpreting, and comparing the proteins' PTMs. This module represents as many maps as proteins for which any integration meets the threshold established by the user. Each map illustrates the change between one condition and another based on the p-value of all calculated integrations, on the y-axis.

```
LPS = -log2(p-value) * sign (Condition2 - Condition1)
```

On the x-axis, the position of each residue of the protein is represented. Specific modifications and hypermodified zones are represented by circles, while partial and total digestion, and zonal changes are represented by rectangles. The size of these markers depends on the frequency of each parameter, in a relative scale depending on the maximum and minimum PSMs frequency of each type of modification. These graphs offer interactivity and enable the visualization of parameter frequency, modified residue, and the specific group of each Δmass.

- **Input:**
  - `.tsv` file with limma p-values and -LPS calculated for each integration.
  - Configuration file (`.ini`): There is a default `.ini` in the “config” folder:
    - PTMMap parameters:
      ```
      pgm_column_name: pgm column name.
      g_column_name: group column name.
      a_column_name: modified residue column name.
      n_column_name: modified position within the protein column name.
      e_column_name: position of the last residue of the peptide within the protein column name.
      p_column_name: peptide column name.
      q_column_name: protein column name.
      d_column_name: Δmass column name.
      qc_column_name: cluster containing overlapping peptides column name.
      pFreq_column_name: peptide frequency column name.
      qcFreq_column_name: qc frequency column name.
      pgmFreq_column_name: pgm frequency column name.
      first_b_column_name: position in protein of the first residue of the peptide.
      description_column_name: description of the protein column name.
      Missing_Cleavages_column_name: number of missing cleavages column name.
      LPS_p2qc_column_name: -LPS of p2qc integration column name.
      LPS_qc2q_column_name: -LPS of qc2q integration column name.
      LPS_pgm2p_column_name: -LPS of pgm2p integration column name.
      LPS_pgm2p_NM_column_name: -LPS of pgm2p integration only with the unmodified peptides column name.
      Filter_pgm2p_NM_column_name: filtering column (FDR, p-value) of pgm2p integration corrected by NMCompare.
      Filter_pgm2p_column_name: filtering column (FDR, p-value) of pgm2p integration calculated only with NM.
      Filter_p2qc_column_name: filtering column (FDR, p-value) of p2qc integration.
      Filter_qc2q_column_name: filtering column (FDR, p-value) of qc2q integration.
      threshold_pgm2p_NM: threshold of filtering column (pgm2p integration only with unmodified peptides).
      threshold_pgm2p: threshold of filtering column (pgm2p integration).
      threshold_p2qc: threshold of filtering column (p2qc integration).
      threshold_qc2q: threshold of filtering column (qc2q integration).
      NM: non-modified peptides group.
      path_plots_with_threshold: folder in which filtered PTM Maps will be saved.
      path_plots_Without_threshold: folder in which PTM maps, without filters, will be saved.
      ```

- **Output:**
  - Two output folders. Both contain only maps of the proteins for which some of their modifications meet the threshold established by the user. In one of the folders (`path_plots_with_threshold`), only the maps featuring modifications that meet the threshold set by the user are represented, while in the other (`path_plots_Without_threshold`), complete maps of the proteins are depicted.


[Go to top](#modules)
___


### 5. qTableReport

This module enables a detailed exploration of significant changes at the protein level in a peptide-centric workflow. It generates an output table summarizing the number of modified and non-modified peptidoforms with significant increases or decreases, along with details on digestion status and qc clusters.


[Go to top](#modules)
___

