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
  - Configuration file (`.yml`)
    - PTMMap parameters:
      ```
      PTMMap:
        # Folder in which filtered PTM Maps will be saved
        path_plots_with_threshold: PTMmaps_FDR
        path_plots_Without_threshold: PTMmaps

        # Required column names
        pgm_column_name: [pgm, LEVEL]
        g_column_name: [g, REL]
        a_column_name: [a, REL]
        n_column_name: [n, REL]
        e_column_name: [e, REL]
        p_column_name: [p, LEVEL]
        q_column_name: [q, LEVEL]
        d_column_name: [d, REL]
        qf_column_name: [qf, LEVEL]
        pFreq_column_name: [pFreq, REL]
        qfFreq_column_name: [qfFreq, REL]
        pgmFreq_column_name: [pgmFreq, REL]
        first_b_column_name: [first_b, REL]
        description_column_name: [description, REL]
        Missing_Cleavages_column_name: [Missing_Cleavage, REL]

        # How Non-modified are named
        NM: NM

        # Group to be analysed
        groups:
          - H-C

        # LPS integrations for: p2qf, qf2q, pfm2p, pgm2p_MN
        # [LowLevel_firstRow, LowLevel_secondRow]
        LPS_ColumnNames:
          - [Z_p2qf_logLimma, LPS]
          - [Z_qf2q_logLimma, LPS]
          - [Z_pgm2p_logLimma, LPS]
          - [Z_pgm2p_dNM_logLimma, LPS]

        # NM comparison for pgm2p integrations
        # [LowLevel_firstRow, LowLevel_secondRow]
        NM_ColumnNames:
          - [Z_pgm2p_dNM_limma, qvalue]
          - [Z_pgm2p_limma_NM_ONLY, qvalue]

        # Filter integrations for: p2qf, qf2q
        # [LowLevel_firstRow, LowLevel_secondRow]
        Filter_ColumnNames:
          - [Z_p2qf_limma, qvalue]
          - [Z_qf2q_limma, qvalue]

        # Threshold of filtering for the given integrations: pgm2p_NM, pgm2p, p2qf, qf2q
        threshold_pgm2p_NM: 0.05
        threshold_pgm2p: 0.05
        threshold_p2qf: 0.05
        threshold_qf2q: 0.05
        pgmFreqThreshold: 0
      ```

- **Output:**
  - Two output folders. Both contain only maps of the proteins for which some of their modifications meet the threshold established by the user. In one of the folders (`path_plots_with_threshold`), only the maps featuring modifications that meet the threshold set by the user are represented, while in the other (`path_plots_Without_threshold`), complete maps of the proteins are depicted.


[Go to top](#modules)
___


### 5. qTableReport

This module enables a detailed exploration of significant changes at the protein level in a peptide-centric workflow. It generates an output table summarizing the number of modified and non-modified peptidoforms with significant increases or decreases, along with details on digestion status and qc clusters.


[Go to top](#modules)
___

