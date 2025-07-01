# Modules

1. [MergeiSanxotPDM](#1-mergeisanxotpdm)
2. [NMpyCompare](#2-nmpycompare)
3. [ReportLimma](#3-reportlimma)
4. [FDRoptimizer](#4-fdroptimizer)
5. [PTMMap](#5-ptmmap)
6. [qReportMaker](#6-qreportmaker)



### 1. MergeiSanxotPDM
This Python script merges the iSanXoT report of protein2proteinall with the PDMTable at PGM level.




### 2. NMpyCompare

This script computes NM-corrected quantitative values from **iSanXoT** [1] reports. It operates by identifying a non-modified (NM) peptidoform within each cluster based on user-defined criteria specified in a configuration file. By default, it uses pgm2p integration, with "pgm" as the lower level and "p" as the higher level.

Once the NM reference is determined, the module appends its associated Z integration values to the report and calculates the difference (ΔZ) between these values and those of the modified peptidoforms for each cluster. The correction is performed per sample and added as new columns to the dataset, labeled as *Z_<lowerLevel>2<higherLevel>_dNM*. In cases where no NM reference is found for a cluster, the original Z value is retained.

NMpyCompare is configurable to accommodate various integration levels and quantification schemes, providing a robust tool for normalization in PTM proteomics workflows.

- **Inputs**
  * `--infile` (`-i`): Path to the input **iSanXoT report** file (in `.tsv` format) with multi-indexed levels.
  * `--config` (`-c`): Path to the **YAML configuration file** *(default: `NMpyCompare.yaml` located in the script directory)* specifying:
  ```
      NMpyCompare:
        # Column name of low level (FirstRowName#SecondRowName)
        low_level: pgm#LEVEL
        # Column name of high level (FirstRowName#SecondRowName)
        high_level: p#LEVEL
        # Column name of integration
        integration: Z_pgm2p
        # ScanFreq
        scanfreq: pgmFreq#REL
        # Columns considered to find non modified (FirstRowName#SecondRowName)
        NM_columns:
          - g#REL:
            - NM
    ```

- **Outputs**
  * `--outfile` (`-o`): Path to the output file, which is a modified iSanXoT report containing additional columns with NM reference values and ΔZ values for each modified peptidoform

[Go to top](#modules)
___




### 3. ReportLimma

Implemented in R, this script performs statistical hypothesis testing for quantitative comparisons across experimental groups using the **limma** [2] package. It operates on integration-level data derived from **iSanXoT** reports and uses a user-defined sample table and configuration file to specify contrasts and analysis parameters.

For each specified contrast, the script computes:
  * **p-values** for differential abundance,
  * **mean differences** between groups,
  * and **LPS (log p-value signed)** scores, defined as the negative logarithm of the p-value multiplied by the sign of the mean difference:
    *LPS = -log10(p-value) × sign(mean difference)*.

This unified metric facilitates interpretation of both statistical significance and direction of change.

- **Inputs**
  * `--input` (`-i`): Path to the **iSanXoT report** file (required). Should contain integration-level quantitative data (e.g., Z-values).
  * `--samples` (`-s`): Path to the **sample table** file (required). Must define sample identifiers and associated experimental groupings for contrast analysis.
  ```
    Control Heteroplasmy
    C1	        H1
    C2	        H2
    C3	        H3
      	        H4
  ```
* `--config` (`-c`): Path to the **YAML configuration file** (required). Specifies:
    ```
    # General parameters
    General:

      # Comparison groups
      groups:
        - [ H, C ]

      # Significance value used across multiple modules. Can be either a qvalue or pvalue. Defaults to pvalue.
      significance_value: 'pvalue'


    # Specific parameters of module
    LimmaCompare:

      # Possible types 
      # - limma
      # - limma_with_duplicates
      # - t-test
      test_type:
        - limma
        - t-test

      # Integrations used to apply Limma
      # [ [LowLevel_firstRow, LowLevel_secondRow], integration ]
      ColumnNames:
        - [ [pgm, LEVEL], Z_pgm2p]
        - [ [pgm, LEVEL], Z_pgm2p_dNM]
        - [ [p, LEVEL], Z_p2qc]
        - [ [qc, LEVEL], Z_qc2q]
        - [ [q, LEVEL], Z_q2all]
    ```

- **Outputs**
  * `--output` (`-o`): Path to the **output file** (required). The output includes original data with additional columns for:
    * Mean difference
    * p-value
    * LPS score for each feature across specified contrasts.

[Go to top](#modules)
___




### 4. FDRoptimizer

**FDRoptimizer** is a Python script designed to optimize spectral count thresholds in integration-level statistical reports-such as those generated from **limma** comparisons in **iSanXoT** reports. The goal of this tool is to maximize statistical significance by scanning a range of spectral count thresholds and selecting the value that yields the highest number of elements below a user-defined false discovery rate (FDR) threshold.

For each group comparison and integration level, the script:

  * Iteratively scans a specified range of minimum spectral count thresholds.
  * Applies the **Benjamini-Hochberg procedure** to compute FDR-adjusted q-values.
  * Identifies the threshold that results in the **maximum number of significant elements** (i.e., those with q-value ≤ specified threshold).
  * Optionally, adds the computed q-values to the original dataset for downstream analysis.

Additionally, FDRoptimizer generates interactive visualizations using **Plotly**, illustrating how the number of significant elements varies with increasing spectral count thresholds. These HTML reports are saved in an output directory to facilitate interpretation.

- **Inputs**
  * `--infile` (`-i`): Path to the input report containing p-values from limma comparisons.
  * `--config` (`-c`): Path to the YAML configuration file specifying:
  ```
    # General parameters
    General:

      # Comparison groups
      groups:
        - [ H, C ]

      # Significance value used across multiple modules. Can be either a qvalue or pvalue. Defaults to pvalue.
      significance_value: 'pvalue'
      

    # Specific parameters of module
    FDRoptimizer:

      # Column 
      # [ [LowLevel_firstRow, LowLevel_secondRow] , [Freqs_firstRow, Freqs_secondRow], Integration]
      ColumnNames:
        - [ [pgm, LEVEL], [pgmFreq, REL], Z_pgm2p_limma-NM(g&REL) ]
        - [ [pgm, LEVEL], [pgmFreq, REL], Z_pgm2p_dNM_limma ]
        - [ [p, LEVEL], [pFreq, REL], Z_p2qc_limma ]
        - [ [qc, LEVEL], [qcFreq, REL], Z_qc2q_limma ]
        - [ [q, LEVEL], [qFreq, REL], Z_q2all_limma ]

      # FDR thresholds applied
      FDR_Thr:
        - 0.01
        - 0.05

      # Scan frequency window
      Window: [0, 50]

      # Add column with FDR of at maximum value
      AddFDR: True
  ```

- **Outputs**
  * `--outfile` (`-o`): Output path to save the updated report with optional q-value annotations:
    * A modified report file (if `AddFDR: true`) with new q-value columns for the optimal threshold.
    * HTML plots showing the number of significant elements vs. scan frequency thresholds for each group and integration.
    * Log file tracking script execution and key statistics.

[Go to top](#modules)
___




### 5. PTMMap

PTMMap is a visualization tool designed to interpret and compare post-translational modifications (PTMs) across proteins. For every protein where at least one integration meets a user-defined threshold, the module generates a corresponding PTM map. Each map visualizes the change between two conditions based on the p-values of all integration events, displayed on the y-axis as LPS (-log2(p-value) * sign (Condition2 - Condition1)).

The x-axis represents the position of each residue in the protein. Specific PTMs and hypermodified regions are shown as circles, while events such as partial or complete digestion and broader zonal changes are represented as rectangles. The size of each marker corresponds to the relative frequency of the associated feature, scaled according to the minimum and maximum PSM frequencies for each modification type. These interactive plots allow users to explore the frequency of each parameter, the modified residue, and the group associated with each Δmass.

- **Inputs**
  * `--infile` (`-i`): Path to the input report containing p-values from limma comparisons and LPS scores for each integration.
  * `--config` (`-c`): Path to the YAML configuration file specifying:
    ```
    # General parameters
    General:

      # Comparison groups
      groups:
        - [ H, C ]

      # Significance value used across multiple modules. Can be either a qvalue or pvalue. Defaults to pvalue.
      significance_value: 'pvalue'
      

    # Specific parameters of module
    PTMMap:

      # Folder in which filtered PTM Maps will be saved
      path_plots_with_threshold: PTMmaps_filtered
      path_plots_Without_threshold: PTMmaps

      # Plot parameters
      font_size: 50
      grid: 'No'
      plot_width: 1700
      plot_height: 850

      # Required column names
      pgm_column_name: [pgm, LEVEL]
      g_column_name: [g, REL]
      a_column_name: [a, REL]
      n_column_name: [n, REL]
      e_column_name: [e, REL]
      p_column_name: [p, LEVEL]
      q_column_name: [q, LEVEL]
      d_column_name: [d, REL]
      qc_column_name: [qc, LEVEL]
      pFreq_column_name: [pFreq, REL]
      qcFreq_column_name: [qcFreq, REL]
      pgmFreq_column_name: [pgmFreq, REL]
      first_b_column_name: [first_b, REL]
      description_column_name: [description, REL]
      Missing_Cleavages_column_name: [Missing_Cleavage, REL]

      # How Non-modified are named
      NM: NM

      # LPS integrations for: p2qc, qc2q, pgm2p, pgm2p_MN
      # [LowLevel_firstRow, LowLevel_secondRow]
      LPS_ColumnNames:
        p2qc:     [Z_p2qc_logLimma, LPS]
        qc2q:     [Z_qc2q_logLimma, LPS]
        pgm2p:    [Z_pgm2p_logLimma, LPS]
        pgm2p_NM: [Z_pgm2p_dNM_logLimma, LPS]

      # NM comparison for pgm2p integrations
      # LowLevel_firstRow
      NM_ColumnNames:
        pgm2p:    Z_pgm2p_dNM_limma
        pgm2p_NM: Z_pgm2p_limma_NM_ONLY

      # Filter integrations for: p2qc, qc2q
      # LowLevel_firstRow
      Filter_ColumnNames:
        p2qc: Z_p2qc_limma
        qc2q: Z_qc2q_limma

      # Threshold of filtering for the given integrations: pgm2p_NM, pgm2p, p2qc, qc2q
      threshold_pgm2p_NM: 0.05
      threshold_pgm2p: 0.05
      threshold_p2qc: 0.05
      threshold_qc2q: 0.05
      pgmFreqThreshold: 0
    ```

- **Outputs**
  * `--outdir` (`-o`): Directory where the protein maps will be saved.
    * The `path_plots_with_threshold` folder includes only maps where at least one modification meets the specified threshold.
    * The `path_plots_Without_threshold` folder contains complete maps for all relevant proteins, regardless of threshold criteria.

[Go to top](#modules)
___




### 6. qReportMaker

This module enables in-depth quantitative and statistical analysis of protein-level changes in peptide-centric proteomics workflows. It processes peptide-level differential expression data, aggregates it at the protein level, and generates detailed reports that highlight proteins with significantly increasing or decreasing post-translational modifications (PTMs) and non-modified (NM) peptidoforms. The module is designed to support quality control and biological interpretation by integrating digestion characteristics and quantification clusters (QC clusters).

The module generates:
  * **qReports**: Excel files summarizing PTM/NM statistics per protein, with separate files for:
    - Upregulated proteins,
    - Downregulated proteins,
    - Combined up/down reports

    Each qReport file includes:
      - Total PTMs, NMs, and digestion-specific metrics,
      - directional statistics (up/down),
      - enrichment p-values via hypergeometric testing, and
      - sortable pivot tables with configurable PTM ordering.

  * **Basal Reports**: TSV files showing total PTM and NM counts across proteins, before significance filtering.
  * **Frequency Tables**: PTM-specific summary tables used to support statistical enrichment calculations.
  * **Links to PTM maps** directly embedded in Excel outputs.



- **Inputs**
  * `--infile` (`-i`): Path to tabular report with peptide-level differential expression results and annotations.
  * `--config` (`-c`): Path to the YAML configuration file specifying:
    ```
    # General parameters
    General:

      # Comparison groups
      groups:
        - [ H, C ]

      # Significance value used across multiple modules. Can be either a qvalue or pvalue. Defaults to pvalue.
      significance_value: 'pvalue'
      

    # Specific parameters of module
    qReportMaker:

      # Number of works in parallel. Othwerise is the 75% of total number of cpu's
      # n_cpu: 8

      # Folder where the outputs will be saved
      outDirName: qReports

      # Folder names of PTMMaps
      path_plots_with_threshold: PTMmaps_filtered
      path_plots_Without_threshold: PTMmaps

      # FDR threshold used
      qvThr: [0.01, 0.05, 0.1, 1]

      #
      # Columns information
      #

      # First row column name
      # Second row column name

      # Column name containing group
      gCol: [g, REL]


      # Column name containing modified aminoacid
      aCol: [a, REL]
          
      # Column name containing peptide position of modification
      mCol: [m, REL]

      # Name of the group corresponding to non-modified
      NMgroup: NM


      #
      # LEVELS
      #

      # Column containing pgm in the following format: PEP[MOD]TIDE
      pdmCol: [pgm, LEVEL]

      # Column containing p
      pCol: [p, LEVEL]

      # Column containing qc
      qcCol: [qc, LEVEL]

      # Column containing q
      qCol: [q, LEVEL]

      # Column containing protein description
      qDescCol: [description, REL]

      #
      # Frequencies
      #

      # Column containing scan frequency of pgm
      pdmFreq: [pgmFreq, REL]

      # Column containing scan frequency of p 
      pFreq: [pFreq, REL]

      # Column containing scan frequency of qc
      qcFreq: [qcFreq, REL]

      # Column containing protein frequency
      qFreq: [qFreq, REL]

      #
      # Mean differences
      #

      # Column containing positive values for "up" pgm and negative values for "down"
      # pgm in Z_pdm2p_dNM. For example, mean difference between Treatment Group and Control Group
      sign: [Z_pgm2p_dNM_dX, dX]

      # Column containing positive values for "up" pdm_NM and negative values for "down"
      # pdm in Z_pdm2p (without NM correction). For example, mean difference between Treatment Group and Control Group
      signNM: [Z_pgm2p_dX, dX]

      sign_p: [Z_p2qc_dX, dX]
          
      sign_qc: [Z_qc2q_dX, dX]
          

      #
      # qValue/pValue
      #

      # Column containing significance value (qvalue/pvalue) for Z_pgm2p_dNM (pgm corrected by their non modified version)
      qvalue_dNM: Z_pgm2p_dNM_limma

      # Column containing qvalue/pvalue for Z_pgm2p considering only non modified pgm
      qvalue_NM: Z_pgm2p_limma_NM_ONLY

      # Column containing qvalue of p
      qvalue_p: Z_p2qc_limma
          
      # Column containing qvalue of qc
      qvalue_qc: Z_qc2q_limma

      # Column containing number of missing cleavages
      missing_cleavages: [Missing_Cleavage, REL]


      #
      # PTM Frequency Table options
      #

      # Window size used to estimate probability
      x: 5

      # q-value column used to filter PTM (aa, dm)
      # Possible values: 'binom1-PSM', 'binom1-PDM', 'binom2-PSM', 'binom2-PDM'
      binom: binom1-PSM

      # Binomial q-value threshold applied to filter PTM
      q_thr: 0.01

      # Values represented in pivot table
      # Possible values: 'x-PSM', 'x-PDM'
      values_pivot: x-PSM

    ```
- **Optional inputs**
  * `--ptmmap` (`-p`): Path to PTMMap plots.
  * `--q2info` (`-q`): Path to report with protein information.

- **Outputs**
  * `--outdir` (`-o`): Directory where the protein reports will be saved:
    * qReports: Excel files summarizing PTM/NM statistics per protein.
    * Basal Reports: TSV files showing total PTM and NM counts across proteins, before significance filtering.
    * Frequency Tables: PTM-specific summary tables used to support statistical enrichment calculations.
    * Links to PTM maps directly embedded in Excel qReports.

[Go to top](#modules)
___




# References

[1] Rodríguez, J. M., Jorge, I., Martinez-Val, A., Barrero-Rodríguez, R., Magni, R., Núñez, E., Laguillo, A., Devesa, C. A., López, J. A., Camafeita, E., & Vázquez, J. (2023). iSanXoT: A standalone application for the integrative analysis of mass spectrometry-based quantitative proteomics data. Computational and structural biotechnology journal, 23, 452–459. https://doi.org/10.1016/j.csbj.2023.12.034

[2] Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47. https://doi.org/:10.1093/nar/gkv007