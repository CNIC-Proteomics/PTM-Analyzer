# Samples

Ensure Python is installed on your system along with the required dependencies.

For further information, consult the [INSTALLATION Guide](INSTALLATION.md).



## Samples 1: iSanXoT Report and PDM table for Mouse Heteroplasmia (Heart)

Input files for PTM-Analyzer, derived from the quantification results of iSanXoT [1] and the PDM table from the PTM-compass [2].

The samples originate from `heart` tissue, based on the study by Bagwan N, Bonzon-Kulichenko E, Calvo E, et al. [3].

### Download the sample files:

+ On Linux:
```bash
cd samples && \
wget https://zenodo.org/records/17141424/files/heteroplasmic_heart.zip?download=1 -O heteroplasmic_heart.zip && \
unzip heteroplasmic_heart.zip && \
cd ..
```

or

+ On Windows:
```batch
@echo off
mkdir samples
cd samples
curl -L -o heteroplasmic_heart.zip https://zenodo.org/records/17141424/files/heteroplasmic_heart.zip?download=1 
powershell -Command "Expand-Archive -Path heteroplasmic_heart.zip -DestinationPath ."
cd ..
@echo on
```

### Execute the programs for the current sample:

0. Prepare workspace:
```bash
mkdir samples/heteroplasmic_heart/results
```

1. mergeiSanxotPGM:
```
python src/mergeiSanxotPGM.py \
  -i samples/heteroplasmic_heart/inputs/q_all.tsv \
  -p samples/heteroplasmic_heart/inputs/DMTable_PeakAssignation_FDRfiltered_DM0S_PA_T_PeakAssignation_SS_Heart_FDR_PDMTable_GM_J_PGM_Table_pgmFreq.tsv \
  -o samples/heteroplasmic_heart/results/quant_pgm.tsv
```
Note: This program is optional. If you already have the PDM table at the PGM level in the iSanXoT report (q_all.tsv file), you do not need to run it.

2. NMpyCompare:
```
python src/NMpyCompare.py \
  -i samples/heteroplasmic_heart/results/quant_pgm.tsv \
  -c samples/heteroplasmic_heart/inputs/params.yml \
  -o samples/heteroplasmic_heart/results/quant_pgm_NM.tsv
```

3. ReportLimma:
You can run the R script directly using your own Rscript installation:
```
Rscript --vanilla src/ReportLimma.R \
  -i samples/heteroplasmic_heart/results/quant_pgm_NM.tsv \
  -c samples/heteroplasmic_heart/inputs/params.yml \
  -s samples/heteroplasmic_heart/inputs/limma_comparisons.tsv \
  -o samples/heteroplasmic_heart/results/quant_pgm_NM_LIMMA.tsv
```

Alternatively, you can use the R-portable program and execute the "ReportLimma" batch script on Windows:

3.1 Open a command prompt window.

3.2 Execute the `ReportLimma` batch script using a R-portable.
```
bin\ReportLimma.bat ^
  -i "%CD%/samples/heteroplasmic_heart/results/quant_pgm_NM.tsv" ^
  -c "%CD%/samples/heteroplasmic_heart/inputs/params.yml" ^
  -s "%CD%/samples/heteroplasmic_heart/inputs/limma_comparisons.tsv" ^
  -o "%CD%/samples/heteroplasmic_heart/results/quant_pgm_NM_LIMMA.tsv"
```

4. FDRoptimizer:
```
python src/FDRoptimizer.py \
  -i samples/heteroplasmic_heart/results/quant_pgm_NM_LIMMA.tsv \
  -c samples/heteroplasmic_heart/inputs/params.yml \
  -o samples/heteroplasmic_heart/results/quant_pgm_NM_LIMMA_FDR.tsv
```

5. PTMMap:
```
python src/PTMMap.py \
  -i samples/heteroplasmic_heart/results/quant_pgm_NM_LIMMA_FDR.tsv \
  -c samples/heteroplasmic_heart/inputs/params.yml \
  -o samples/heteroplasmic_heart/results
```

6. qTableReport:
```
python src/qReportMaker.py \
  -i samples/heteroplasmic_heart/results/quant_pgm_NM_LIMMA_FDR.tsv \
  -q samples/heteroplasmic_heart/inputs/myMitocarta.tsv \
  -c samples/heteroplasmic_heart/inputs/params.yml \
  -o samples/heteroplasmic_heart/results
```


## Samples 2: iSanXoT Report and PDM table for Mouse Heteroplasmia (Liver)

You can download the input files for this `liver` sample, derived from the study by Bagwan N, Bonzon-Kulichenko E, Calvo E, et al. [1] at the following URL:

https://zenodo.org/records/17141424/files/heteroplasmic_liver.zip?download=1

To execute the pipeline, follow the same steps as in Sample 1.


## Samples 3: iSanXoT Report and PDM table for Mouse Heteroplasmia (Muscle)

You can download the input files for this sample from the following URL:

https://zenodo.org/records/17141424/files/heteroplasmic_muscle.zip?download=1

To execute the pipeline, follow the same steps as in Sample 1.




# References

[1] Rodríguez, Jose Manuel et al. *iSanXoT: A standalone application for the integrative analysis of mass spectrometry-based quantitative proteomics data.* Computational and structural biotechnology journal vol. 23 452-459. 26 Dec. 2023, doi: https://doi.org/10.1016/j.csbj.2023.12.034

[2] Cristina A. Devesa, Rafael Barrero-Rodríguez, Andrea Laguillo-Gómez, et al. *Integrative multi-layer workflow for quantitative analysis of post-translational modifications.* bioRxiv 2025.01.20.633864; doi: https://doi.org/10.1101/2025.01.20.633864

[3] Bagwan N, Bonzon-Kulichenko E, Calvo E, et al. *Comprehensive Quantification of the Modified Proteome Reveals Oxidative Heart Damage in Mitochondrial Heteroplasmy*. Cell Reports. 2018;23(12):3685-3697.e4. https://doi.org/10.1016/j.celrep.2018.05.080