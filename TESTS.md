# Samples

Ensure Python is installed on your system along with the required dependencies.

For further information, consult the [INSTALLATION Guide](INSTALLATION.md).


## Samples 1: iSanXoT Report and PDM table for Mouse Heteroplasmia (Heart)

### Download the sample files:

+ On Linux:
```bash
cd samples && \
wget https://zenodo.org/records/15090841/files/heteroplasmic_heart.zip?download=1 -O heteroplasmic_heart.zip && \
unzip heteroplasmic_heart.zip && \
cd ..
```

or

+ On Windows:
```batch
@echo off
mkdir samples
cd samples
curl -L -o heteroplasmic_heart.zip https://zenodo.org/records/15090841/files/heteroplasmic_heart.zip?download=1 
powershell -Command "Expand-Archive -Path heteroplasmic_heart.zip -DestinationPath ."
cd ..
```

### Execute the programs for the current sample:

1. MergeiSanxotPDM:
```
python 1_MergeiSanxotPDM/MergeiSanxotPDM.py \
  -i samples/heteroplasmic_heart/inputs/q2all.tsv \
  -p samples/heteroplasmic_heart/inputs/experiment_PDMTable_GM_J.txt \
  -o samples/heteroplasmic_heart/results/isanxot_report_q2all_pdm.tsv
```

2. NMpyCompare:
```
python 2_NMpyCompare/NMpyCompare.py \
  -i samples/heteroplasmic_heart/results/isanxot_report_q2all_pdm.tsv \
  -c samples/heteroplasmic_heart/inputs/params.yml \
  -o samples/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM.tsv
```

3. ReportLimma:
```
Rscript --vanilla 3_ReportLimma_wo_GUI/app_wo_GUI.R \
  -i samples/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM.tsv \
  -c samples/heteroplasmic_heart/inputs/params.yml \
  -s samples/heteroplasmic_heart/inputs/limma_comparisons.tsv \
  -o samples/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM_LIMMA.tsv
```

4. FDRoptimizer:
```
python 4_FDRoptimizer/FDRoptimizer.py \
  -i samples/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM_LIMMA.tsv \
  -c samples/heteroplasmic_heart/inputs/params.yml \
  -o samples/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM_LIMMA_FDR.tsv
```

5. PTMMap:
```
python 5_PTMMap/PTMMap.py \
  -i samples/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM_LIMMA_FDR.tsv \
  -c samples/heteroplasmic_heart/inputs/params.yml \
  -o samples/heteroplasmic_heart/results
```

6. qTableReport:
```
python 6_qTableReport/qReportMaker.py \
  -i samples/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM_LIMMA_FDR.tsv \
  -q samples/heteroplasmic_heart/inputs/myMitocarta.tsv \
  -c samples/heteroplasmic_heart/inputs/params.yml \
  -o samples/heteroplasmic_heart/results
```

## Samples 2: iSanXoT Report and PDM table for Mouse Heteroplasmia (Liver)

You can download the input files for this sample from the following URL:

https://zenodo.org/records/15090841/files/heteroplasmic_liver.zip?download=1

To execute the pipeline, follow the same steps as in Sample 1.


## Samples 3: iSanXoT Report and PDM table for Mouse Heteroplasmia (Muscle)

You can download the input files for this sample from the following URL:

https://zenodo.org/records/15090841/files/heteroplasmic_muscle.zip?download=1

To execute the pipeline, follow the same steps as in Sample 1.