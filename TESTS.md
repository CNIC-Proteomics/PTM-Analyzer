# TESTS

## Tests 1: iSanXoT report from mouse heteroplasmia (heart)

### 1. NMpyCompare
```
python 1_NMpyCompare/NMpyCompare.py \
  -i tests/heteroplasmic_heart/inputs/isanxot_report_q2all_pdm.tsv \
  -c tests/heteroplasmic_heart/inputs/params.yml \
  -o tests/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM.tsv
```

### 2. ReportLimma
```
Rscript --vanilla 2_ReportLimma_wo_GUI/app_wo_GUI.R \
  -i tests/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM.tsv \
  -c tests/heteroplasmic_heart/inputs/params.yml \
  -s tests/heteroplasmic_heart/inputs/limma_comparisons.tsv \
  -o tests/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM_LIMMA.tsv
```

### 3. FDRoptimizer
```
python 3_FDRoptimizer/FDRoptimizer.py \
  -i tests/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM_LIMMA.tsv \
  -c tests/heteroplasmic_heart/inputs/params.yml \
  -o tests/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM_LIMMA_FDR.tsv
```

### 4. PTMMap
```
python 4_PTMMap/PTMMap.py \
  -i tests/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM_LIMMA_FDR.tsv \
  -c tests/heteroplasmic_heart/inputs/params.yml \
  -o tests/heteroplasmic_heart/results/PTMMaps/H-C
```

### 5. qTableReport
```
python 5_qTableReport/qReportMaker.py \
  -i tests/heteroplasmic_heart/results/isanxot_report_q2all_pdm_NM_LIMMA_FDR.tsv \
  -q tests/heteroplasmic_heart/inputs/myMitocarta.tsv \
  -p tests/heteroplasmic_heart/results/PTMMaps \
  -c tests/heteroplasmic_heart/inputs/params.yml \
  -o tests/heteroplasmic_heart/results/qReports
```

