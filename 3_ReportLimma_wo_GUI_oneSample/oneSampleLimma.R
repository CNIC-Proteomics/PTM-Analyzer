library(limma)

'
https://bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
16.1 Swirl Zebrafish: A Single-Group Experiment
pg. 76
'

setwd("S:\\U_Proteomica\\LABS\\LAB_DF\\Comet_PTM\\FA_Biopsias\\2_ReportLimma_wo_GUI_oneSample")

groups <- c('B0', 'F', 'I')
integrations <- c('Z_pdm2pgm', 'Z_pgm2p', 'Z_pgm2p_dNM', 'Z_p2qf', 'Z_qf2q', 'Z_q2all')

intgr <- integrations[1]

for (g in groups) {
  

  for (intgr in integrations) {
  
    df <- read.csv(paste0(g, "\\", intgr, "\\", "input.tsv"), header=T, sep="\t")
    rownames(df) <- df$LEVEL
    df$LEVEL <- NULL
    
    samples <- colnames(df)
    design <- matrix(rep(1, length(samples)))
    rownames(design) <- samples
    colnames(design) <- g
    
    fit <- lmFit(df,design)
    fit <- eBayes(fit)
    
    pvalue <- as.data.frame(fit$p.value)
    meanRow <- as.data.frame(fit$Amean)
    
    output <- merge(meanRow, pvalue, 'row.names', all=T)
    colnames(output) <- c('LEVEL', 'Mean', 'pvalue')
    
    write.table(output, paste0(g, "\\", intgr, "\\", "output.tsv"), sep="\t", row.names=F, col.names=T)
  
  }
}