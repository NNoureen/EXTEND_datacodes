

###################### Co-expressed Genes Identification ###########################################################

GeneExp_Tumor_Pancan = read.table('Data.txt',sep='\t',head=T)

TERT_Group = GeneExp_Tumor_Pancan[which(GeneExp_Tumor_Pancan$Status == "TERT"),]
TERT_Group = TERT_Group[,-1]

## Find correlation between genes and pick up the upper triangle of the correlation matrix. 
res.cor <- TERT_Group %>%  # (1)
    correlate() %>%                  # (2)
  shave(upper = TRUE) %>%            # (3)
  stretch(na.rm = TRUE)         # (4)

TERTGroup = res.cor[which(res.cor$x == "TERT" | res.cor$y == "TERT"),]  
write.table(TERTGroup,'TERTGroup_LearntFromTERTSamples.txt',sep='\t',quote= FALSE,row.names = FALSE)  

############################ DEGS Identification #####################################################

Data = read.table('Data.txt',head=T, sep='\t',skip = 2)
row.names(Data)=Data$Name
Data = Data[,-c(1,2)]

TelMat=round(Data)  
  CPM = cpm(TelMat)
  cpm_filtered_matrix=TelMat
  DGE<-DGEList(cpm_filtered_matrix)
  DGE<-calcNormFactors(DGE,method =c("TMM")) 
  group <- as.factor(c(rep('TP',81),rep('TN',123)))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  cont.matrix <- makeContrasts(D1_comp= TP - TN ,levels=design)
  y <- estimateDisp(DGE, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  res <- glmQLFTest(fit, contrast=cont.matrix)
  restop = topTags(res, n=20000)
  write.table(restop,file= 'LGG_DEGS.txt',sep='\t',quote=FALSE,  row.names = TRUE)

#################################################################################


GeneSet1 = scan('TERC_CoexpressedGenes.txt','')
GeneSet2 = scan('TERT_Coexpressed.txt','')
Common = intersect(GeneSet1,GeneSet2)
writeLines(Common,'TERTCorr_DEGsCommon.txt',sep='\n')