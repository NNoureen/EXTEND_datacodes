library('pracma')
library('ComplexHeatmap')
library('ggplot2')
library('gplots')
library('dplyr')
library('ggpubr')
stringsAsFactors=FALSE
library(circlize)
library(stringr)
library(lemon)
library(forcats)
library(corrr) 
library(gridBase)
library(gridGraphics)
library(ggpubr)
library(ggforce) ## geom_sina()
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")



##################################   FIGURE 1 ######################################################################################################################################################################


GSCData = read.table('GSC_Data_Fig1a.txt',sep='\t',head=T)   ### Data Figure 1a
TelomeraseScores_BCL_TERT = read.table('BCL_Data_Fig1c.txt',sep='\t',head=T)  ### Data Figure 1b
TelomeraseScores_LungAll = read.table('NSLC_Data_Fig1b.txt',sep='\t',head=T)   ### Data Figure 1c
Liposarcomas_WithNormals2 = read.table('Liposarcomas_Data_Fig1d.txt',sep='\t',head=T)  ### Data Figure 1d
NBM_Final = read.table('NBM_Data_Fig1e.txt',sep='\t',head=T)   ### Data Figure 1e

my_comparisons1 = list(c("ALT CellLine","Telomerase CellLine"),c("ALT CellLine","ALT Liposarcoma"),c("ALT Liposarcoma","Telomerase Liposarcoma"),c("Telomerase CellLine","Telomerase Liposarcoma"))


pdf('Figure1_Allparts.pdf')

########### Figure 1 Part a #############

P1 = ggscatter(GSCData, x = "RTA", y = "EXTENDScores",color = "red",size = 1.3,label = 'Samples',font.label=c(5,"black"),repel=TRUE,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)


########### Figure 1 Part b #############

P2 = ggscatter(TelomeraseScores_BCL_TERT, x = "TelomeraseAssay", y = "EXTENDScores",color = "red",size = 1.3,label = 'cellline',font.label=c(5,"black"),repel=TRUE,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE,conf.int.level = 0.65 # Add confidence interval
)


########### Figure 1 Part c #############

P3 = ggscatter(TelomeraseScores_LungAll, x = "Telomerase", y = "EXTENDScores",color = "red",size = 1.3,label = 'celltype',font.label=c(5,"black"),repel=TRUE,
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)


########### Figure 1 Part d and e #############

P4 = ggboxplot(Liposarcomas_WithNormals2,x="Mechanism",y="EXTENDScores",size =  0.2, width = 0.5,outlier.shape=NA,color = "Mechanism", add = "jitter",add.params= list(size=0.5),legend="none",order = c("ALT CellLine","Telomerase CellLine","ALT Liposarcoma","Telomerase Liposarcoma"))+ theme(axis.text.x=element_text(angle=30,size=8,vjust=1,hjust=1))+
		   stat_compare_means(comparisons = my_comparisons1, method= "t.test",size=1.5)

P5 = ggboxplot(NBM_Final, x = "Mechanism", y = "EXTENDScores",size =  0.2, width = 0.5,outlier.shape=NA,color="Mechanism",add = "jitter",add.params= list(size=0.2),order = c("NO_TMM","ALT","MYCN_Amp","TERT_High","TERT_Rearrangement"))+ theme(axis.text.x=element_text(angle=30,size=5,vjust=1,hjust=1),legend.position= "none")+guides(size = FALSE)


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,ylab="EXTEND Scores",xlab="RTA",title="a",font.title=8,subtitle="Rho = 0.48; P = 0.0096",font.subtitle=6),vp = define_region(row = 1, col = 1:2))
print(ggpar(P2,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,ylab="EXTEND Scores",xlab="log2_TEA_S6",title="b",font.title=8,subtitle="Rho = 0.72; P = 0.01",font.subtitle=6),vp = define_region(row = 1, col = 4:5))
print(ggpar(P3,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,ylab="EXTEND Scores",xlab="Telomerase(ddTRAP)",title="c",font.title=8,subtitle="Rho = 0.65; P = 0.008",font.subtitle=6),vp = define_region(row = 2, col = 1:2))
print(ggpar(P4,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,ylab="EXTEND Scores",xlab="Groups",title="d",font.title=8),vp = define_region(row = 2, col = 4:5))
print(ggpar(P5,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,ylab="EXTEND Scores",xlab="TMM Mechanism",title="e",font.title=8),vp = define_region(row = 3, col = 1:2))
dev.off()



##################################   FIGURE 2 ######################################################################################################################################################################



GTEX_All = read.table('GTEXData_Fig2a.txt',sep='\t',head=T)     #### Figure 2a DATA 



#### Figure 2 Part a ##############


pdf('GTEX_Fig2a.pdf',onefile=FALSE)
P = ggplot(GTEX_All, aes(x=fct_reorder(Tissues,EXTEND,.desc =FALSE)))
P = P+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,aes(y = TERT,fill="TERT"))
P = P+ geom_boxplot(size= 0.1,outlier.size = 0.0,outlier.shape=NA,aes(y = EXTEND*5.9,fill="EXTEND"))
P <- P + scale_y_continuous(sec.axis = sec_axis(~./5.9, name = "EXTEND Scores"))
P <- P + labs(y = "TERT",x = "GTEX Tissues",fill="Scores")+theme_classic()+theme(axis.text.x=element_text(angle=40,size=5,vjust=1,hjust=1),legend.key.size= unit(0.5,"cm"),legend.text=element_text(size=6),legend.title=element_text(size=8),legend.key.width = unit(0.5,"cm"))+guides(size = FALSE)
P <- P + scale_fill_manual(values = c("EXTEND" = "hotpink","TERT"="aquamarine3"))

grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 




print(ggpar(P,font.xtickslab =5,font.ytickslab =6,font.x = 7,font.y=7),vp = define_region(row = 1, col = 1:5))
dev.off()






#### Figure 2 Parts b and c ##############


Heart = read.table('HeartDevelop_Fig2b.txt',sep='\t',head=T)  #### Figure 2b DATA 
Liver = read.table('LiverDevelop_Fig2c.txt',sep='\t',head=T)  #### Figure 2c DATA 

pdf('Figure2b_c.pdf',onefile=FALSE)

P2 <- ggplot(Heart, aes(x = Time))
P2 <- P2 + geom_point(aes(y = TERT, colour = "TERT"))
P2= P2+geom_line(aes(x=Time,y=MeanTERT,group = 1,colour = "TERT"))
P2 <- P2 + geom_point(aes(y = EXTEND*0.42, colour = "EXTEND"))
P2 <- P2 + scale_y_continuous(sec.axis = sec_axis(~./0.42, name = "EXTEND Scores"))
P2 = P2+geom_line(aes(x=Time,y=MeanEXTEND*0.42,group = 1,colour = "EXTEND")) 
  # modifying colours and theme options
  P2 <- P2 + scale_colour_manual(values = c("EXTEND"="red","TERT"= "blue"))
  P2 <- P2 + labs(y = "TERT",x = "Heart Development TimePoints", colour = "Scores",size=6)+ theme_classic()+theme(axis.text.x=element_text(angle=50,vjust=1,hjust=1),legend.key.size= unit(0.01,"cm"),legend.text=element_text(size=6),legend.position="right",legend.title=element_text(size=7),legend.key.width = unit(0.05,"cm"))+guides(size = FALSE)

  
  
  
P3 <- ggplot(Liver, aes(x = Time))
P3 <- P3 + geom_point(aes(y = TERT, colour = "TERT"))
P3= P3+geom_line(aes(x=Time,y=MeanTERT,group = 1,colour = "TERT"))
P3 <- P3 + geom_point(aes(y = EXTEND*1.83, colour = "EXTEND"))
P3 <- P3 + scale_y_continuous(sec.axis = sec_axis(~./1.83, name = "EXTEND Scores"))
  # modifying colours and theme options
P3 = P3+geom_line(aes(x=Time,y=MeanEXTEND*1.83,group = 1,colour = "EXTEND")) 
P3 <- P3 + scale_colour_manual(values = c("EXTEND"="red","TERT"= "blue"))
P3 <- P3 + labs(y = "TERT",x = "Liver Development TimePoints", colour = "Scores",size=6)+ theme_classic()+theme(axis.text.x=element_text(angle=50,vjust=1,hjust=1),legend.key.size= unit(0.01,"cm"),legend.text=element_text(size=6),legend.position="right",legend.title=element_text(size=7),legend.key.width = unit(0.05,"cm"))+guides(size = FALSE)
  
  
  
  
  
grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 



print(ggpar(P2,font.xtickslab =5,font.ytickslab =6,font.x = 9,font.y=9),vp = define_region(row = 1, col = 1:4))
print(ggpar(P3,font.xtickslab =5,font.ytickslab =6,font.x = 9,font.y=9),vp = define_region(row = 2, col = 1:4))

dev.off()







########### Figure 2 Part d #############


GSE81507_Scores = read.table('GSE81507_Data_Fig1f.txt',sep='\t',head=T)  ### Data Figure 1f
GSE81507df2 <- data_summary(GSE81507_Scores, varname="Scores", groupnames=c("Samples", "Labels"))
GSE81507df1 = GSE81507df2[c(1,4,7,10),]
GSE81507df3 = GSE81507df2[c(2,5,8,11),]
GSE81507df4 = GSE81507df2[c(3,6,9,12),]


	
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

pdf('Figure2d.pdf')   
		   
P6 = ggplot(GSE81507df4, aes(x=Samples,y=Scores,fill=Labels))
P6 = P6+ geom_bar(stat="identity", color="black",position=position_dodge(),width=0.6,fill="pink")+geom_errorbar(aes(ymin=Scores, ymax=Scores+sd), width=.2,
                 position=position_dodge(.9))
P6 <- P6 + labs(y = "EXTEND Scores",x = "Samples",fill="Labels")+theme_classic()+theme(axis.text.x=element_text(angle=30,size=8,vjust=1,hjust=1),legend.key.size= unit(0.5,"cm"),legend.key.width = unit(0.3,"cm"),legend.position = "none")

grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 



print(ggpar(P6,font.xtickslab =5,font.ytickslab =6,font.x = 9,font.y=9),vp = define_region(row = 1, col = 1:4))

dev.off()




##################################   FIGURE 3 ######################################################################################################################################################################


#### Figure 3 Part a ##############

Pancan_TumorNormal = read.table('TCGA_Data_Figure3a.txt',sep='\t',head=T)    #### Figure 3a DATA #######

pdf('Figure3a.pdf',onefile=FALSE)

P = ggviolin(Pancan_TumorNormal,x="Cancer",y="EXTENDScores",size =  0.1, width = 1.5,draw_quantiles=c(0.5),outlier.shape=NA,color = "black",fill="Groups",palette =c("lightblue","pink"),
order = c('PAAD','KIRC','THCA','PRAD','KIRP','LUAD','KICH','LIHC','BRCA','HNSC','LUSC','STAD','UCEC','ESCA','CRC','BLCA'))+theme_classic()+theme(axis.text.x=element_text(angle=50,size=7,vjust=0.5),legend.position = "top",legend.key.size= unit(0.3,"cm"),legend.key.width = unit(0.3,"cm"),legend.title = element_text(size=7),legend.text =element_text(size=6))


grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow =7, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P,font.xtickslab =c("black",6),font.ytickslab =c("black",6),font.x = 7,font.y=7,font.legend=6,ylab="EXTEND Scores",xlab="Cancer"),vp = define_region(row = 1:2, col = 1:3))


dev.off()


#### Figure 3 Part b ##############



AllTumors = read.table('TCGA_DiseaseStage_Fig3bData.txt',sep='\t',head=T)
Cancers = c('CRC','HNSC','KIRC','KIRP','LUAD','STAD','THCA')
AllTumors = AllTumors[which(AllTumors$Cancer %in% Cancers),]

pdf('Figure3b.pdf',onefile=FALSE)

P2 = ggboxplot(AllTumors,x="Cancer",y="EXTEND",size =  0.2, width = 0.5,outlier.shape=NA,color="black",fill = "Level",palette =c("cyan","green","pink","grey"),
order = c('THCA','KIRC','KIRP','LUAD','HNSC','CRC','STAD'))+theme_classic()+theme(axis.text.x=element_text(angle=50,size=7,vjust=0.5),legend.position = "top",legend.key.size= unit(0.3,"cm"),legend.key.width = unit(0.3,"cm"),legend.title = element_text(size=7),legend.text =element_text(size=6))



grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow =9, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P2,font.xtickslab =c("black",6),font.ytickslab =c("black",6),font.x = 7,font.y=7,font.legend=6,xlab="AllTumors"),vp = define_region(row = 4:6, col = 1:2))


dev.off()





#### Figure 3 Part c ##############

AllTumors = read.table('SKCM_Stage_Levels_Figure3cData.txt',sep='\t',head=T)


pdf('Figure3c.pdf',onefile=FALSE)

P4 = ggviolin(AllTumors,x="TumorTissueSite",y="EXTEND",size =  0.1, order = c('primary tumor','RegionalCut_SubCut','regional lymph node','distant metastasis'),width = 0.9,draw_quantiles=c(0.5),outlier.shape=NA,color = "black",fill="maroon",alpha=0.5)+theme_classic()+theme(axis.text.x=element_text(angle=30,size=7,vjust=1,hjust=1),legend.position = "top",legend.key.size= unit(0.3,"cm"),legend.key.width = unit(0.3,"cm"),legend.title = element_text(size=7),legend.text =element_text(size=6))+stat_compare_means(comparisons = my_comparisons5, method= "t.test",size=1.9)

grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow =2, ncol =2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P4,font.xtickslab =c("black",6),font.ytickslab =c("black",6),font.x = 7,font.y=7,font.legend=6,xlab="SKCM TissueSite",ylab="EXTEND Scores"),vp = define_region(row = 1, col = 1))

dev.off()



#### Figure 3 Part d ##############

HR_Univariate_OS = read.table('HR_OS_FinalData_Fig3d.txt',sep='\t',head=T)



Cancer <- HR_Univariate_OS$Cancer
HazardRatio  <- HR_Univariate_OS$Exp_Coef
LowerHR <- HR_Univariate_OS$Lower_95
UpperHR <- HR_Univariate_OS$Upper_95
Pval = HR_Univariate_OS$Pvalue
FDR = p.adjust(Pval,method="bonferroni")

dfOS <- data.frame(Cancer, HazardRatio, LowerHR, UpperHR,Pval,FDR)

dfOS$HazardRatio = log(dfOS$HazardRatio)
dfOS$LowerHR = log(dfOS$LowerHR)
dfOS$UpperHR = log(dfOS$UpperHR)


dfOS$HazardRatio = round(dfOS$HazardRatio,2)
dfOS$LowerHR = round(dfOS$LowerHR,2)
dfOS$UpperHR = round(dfOS$UpperHR,2)
dfOS$Pval = round(dfOS$Pval,2)
dfOS$FDR = round(dfOS$FDR,2)


# reverses the factor level ordering for labels after coord_flip()
dfOS$Cancer <- factor(dfOS$Cancer, levels=rev(dfOS$Cancer))


dfOS$Sig = NA

for(i in 1:nrow(dfOS)){

if(dfOS$Pval[i]< 0.055){
dfOS$Sig[i] = "*"}

}

pdf('Fig3d.pdf')


P6 = ggplot(data=dfOS, aes(x=fct_reorder(Cancer,HazardRatio,.desc =FALSE), y=HazardRatio, ymin=LowerHR, ymax=UpperHR)) +
        geom_point(size=1, color="black", fill="pink") + 
		geom_errorbar(width=0.2,size=0.1)+
        geom_hline(yintercept=0, lty=2,color="blue",size=0.1) +  # add a dotted line at x=1 after flip
		xlab("Cancer") + ylab("Hazard Ratio OS(95% CI)") +
        theme_classic()+theme(axis.text.x=element_text(angle= 50,size=8,vjust=1,hjust=1))  # use a white backgroundcoord_flip() +  # flip coordinates (puts labels on y axis)
		
P6 = P6+geom_text(aes(label=Sig), hjust = 0,vjust=0.5,nudge_y = 2,size=4)	






grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow =3, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P6,font.xtickslab =c("black",6),font.ytickslab =c("black",5),font.x = 7,font.y=7,font.legend=5),vp = define_region(row = 2, col = 1:2))
dev.off()



#### Figure 3 Part e ##############


FDRPancan =read.table('SignalingPathways_Fig3E_data.txt',sep='\t',head=T)

pdf('Fig_3E_Pathways.pdf')
   
P3 = ggplot(FDRPancan, aes(x = Cancer, y = GeneOrder2,color=Corr))+geom_point(size=FDRPancan$Significance)+scale_color_gradient2(midpoint=0, low="navyblue", mid="white",high="maroon")+theme_classic()+theme(axis.text.x=element_text(angle= 50,size=8,vjust=1,hjust=1),legend.key.height=unit(0.1,"cm"),legend.position="bottom")



grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 10, ncol = 10)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P3,font.xtickslab =4,font.ytickslab =5,font.x = 7,font.y=7, font.legend=5),vp = define_region(row = 1:7, col = 6:10))

dev.off()




Barplot = read.table('GeneOrder_Barplot_Figure3E_data.txt',sep='\t',head=T)



pdf('Fig_3E_BarplotGenes.pdf')
P1 = ggplot(Barplot, aes(x=GeneOrder, y=NoofCancerTypes,fill=Status)) + 
  geom_bar(stat = "identity",width=0.7)+scale_fill_manual(values=c("UP"="red","DN"="blue"))+theme_classic()+theme(legend.key.size=unit(0.1,"cm"))+coord_flip()

#P2 = ggplot(Dngenes, aes(y=logFDR, x=Pathway)) + 
#  geom_bar(stat = "identity",fill="blue",width=0.7)+theme_classic()+theme(legend.key.size= unit(0.05,"cm"))+coord_flip()

  
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 6, ncol = 4)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(P1,font.xtickslab =6,font.ytickslab =4,font.x = 7,font.y=7,font.title=6,font.legend=6),vp = define_region(row = 1:3, col = 1))  
#print(ggpar(P2,font.xtickslab =6,font.ytickslab =4,font.x = 7,font.y=7,font.title=6,font.legend=6,),vp = define_region(row = 1, col = 2))  

dev.off()




######################################################################################### Figure 4 #######################################################################################

#### Figure 4 Part a ##############


Pancan_Stemness_All_WithFDR = read.table('Stemness_Fig4aData.txt',sep='\t',head=T)

Pancan_Stemness_All_WithFDR$Corr = as.factor(Pancan_Stemness_All_WithFDR$Corr)

pdf('Fig4a.pdf')
P1 = ggscatter(Pancan_Stemness_All_WithFDR, x="Stemness", y="EXTEND",size ="Corr",color = "red",label="Cancer",font.label=c(5,"black"),repel=TRUE,alpha=0.3,add = "reg.line", 
add.params = list(color = "blue", fill = "grey",size=0.2), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)


grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow =7, ncol = 8)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 



print(ggpar(P1,font.xtickslab =c("black",7),font.ytickslab =c("black",7),font.x = 9,font.y=9,font.legend=7,ylab= "EXTEND Scores"),vp = define_region(row = 1:3, col = 4:6))

dev.off()




#### Figure 4 Part b ##############


EXTEND_stemness_ssGBM = read.table('EXTEND_stemness_ssGBM_Data4b_4e.txt',sep='\t',head=T)

  
pdf('Fig4b.pdf')
  
x1 <- EXTEND_stemness_ssGBM$Stemness
x2 <- EXTEND_stemness_ssGBM$EXTENDScores
df <- data.frame(x1,x2)

## Use densCols() output to get density at each point
x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]
data=df[order(df$dens),]


sp2 = ggscatter(data, x="x1", y="x2",size =0.6,color = "dens", add = "reg.line", 
add.params = list(color = "blue", fill = "grey",size=0.2), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ scale_color_gradientn(colors = rainbow(10))+ stat_cor(method = "spearman",size=2)  +theme_classic()+theme(legend.key.size= unit(0.4,"cm"))
 


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(sp2,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,font.title=6,font.legend=6,xlab="Stemness",ylab="EXTEND Scores"),vp = define_region(row = 1, col = 1))  

dev.off()






#### Figure 4 Part c ##############




TM_Stemness = read.table('scHNSC_Fig4c_Data.txt',sep='\t',head=T)

x1 <- TM_Stemness$Stemness
x2 <- TM_Stemness$NormEXTENDScores
df <- data.frame(x1,x2)

## Use densCols() output to get density at each point
x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]
data=df[order(df$dens),]

 
pdf('Fig4c.pdf')

sp2 = ggscatter(data, x="x1", y="x2",size =0.6,color = "dens", add = "reg.line", 
add.params = list(color = "blue", fill = "grey",size=0.2), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ scale_color_gradientn(colors = rainbow(10))+ stat_cor(method = "spearman",size=2)  +theme_classic()+theme(legend.key.size= unit(0.4,"cm"))
   
 grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 3 ,ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(sp2,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,font.title=6,font.legend=6,xlab="Stemness",ylab="EXTEND Scores"),vp = define_region(row = 1, col = 1))  

dev.off()



#### Figure 4 Part e ##############

EXTEND_stemness_ssGBM = read.table('EXTEND_stemness_ssGBM_Data4b_4e.txt',sep='\t',head=T)


my_comparisons1 = list(c("G1_S","G2_M"),c("G1_S","NonCycling"),c("G2_M","NonCycling"))


pdf('Fig4e.pdf')

P1 = ggplot(data = EXTEND_stemness_ssGBM, 
         aes(x = Phase, y = EXTENDScores, fill = "aquamarine3")) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0),fill="aquamarine3") +
  geom_point(aes(y = EXTENDScores),color="green4", 
             position = position_jitter(width = .15), size = .2) +
  geom_boxplot(width = .1, outlier.shape = NA,fill="white") +
  #expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.position="none")+stat_compare_means(comparisons = my_comparisons1, method= "t.test",size=2)# coord_flip() + # flip or not
  

grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 6, ncol = 6)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7),vp = define_region(row = 1:2, col = 1:2))
 
dev.off()




#### Figure 4 Part f ##############


TM_Stemness = read.table('scHNSC_Fig4c_Data.txt',sep='\t',head=T)
my_comparisons1 = list(c("G1_S","G2_M"),c("G1_S","NonCycling"),c("G2_M","NonCycling"))


pdf('Fig4f.pdf')
  
P1 = ggplot(data = TM_Stemness, 
         aes(x = Phase, y = NormEXTENDScores, fill = "aquamarine3")) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0),fill="aquamarine3") +
  geom_point(aes(y = NormEXTENDScores),color="green4", 
             position = position_jitter(width = .15), size = .2) +
  geom_boxplot(width = .1, outlier.shape = NA,fill="white") +
  #expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.position="none")+stat_compare_means(comparisons = my_comparisons1, method= "t.test",size=2)# coord_flip() + # flip or not

  
  
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 6, ncol = 6)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,ylab="EXTEND Scores"),vp = define_region(row = 1:2, col = 1:2))
dev.off()



  


#### Figure 4 Part g ##############



Stemness_FDR = read.table('Prolif_Fig4gData.txt',sep='\t',head=T)
Stemness_FDR$Corr = as.factor(Stemness_FDR$Corr)


pdf('Fig4g.pdf')
P1 = ggscatter(Stemness_FDR, x="MKI67", y="EXTEND",size ="Corr",color = "SignificanceMKI67",palette =c("Significant"="red","NotSign"="mediumblue"),label="Cancer",font.label=c(5,"black"),repel=TRUE,alpha=0.3,add = "reg.line", 
add.params = list(color = "blue", fill = "grey",size=0.2), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "spearman",size=2)



grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow =7, ncol = 8)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 



print(ggpar(P1,font.xtickslab =c("black",7),font.ytickslab =c("black",7),font.x = 9,font.y=9,font.legend=7,ylab="EXTEND Scores"),vp = define_region(row = 1:3, col = 4:6))

dev.off()




