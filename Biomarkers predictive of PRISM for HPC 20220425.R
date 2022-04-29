#which genes are associated with sensitivity to each drug mechanism-of-action?
#BG 20211207; last edit: BG 20220428

rm(list=ls(all=TRUE))

path.base <- paste0(getwd(),"/")
path.inputs <- paste0(path.base,"Inputs/")
path.outputs.RNAseq <- paste0(path.base,"RNAseq/")
path.outputs.proteomic <- paste0(path.base,"proteomic/")

dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE) #creates personal library
.libPaths(Sys.getenv("R_LIBS_USER")) #add to the path

if (!require(devtools)){install.packages(dev.tools)}
devtools::install_github('BelindaBGarana/DMEA')
if(require(rlang)){remove.packages("rlang")}
if(require(plyr)){remove.packages("plyr")}
if(require(dplyr)){remove.packages("dplyr")}
install.packages(c("rlang","GSA","plyr","dplyr","data.table","ggplot2","gridExtra","sjmisc","parallel","snow","doSNOW","viridis","tibble","stringr"), repos = "http://cran.us.r-project.org");
library(DMEA);
library(GSA);library(rlang);library(plyr);library(dplyr);library(data.table);library(ggplot2);

##Step 0: Import datasets and reduce for common gene names and adherent cancer cell lines
setwd(path.inputs)

#import PRISM data, info, and gmt
PRISM.AUC <- read.csv(file="PRISM_drug_mean_AUC_6-23-21.csv",header=T) #481 cell lines
PRISM.AUC$X <- NULL

#Import drug info
drug.info <- read.csv(file = "PRISM_secondary-screen-replicate-treatment-info.csv", header=T)
colnames(drug.info)[colnames(drug.info)=="name"] <- "Drug"
drug.info$Drug <- gsub("[[:punct:][:blank:]]",".",drug.info$Drug)
drug.moa <- na.omit(distinct(drug.info[,c("Drug","moa")]))

#Import gmt
gmt <- GSA.read.gmt(file="MOA_gmt_file_n6 wo special chars.gmt")
moa <- gmt$geneset.names

#import RNAseq data
if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos = "http://cran.us.r-project.org")}

# The following initializes usage of Bioc devel
BiocManager::install(version='3.14') #using '3.14' instead of 'devel' because USC HPC doesn't have R v4.2
BiocManager::install("depmap")
library("plyr");library("dplyr");library("ggplot2");library("viridis");
library("tibble");library("gridExtra");library("stringr");library("depmap");library("ExperimentHub");
eh <- ExperimentHub()
query(eh, "depamp")

RNA.data <- depmap_TPM()
RNA.data <- RNA.data[RNA.data$cell_line %in% overlap, ]
RNA.df <- as.data.frame(RNA.data)
RNA.df <- reshape2::dcast(RNA.df, cell_line ~ gene_name, value.var="rna_expression")

prot.data <- depmap_proteomic()
prot.df <- as.data.frame(prot.data)
prot.df <- ddply(prot.df, .(cell_line, gene_name), summarize, 
                 avg_expr = mean(protein_expression, na.rm=TRUE),
                 N_unique_id = length(unique(protein_id)))
prot.df <- reshape2::dcast(prot.df, cell_line ~ gene_name, value.var="avg_expr") #12,197 gene names
#prot.df.3 <- reshape2::dcast(prot.df[prot.df$N_unique_id>3,], cell_line ~ gene_name, value.var="avg_expr") #just 11 gene names, one of which is NA

#get info for adherent cancer cell lines
info <- read.csv(file="CCLE_sample_info.csv",header=T)
prism.info <- info[info$CCLE_Name %in% PRISM.AUC$CCLE_ID,] #480
info.adherent <- prism.info[prism.info$culture_type=="Adherent",] #320
overlap.CCLE <- info.adherent$CCLE_Name[info.adherent$CCLE_Name %in% RNA.data$CCLE_ID] #320 overlapping adherent cancer cell lines

#import DMEA results across all moa
all.DMEA <- read.csv(file="PRISM_DMEA_AUC_per_cell_line.csv",header=T)
all.DMEA$X <- NULL

##run correlations for each moa: DMEA NES vs. RNAseq
#split up samples
samples.grouped <- read.csv(file="Adherent_cancer_cell_lines_PRISM_CCLE_RNAseq_5groups.csv",header=T)

group1 <- samples.grouped[samples.grouped$group==1,]$CCLE_ID
group2 <- samples.grouped[samples.grouped$group==2,]$CCLE_ID
group3 <- samples.grouped[samples.grouped$group==3,]$CCLE_ID
group4 <- samples.grouped[samples.grouped$group==4,]$CCLE_ID
group5 <- samples.grouped[samples.grouped$group==5,]$CCLE_ID

#create training sets
train1 <- unique(c(group2, group3, group4, group5)) #leaves out group1
train2 <- unique(c(group1, group3, group4, group5)) #leaves out group2
train3 <- unique(c(group1, group2, group4, group5)) #leaves out group3
train4 <- unique(c(group1, group2, group3, group5)) #leaves out group4
train5 <- unique(c(group1, group2, group3, group4)) #leaves out group5

#load function for correlations
DMEA.biomarkers <- function(overlap.CCLE, biomarker.data, biomarker.type="RNAseq", moa, all.DMEA, training.samples, exp){
  print(paste0("Running correlations for experiment ",exp," using ",biomarker.type," data"))
  biomarker.corr.list <- list()
  corr.input.list <- list()
  for(i in 1:length(moa)){
    moa.data <- all.DMEA[all.DMEA$Drug.Set==moa[i],c("CCLE_ID","KS_Normalized")]
    train.moa.data <- moa.data[moa.data$CCLE_ID %in% training.samples,]
    corr.input.list[[moa[i]]] <- train.moa.data
    if(nrow(moa.data)>2){
      #merge NES values for 1 moa with biomarker data for all genes
      moa.biomarker <- merge(train.moa.data,biomarker.data,by="CCLE_ID")
      
      corr.biomarker <- rank.corr(moa.biomarker,value=biomarker.type,plots=FALSE)
      biomarker.corr.list[[moa[i]]] <- corr.biomarker$result
      #if(!is.na(corr.biomarker$scatter.plots[1])){ggsave(filename=paste0("Exp_",exp,"_",moa[i],"_Spearman_correlations_",biomarker.type,"_vs_drug_AUC_NES_adherent_FDR0.05.pdf"),corr.biomarker$scatter.plots)}
    }
  }
  all.biomarker.corr <- rbindlist(biomarker.corr.list, use.names=TRUE, idcol="moa")
  write.csv(all.biomarker.corr,file=paste0("Exp_",exp,"_all_adherent_correlations_",biomarker.type,"_vs_AUC_NES_w_moa.csv"))
}

#load function for DMEA
DMEA.biomarker.validation <- function(all.biomarker.corr, biomarker.type, n.genes=10, type.genes="min", moa, expr.data, expr.type, sample="CCLE_ID", overlap.CCLE,type="pearson",training.samples,exp){
  top.list <- list()
  WV.list <- list()
  corr.list <- list()
  all.DMEA.list <- list()
  subset.self.DMEA <- as.data.frame(moa)
  for(i in 1:length(moa)){
    biomarker.corr <- all.biomarker.corr[all.biomarker.corr$moa==moa[i],]
    if(nrow(biomarker.corr)>0){
      #remove any overlapping cell lines
      expr.data.nooverlap <- expr.data[!expr.data$CCLE_ID %in% training.samples, ]
      rownames(expr.data.nooverlap) <- expr.data.nooverlap$CCLE_ID
      
      #extract top correlated genes for DMEA
      if(type.genes=="min"){
        if(type=="pearson"){
          top <- biomarker.corr %>% slice_min(Pearson.est,n=n.genes)
          top$abs.Pearson.est <- abs(as.numeric(top$Pearson.est))
          weight.val <- "abs.Pearson.est"
          est <- "Pearson.est"
        }else if(type=="spearman"){
          top <- biomarker.corr %>% slice_min(Spearman.est,n=n.genes)
          top$abs.Spearman.est <- abs(as.numeric(top$Spearman.est))
          weight.val <- "abs.Spearman.est"
          est <- "Spearman.est"
        }else{
          print("Parameter 'type' must be either 'pearson' (default) or 'spearman'.")
          break
        }
      }else if(type.genes=="max"){
        if(type=="pearson"){
          top <- biomarker.corr %>% slice_max(Pearson.est,n=n.genes)
          weight.val <- "Pearson.est"
          est <- "Pearson.est"
        }else if(type=="spearman"){
          top <- biomarker.corr %>% slice_max(Spearman.est,n=n.genes)
          weight.val <- "Spearman.est"
          est <- "Spearman.est"
        }else{
          print("Parameter 'type' must be either 'pearson' (default) or 'spearman'.")
          break
        }
      }else if(type.genes=="min&max"){
        if(type=="pearson"){
          biomarker.corr$abs.Pearson.est <- abs(as.numeric(biomarker.corr$Pearson.est))
          top <- biomarker.corr %>% slice_max(abs.Pearson.est,n=n.genes)
          weight.val <- "abs.Pearson.est"
          est <- "Pearson.est"
        }else if(type=="spearman"){
          biomarker.corr$abs.Spearman.est <- abs(as.numeric(biomarker.corr$Spearman.est))
          top <- biomarker.corr %>% slice_max(abs.Spearman.est,n=n.genes)
          weight.val <- "abs.Spearman.est"
          est <- "Spearman.est"
        }else{
          print("Parameter 'type' must be either 'pearson' (default) or 'spearman'.")
          break
        }
      }else{
        print("Parameter 'type.genes' must be either 'min' (default), 'max', or 'min&max'.")
        break
      }
      top.list[[moa[i]]] <- top
      
      #remove NAs from relevant expr.data
      expr.data.nooverlap <- expr.data.nooverlap[expr.data.nooverlap$CCLE_ID %in% overlap.CCLE,]
      filtered.expr.data.nooverlap <- expr.data.nooverlap %>% select(c(sample,all_of(top$Gene)))
      filtered.expr.data.nooverlap <- na.omit(filtered.expr.data.nooverlap)
      
      if(nrow(filtered.expr.data.nooverlap)>2){
        #Run DMEA
        DMEA.test <- DMEA(drug.sensitivity = PRISM.AUC, gmt=gmt, expression = filtered.expr.data.nooverlap, weights=top, gene.names="Gene",weight.values = weight.val,estimate=est,scatter.plots=FALSE,scatter.plot.type = type)
        WV.list[[moa[i]]] <- DMEA.test$WV.scores
        corr.list[[moa[i]]] <- merge(DMEA.test$corr.result,drug.moa,by="Drug")
        #if(!is.na(DMEA.test$corr.scatter.plots[1])){ggsave(filename=paste0("Exp_",exp,"_",moa[i],"_",type,"_correlations_top_",n.genes,"_genes_",biomarker.type,"_",expr.type,"_corr_WGV_v_drug_AUC_adherent_FDR0.05.pdf"),DMEA.test$corr.scatter.plots)}
        DMEA.results <- DMEA.test$DMEA.result
        all.DMEA.list[[moa[i]]] <- DMEA.results
        self.DMEA <- DMEA.results[DMEA.results$Drug.Set==moa[i],] 
        subset.self.DMEA[subset.self.DMEA$moa==moa[i],2:6] <- self.DMEA[,3:7]
        if(self.DMEA$p_value<0.05 & self.DMEA$FDR_q_value<0.25){
          for(j in 1:length(DMEA.test$DMEA.mtn.plots)){
            ggsave(file=paste0("Exp_",exp,"_",moa[i],"_DMEA_top_",n.genes,"_genes_",biomarker.type,"_",expr.type,"_adherent_",type,"_mtnplot_",j,".pdf"),DMEA.test$DMEA.mtn.plots[[j]])
          }
          
          #Volcano plot
          #Categorize data by significance level
          DMEA.results$Significance <- "FDR > 0.25"
          DMEA.results[DMEA.results$FDR_q_value < 0.25,]$Significance <- "FDR < 0.25"
          DMEA.results$Significance <- factor(DMEA.results$Significance,levels=c("FDR < 0.25","FDR > 0.25"))
          if(nrow(DMEA.results[DMEA.results$p_value==0,])>0){DMEA.results[DMEA.results$p_value==0,]$p_value <- 0.00099}
          #Define x, y limits
          limit.x <- 3
          limit.y <- 3.5
          #Plot
          volc_DMEA <- ggplot(data = DMEA.results, aes(x = KS_Normalized, y = -log(p_value,10), color=Significance)) +
            geom_point(size = 6) + ggrepel::geom_text_repel(data=subset(DMEA.results,Significance=="FDR < 0.25"),mapping=aes(label=Drug.Set,size=I(7.5)),nudge_y = 0.25) +
            scale_color_manual(values=c("red","azure4"),name="Significance",breaks=c("FDR < 0.25","FDR > 0.25")) +
            xlim(-limit.x,limit.x) + ylim(0,limit.y) + xlab("Normalized Enrichment Score") + ylab("-Log(p-value)") +
            geom_vline(xintercept=0,linetype="solid",color="grey",size=1) +
            #geom_hline(yintercept=-log(0.25,10), linetype="dashed",color="red",size=1) +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = 'black', size = 0.65),
                  legend.text=element_text(size=18),axis.text=element_text(size=22),axis.title=element_text(size=25,face="bold"),
                  panel.background = element_rect(fill="white", colour="white", size=0.5,linetype="solid", color="black"), text = element_text(size = 20),
                  legend.position = c(0.5,0.75), legend.background=element_rect(fill="white",color="black"))
          ggsave(file=paste0("Exp_",exp,"_",moa[i],"_volcano_top_",n.genes,"_genes_",biomarker.type,"_",expr.type,"_NES_vs_minus_logp_",type,".pdf"), volc_DMEA, width=11,height=8.5)
        }
      }
    }
  }
  all.top <- rbindlist(top.list, use.names=TRUE, idcol="moa")
  write.csv(all.top,file=paste0("Exp_",exp,"_","All_adherent_top_",n.genes,"_weights_",biomarker.type,"_",expr.type,"_corr_AUC_NES_",type,".csv"))
  all.WV <- rbindlist(WV.list, use.names=TRUE, idcol="moa")
  write.csv(all.WV,file=paste0("Exp_",exp,"_","All_adherent_CCLE_adherent_WGV_top_",n.genes,"_genes_",biomarker.type,"_",expr.type,"_corr_AUC_NES_",type,".csv"))
  all.corr <- rbindlist(corr.list, use.names=TRUE, idcol="moa")
  write.csv(all.corr,file=paste0("Exp_",exp,"_","All_adherent_correlations_top_",n.genes,"_genes_",biomarker.type,"_",expr.type,"_corr_WGV_vs_AUC_w_moa_",type,".csv"))
  all.DMEA.results <- rbindlist(all.DMEA.list, use.names=TRUE, idcol="moa")
  write.csv(all.DMEA.results,file=paste0("Exp_",exp,"_","All_adherent_DMEA_top_",n.genes,"_genes_",biomarker.type,"_",expr.type,"_corr_WGV_",type,".csv"))
  write.csv(subset.self.DMEA,file=paste0("Exp_",exp,"_","All_adherent_validation_drugSEA_top_",n.genes,"_genes_",biomarker.type,"_",expr.type,"_corr_WGV_",type,".csv"))
}

#run 5 experiments for cross-validation
training.samples <- c(train1, train2, train3, train4, train5)
n.genes <- c(1, 5, 10, 25, 50, 100, 250, 500)
corr.types <- c("spearman","pearson")
types.genes <- c("min","max","min&max")
paths.outputs <- c(path.outputs.RNAseq, path.outputs.proteomic)
biomarker.types <- c("RNAseq","proteomic")
biomarker.df <- list(RNA.df, prot.df)
for(k in 1:length(biomarker.types)){
  for(i in 1:5){ #5 experiments
    setwd(paths.outputs[k])
    DMEA.biomarkers(overlap.CCLE=overlap.CCLE, biomarker.data=biomarker.df[[k]], biomarker.type=biomarker.types[k], moa=moa, all.DMEA=all.DMEA, training.samples=training.samples[i], exp=i)
    #try DMEA with Pearson or Spearman estimates and various number of genes
    all.biomarker.corr <- read.csv(file=paste0("Exp_",i,"_all_adherent_correlations_",biomarker.type.list[k],"_vs_AUC_NES_w_moa.csv"))
    all.biomarker.corr$X <- NULL
    for(j in 1:length(n.genes)){
      for(m in 1:length(corr.types)){
        for(l in 1:length(types.genes)){
          print(paste0("Running DMEA for ",n.genes[j]," using the ",types.genes[l]," ",corr.types[m]," correlated genes"))
          setwd(paste0(paths.outputs[k],n.genes[j],"_genes/",corr.types[m],"/",types.genes[l],"/"))
          DMEA.biomarker.validation(all.biomarker.corr=all.biomarker.corr, biomarker.type=biomarker.types[k], n.genes=n.genes[j], type.genes=type.genes[l], expr.data=biomarker.df[[k]],expr.type=biomarker.types[k],moa=moa,overlap.CCLE=overlap.CCLE,type=corr.types[m],training.samples=training.samples[i],exp=i)
        }
      }
    }
  }
}