####GSVA on TLR Signature
4/30/2020

#Scott Jelinsky
#load library
library(GSVA)
library(GSVAdata)
library(GSEABase)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)

###############
#Set up GeneSetColletion objects from Hallmark genesets
gsc <- getGmt("/app/PathwayAnalysis/Data/h.all.v5.2.symbols.gmt",
              collectionType=BroadCollection(category="c3"),
              geneIdType=SymbolIdentifier())
gsc <- getGmt("/app/Shiny/AMP_Phase1_RA/Data/TLR.gmt",
              collectionType=BroadCollection(category="c3"),
              geneIdType=SymbolIdentifier())
exprs = exprs(lowinputGSET)
  row.names(exprs) <- make.names(fData(lowinputGSET)$gene_name, unique = T)
  
gset1 <- new("ExpressionSet", exprs = exprs)
        fData(gset1) <- fData(lowinputGSET)
        pData(gset1) <- pData(lowinputGSET)
    

#remove low level expressed genes
        LowExpRem <- function(exprs=exprs, gset1=gset1){
          maxexpr <- data.frame(apply(exprs,1,max))
          indexll <- data.frame(apply(exprs,1,max) >1)[,1]
          gset1.f <- gset1[indexll,] 
        return(gset1.f)
        }
gset1.f <- LowExpRem(exprs=exprs, gset1=gset1)        
      

#Explore expression of marker genes
library(dplyr)
for (i in 1:7){
  GeneList1 <- gsc@.Data[[i]]@geneIds #get gene list for cluster 1
      xx <- min(apply(exprs[row.names(exprs)%in%GeneList1,][,],1,max) )
      print(xx)
}
  
####All markers have high expression and should be included
  
###Create Unique GENEset collections
        
        
                
###############
#Single cell GSVA
res  <- gsva(expr = exprs(gset1), gsc, min.sz=10, max.sz=500, verbose=TRUE) #KI
res  <- gsva(expr = exprs(gset1.f), gsc, min.sz=10, max.sz=500, verbose=TRUE) #KI

#### Determine Regulated Genesets
  RegGenSet <- function(gset1.f=gset1.f){
    design <- model.matrix(~0+factor(gset1.f$DiseaseTissue))
    colnames(design) <- gsub("factor(gset1.f$DiseaseTissue)", "", colnames(design), fixed = T)
    
      #fit <- lmFit(res, design)
      fit <- eBayes(fit)
    CONtrasts <- c("RA_Mono  - OA_Mono",
                   "RA_Tcell - OA_Tcell",
                   "RA_Bcell - OA_Bcell",
                   "RA_Fibro - OA_Fibro")
    cont.wt1 <- makeContrasts(contrasts=CONtrasts,
                              levels=design )
    fit1.wt <- (contrasts.fit(fit, cont.wt1[,]))  
    efit1 <- eBayes(fit1.wt)
      Mono <- topTable(efit1, coef=1, number=Inf,sort.by = "none")
      TCell <- topTable(efit1, coef=2, number=Inf,sort.by = "none")
      BCell <- topTable(efit1, coef=3, number=Inf,sort.by = "none")
      Fibro <- topTable(efit1, coef=4, number=Inf,sort.by = "none")
   RegSet <- list(Mono=Mono, Tcells=TCell, Bcells=BCell, Fibroblasts=Fibro)
   return(RegSet)
  }

  RegSet <- RegGenSet(gset1.f = gset1.f)
RegSet


        ###############
        #Function to format data for plotting and 
        # to select for specific pathways        
        res.melt.function <- function(ress=res, eset=gset1){
          res.melt <- melt(ress)
          #res.melt$Var2 <- make.names(res.melt$Var2)
          row.names(pData(eset)) <- make.names(row.names(pData(eset)))
          row.names(res.melt) <- make.names(row.names(res.melt))
          res.melt <- merge(res.melt, pData(eset), by.x="Var2", by.y="Sample")
          #index <- grep("COAG|IL6|ALPHA|INFLAMMATORY|MYC_TARGETS_V2|COMPLEMENT", res.melt$Var1)
          res.melt.coag <- res.melt#[index,]
          #index2 <- grep("IL6|_3|_4", res.melt.coag$Group,invert = T)
          #res.melt.coag <- res.melt.coag[index2,]
          return(res.melt.coag)
        }
        
        res.melt.coag <- res.melt.function(ress=res, eset=gset1.f)  #Tofa
        
        ###################
        #setup comparisons for ploting Add Mean Comparison P-Values To A Ggplot    
        my_comparisons <- list( c("RA_Mono", "OA_Mono")
        )
 
        # New facet label names for dose variable
        # Var1.labs <- c("IL6_Jak", "IFNa", "Compl", "MYC", "Inflamm", "Coagulation")
        # names(Var1.labs) <- c( "HALLMARK_IL6_JAK_STAT3_SIGNALING",   "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        #                        "HALLMARK_COMPLEMENT",                "HALLMARK_MYC_TARGETS_V2",           
        #                        "HALLMARK_INFLAMMATORY_RESPONSE",     "HALLMARK_COAGULATION")
        # 
        ###############
        #PLOTing function 
        plotfunction <- function(RES= res.melt.coag){
          require(cowplot)
          require(ggpubr)
          require(ggpval)
            #RES<- RES[RES$Cell.type=="Mono",]
              my_comparisons <- list( c("OA", "RA"))
                a <- ggplot(RES,aes(x=Disease, y=value))+ 
                  geom_boxplot(outlier.shape = NA) + geom_jitter(alpha=0.4)
                a <- a+ facet_grid(Cell.type~Var1)
                
                a
                #a <- a + stat_compare_means(comparisons = my_comparisons,label.y = 0.8) # Add pairwise comparisons p-value
                a <- a + ylab("ssgsva Pathway enrichment score") + 
                  #scale_y_discrete(labels=c("1" = "Baseline", "2" = "3mths")) +
                  theme(axis.text.x = element_text( angle = 90), 
                  axis.title.x = element_blank()) 
                RegMon = dplyr::bind_rows(RegSet, .id = "variable")
                RegMon$p.format <- format(RegMon$P.Value,digits=2)
                RegMon$Group <- rownames(RegSet[[1]])
                labels   <- data.frame(Var1=RegMon$Group, Cell.type= RegMon$name, label=RegMon$p.format, Disease="RA", value=0)
                  labels$Cell.type <- gsub("s|blast", "", labels$Cell.type)
                  labels$Cell.type <- gsub("cell", " cell", labels$Cell.type)
                a <- a + geom_text(data=labels, aes(x = -Inf, y = Inf, label = label, group=Var1),
                          size = 5,
                          hjust = -0.5,
                          vjust = 1.4,
                          inherit.aes = FALSE)
                a <- a+ geom_rect(data = subset(labels, as.numeric(as.character(labels$label)) <0.05),aes(fill = Var1),xmin = -Inf,xmax = Inf,
                          ymin = -Inf,ymax = Inf,alpha = 0.3)
          a
          return(a)
        }
        

        
        
        plotfunction(RES = res.melt.coag)  #KI
        
      
