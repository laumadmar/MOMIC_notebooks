##########################################################################################
###### Custom script for differential analysis with limma, allowing subsets and   ########
###### covariables. Contrast are also defined by the user as input parameters.    ########
###### Heatmap, pca and volcano plot are created. [Part of multi-omics pipeline]  ########
###### Author: Laura Madrid Marquez Sep 2018                                      ########
###### Notes: copy and rename the file to do your own modifications               ########
##########################################################################################

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(NMF))

# Reads clinical and expression datasets, subset by phenotype indicated as input parameter and call limma function.
# 1st parameter: path where the clinical and expression sets are. Optional is the input datasets are R dataframe objects.
# 2nd parameter: clinicalset dataframe or file name that points to the dataframe in the provided path
# 3rd parameter: expression dataframe or file name that points to the dataframe in the provided path. Probes as rows, samples as columns
## 4-7th parameters: parameters used in limma for Diferential Expression. 4th and 5th will be used as contrast: condition1 - condition2 (i.e. case-control)
# 4th: first condition that satisfied the first group your are comparing. i.e "case"
# 5th: second condition that satisfied the second group your are comparing. i.e "control"
# 6th: column name from the clinical dataframe containing the condition 
# 7th: character string separated by commas indicating other variables (column names from clinical data) to include as covariates in limma 
## 8-9th parameters are optional and needed if you want to do diff analysis on a subset of the original data (i.e. only in males)
# 8th parameter: column name from the clinical dataframe to be used to subset both clinical and expr dataframes
# 9th parameter: value of column name from the clinical df to be used to subset
# 10th parameter: optional - output path to store limma table as a result (default to current directory)
# 11th parameter: optional - output file name
# 12th parameter: platform/chip 
# 13th parameter: optional - number of probes to be represented in the heatmap plot
# DELETED column.match parameter: column name from clinicalset that matches colnames (samples) in expression dataset to be used to subset
diff_analysis_function <- function(input.path=NA,clinicalset,expression,DE.cond1,DE.cond2,DE.variable, DE.covariates=NA,subset.variable=NA,subset.value=NA,output.path=NA,out.name=NA,platform,nheatmap=50){
    cat("\n#### ",DE.cond1," versus ",DE.cond2," ####\n")
    if(!is.na(subset.variable)){ cat("\n## Subset ",paste(subset.variable,subset.value,sep=" ")," ##\n")}
    
    # read datasets
    if(class(clinicalset)=="character"){
        clinicalset=read.table(paste(input.path,clinicalset,sep="/"),check.names = F)
    }   
    if(class(expression)=="character"){
        expression=read.table(paste(input.path,expression,sep="/"),check.names = F)
    }
    cat("\ndim(clinicalset) ",dim(clinicalset))
    cat("\ndim(expr) ",dim(expression))
    
    # subset original data keeping only subjects that belongs to cond1 and cond2
    cat("\nSubseting datasets to keep those subjects belonging to provided conditions: ",DE.cond1," and ",DE.cond2)
    clinicalset=clinicalset[ clinicalset[[DE.variable]] %in% c(DE.cond1,DE.cond2), ]
    cat("\ndim(clinical.subset) ",dim(clinicalset))

    #expression=expression[, colnames(expression) %in% clinicalset[[column.match]]]
    expression=expression[, colnames(expression) %in% rownames(clinicalset)]
    cat("\ndim(expr.subset) ",dim(expression))

    # extra optional subseting by variable
    if(!is.na(subset.variable)){
        cat("\nSubseting datasets by ",subset.value)
        clinicalset = subset(clinicalset,clinicalset[[subset.variable]]==subset.value)
        cat("\ndim(clinical.subset) ",dim(clinicalset))

        #expression=expression[, colnames(expression) %in% clinicalset[[column.match]]]
        expression=expression[, colnames(expression) %in% rownames(clinicalset)]
        cat("\ndim(expr.subset) ",dim(expression))
    }
    
    # if both input and output path are not specified, the output will be stored in current directory under out.name name
    if(is.na(out.name)){ # create a meaningfull output file name if not specified
        out.name=paste("limma",DE.cond1,DE.cond2,sep="_")
        if(!is.na(subset.value)){
            out.name=paste(out.name,subset.value,sep="_")
        }
    }
    
    if(!is.na(output.path)){
        out.name=paste(output.path,out.name,sep="/")
    } else if(!is.na(input.path)){
        out.name=paste(input.path,out.name,sep="/")
    } 
    # call limma function and plot heatmap and pca.        
    top=limma_param_function(clinicalset,expression,DE.cond1,DE.cond2,DE.variable,DE.covariates,out.name,platform,nheatmap)
    heatmapset<-expression[rownames(top) ,] #top 50 DE probes
    heatmap_function(heatmapset,sapply(strsplit(out.name,"/"), tail, 1))
    #pca_function(expression,clinicalset,DE.variable) 
    
}

# limma function to do diff expression. It returns a limma dataset (toptable). 
# 1st param: clinical/metadata dataset
# 2nd param: expression dataset
# 3rd: first condition that satisfied the first group your are comparing. i.e "case"
# 4th: second condition that satisfied the second group your are comparing. i.e "control"
# 5th: column name from the clinical dataframe containing the condition 
# 7th: character string separate by commas indicating other variables (column names from clinical data) to include as covariates in limma
# 8rd parm: output file name
# 9th param: number of rows to be returned
limma_param_function<-function(clinicalset, expresion_df,DE.cond1,DE.cond2,DE.variable,DE.covariates=NA,filename,platform,number){ 
    cat("\n\n== Doing limma == \n")
    # design limma table
    cov.model=paste("~0 + clinicalset[[DE.variable]]")
    #if(!is.na(DE.covariates)){ 
    if(length(DE.covariates)!=0){ 
        ncov=length(DE.covariates)
        for (i in 1:ncov){
            cov.model=paste(cov.model,paste("clinicalset[[DE.covariates[",i,"]]]",sep=""),sep="+")
        }
    }  
    cat("\nmodel formula: ", cov.model)
    design.matrix = model.matrix(as.formula(cov.model))
    colnames(design.matrix)=make.names(colnames(design.matrix)) # make.names removes strange characters
    colnames(design.matrix)[grepl(DE.cond1,colnames(design.matrix))]=DE.cond1
    colnames(design.matrix)[grepl(DE.cond2,colnames(design.matrix))]=DE.cond2
    cat("\n\ndesign matrix: \n")
    print(head(design.matrix,n=3))
    
    # limma DE if there are enough samples on each condition (i.e >0)
    topn=NA
    if(sum(design.matrix[,1]==1)>0 && sum(design.matrix[, 2]==1)>0) {
        contrast=paste(eval(DE.cond1),eval(DE.cond2),sep="-")
        cat("\n\ncontrast: ",contrast,"\n")
        contrast.matrix <- makeContrasts(contrasts = contrast, levels=design.matrix)
        #cat("\n\ncontrast matrix: ",dim(contrast.matrix))
        fit = lmFit(expresion_df, design.matrix)
        fit <- contrasts.fit(fit, contrast.matrix)
        fit <- eBayes(fit)
        
        toplist=topTable(fit,sort.by="P",confint=TRUE)
        cat("\n\nlimma topTable: \n")
        options(width=120); print(toplist)  

        topn=topTable(fit,sort.by="P",n=number)

        alltable=topTable(fit,sort.by="P",n=Inf,confint=TRUE)
        cat("\n\nDE summary: \n")
        print(summary(decideTests(fit)))
        write.table(alltable,filename)
        
        alltable.annot<-annotate_function(alltable,platform)
        write.table(alltable.annot,paste(filename,"annot",sep="_"))
        volcano_plot_function(alltable.annot)
        cat("\n\nannotated limma toptable: \n")
        print(head(alltable.annot,n=10))

    } else {cat("\n\n*Not enough sample on for any of the conditions*")}    
    return (topn) # for heatmap
}

# plot a coloured volcano
# 1st parameter: limma table
volcano_plot_function<-function(dataset){
    # volcano plot
    with(dataset, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", ylab="-log10 p-value", xlab="log2 fold change"))

    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(dataset, adj.P.Val<.05 ), points(logFC, -log10(P.Value), pch=20, col="firebrick"))
    with(subset(dataset, abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="darkorange4"))
    set<-subset(dataset, adj.P.Val<.05 & abs(logFC)>1)
    with(set, points(logFC, -log10(P.Value), pch=20, col="forestgreen"))

    # Label points with the textxy function from the calibrate plot
    suppressPackageStartupMessages(library(calibrate))
    with(set, textxy(logFC, -log10(P.Value), labs=rownames(set), cex=.4))
    # if I use as labels rownames(alltable.annot) it is ilegible
}    
# translate probe IDs to gene symbols
# 1st parameter: limma table
# 2nd parameter: chip
annotate_function<-function(dataset,platform){
    probes <- rownames(dataset)
    annot.path="/mnt/data/GWES/Microarray/ref_files/annot_files"
    if(platform=="570"){
        annot=read.table(paste(annot.path,"GPL570-platform.txt",sep="/"),header=TRUE, sep="\t")
        names(annot)<-c("platform","PROBEID","gene_symbol","gene_id","organism")
    }  else if (platform=="2700"){
        annot=read.csv(paste(annot.path,"GPL2700-platform.csv",sep="/"),header=TRUE, sep="\t")
        names(annot)<-c("platform","PROBEID","gene_symbol","gene_id","organism")
    } else if(platform=="affy"){
        suppressPackageStartupMessages(library(hgu133plus2.db))
        suppressPackageStartupMessages(library(hgu133a.db))
        annot <- select(hgu133plus2.db, probes, c("SYMBOL", "ENTREZID")) #c("SYMBOL")
        # in case of hgu133 there exists another platform - merge both first
        annot2<- select(hgu133a.db, probes, c("SYMBOL", "ENTREZID"))
        annot<-merge(annot,annot2,all=T)
        names(annot)<-c("PROBEID","gene_symbol","ENTREZID")
    }  else if (platform=="illumina"){
        annot=read.csv(paste(annot.path,"HumanHT-12_V3_0_R3_11283641_A.csv",sep="/"),header=TRUE, sep="\t")
        names(annot)<-c("ENTREZID","gene_symbol","PROBEID")
    }  else if (platform=="agilent"){
        annot.file=read.csv(paste(annot.path,"SurePrintG3HumanGEv38x60KMicroarrayGeneList072363_D_GeneList_20150612.csv",sep="/"),sep=",",header=T)
        annot=annot.file[, c("ProbeID","GeneSymbol")]
        names(annot)=c("PROBEID","gene_symbol")
    }  
    df.merge<-merge(dataset,annot,by.x="row.names",by.y="PROBEID",all.x=TRUE)
    df.annotated<-df.merge[with(df.merge, order(P.Value)), ]
    df.clean<-df.annotated[!duplicated(df.annotated$gene_symbol), ]
    df.clean<-df.clean[!(is.na(df.clean$gene_symbol) | df.clean$gene_symbol==""), ]
    df.final<-df.clean[, -c(which(names(df.clean)=="gene_symbol"),which(names(df.clean)=="Row.names"))]
    rownames(df.final)=df.clean$gene_symbol
    cat("\n dim alltable annotated ",dim(df.final),"\n")
    
    return (df.final)
}

# Plots a heatmap with the input dataset
heatmap_function<-function(heatmapset,title){
    #library(gplots)
    #heatmap.2(as.matrix(heatmapset), scale = "row", trace="none",cexRow=0.5,cexCol=0.5,main=title)    
    NMF::aheatmap(as.matrix(heatmapset) ,Rowv = TRUE , Colv = TRUE ,distfun = "euclidean", hclustfun = "average",scale = "row",main=title)
}

# Plots a pca (if there are any NA in the expression set it removes them first)
pca_function<- function(expression_df,clinicalset,variable) {
    pca <- prcomp(t(na.omit(expression_df)))
    pca2 <- as.data.frame(pca$x)
    pca2$group <- clinicalset[[variable]] # if group is a number, the legend in plot is a color gradient
        
    theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
    p<-ggplot(pca2,aes(x=PC1,y=PC2,color=group))
    p<-p+geom_point() +theme
    plot(p)
}    

