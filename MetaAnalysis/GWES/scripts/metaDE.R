##########################################################################################
###### Custom script for meta analysis with metaDE R library                      ########
###### Combines estimators from differential analysis results                     ########
###### Author: Laura Madrid Márquez Sep 2018                                      ########
##########################################################################################

suppressPackageStartupMessages(library(MetaDE))
source("scripts/metaDE.ES_custom.R")

# Creates a dataframe cointaing estimators for using in metaDE
# Estimators are Effect Size (ES), Variance (Var) and P value.
# 1st parameter: path where the input dataset is located
# 2nd parameter: filename of input dataset. The one containing limma/differential analysis results
# 3rd parameter: type of dataframe to be prepared. "onlyP" if you want to consider only P values to do the meta analysis. By default it calculates all estimators.
prepare_matrix_function <- function(path,filename,type=NA){
    datafile=paste(path,filename,sep="/")
    dataset <- read.table (datafile,check.names = FALSE) 
    
    # calculates ES and Var if limma dataset as input
    if(is.na(type)){
        dataset$ES<-as.numeric(dataset$logFC)
        stderr<-(as.numeric(dataset$CI.R)-as.numeric(dataset$CI.L))/3.92
        dataset$Var <- stderr^2
        matrix<-dataset[ , c("ES","Var","P.Value")]
        matrix$gene=rownames(dataset)
    } else { # when only P values are available
        matrix<-data.frame(gene=rownames(dataset),P.Value=dataset$minfdr)
    } 
    return (matrix)
}

# Performs a meta analysis using the MetaDE.ES algorithm from metaDE R library. This method is customised
# so that it gives also an estimator when the gene is not present in all datasets.More info in /mnt/data/scripts/MetaDE.ES_custom.R
# It ranks the resulting genes in ascending order of Fisher P value.
# 1st parameter: single dataframe merging datasets obtained from prepare_matrix_function
# 2nd parameter: optional - key name to add to logFC column in result table (e.g case-control)
# 3rd parmeter: optional - path where to store the output file (default to current directory)
# 4rd parameter: optional - output file name 
meta_function <- function(allmatrix,keyname=NA,output.path=NA,output.file=NA){
    # ES matrix with cbind
    EScolumns=grep("ES", colnames(allmatrix)) # index of columns containing ES in their name
    ESmatrix<-allmatrix[, EScolumns]
    rownames(ESmatrix)<-allmatrix$gene
    #names(ESmatrix)<-names(allmatrix[EScolumns])
    options(width=120); print(head(ESmatrix,n=3))
    
    # Var matrix with cbind
    VARcolumns=grep("Var", colnames(allmatrix)) # index of columns containing Var in their name
    VARmatrix=allmatrix[, VARcolumns] # to create and inizialise the matrix to the same number of rows as allmatrix
    rownames(VARmatrix)<-allmatrix$gene
    options(width=120); print(head(VARmatrix,n=3))
    
     # Meta custom DE function
    cat("\n\n* doing Meta DE\n")
    ind.res <- list(ESmatrix,VARmatrix)
    names(ind.res) <- c("ES","Var")
    metaES <- MetaDE.ES.custom(ind.res,meta.method='REM')
    
    # P matrix with cbind
    Pcolumns=grep("P.Value", colnames(allmatrix)) # index of columns containing P.value in their name
    Pmatrix=allmatrix[, Pcolumns] # to create and inizialise the matrix to the same number of rows as allmatrix
    rownames(Pmatrix)<-allmatrix$gene
    options(width=120); print(head(Pmatrix,n=3))
    
    # Meta P function
    cat("\n\n* doing Meta P\n")
    Pind.res <- list(Pmatrix)
    names(Pind.res) <- c("p")
    metaP<-MetaDE.pvalue(Pind.res,meta.method='Fisher')
    
    #Prepare final dataset
    final <- data.frame(metaES$mu.hat,metaES$mu.var,metaES$Qpval,metaES$pval,metaES$FDR,metaP$meta.analysis$pval,metaP$meta.analysis$FDR,metaES$n)
    names(final) <-  c(paste("logFC",keyname,sep="."),"Var","Qpvalue","REM.Pvalue","REM.FDR","Fisher.Pvalue","Fisher.FDR","n estimators")
    final <- final[order(final$Fisher.Pvalue), ]

    rank<-rank(final$Fisher.Pvalue,na.last = TRUE, ties.method = "min") # rank results attending to fisher p value asc
    ranked<-cbind(rank,final)
    NonNAfisher <- which(!is.na(ranked$Fisher.Pvalue))
    lastNonNA <- max(NonNAfisher) #assign same rank if Fisher p value is NA
    ranked$rank[is.na(ranked$Fisher.Pvalue)]<-lastNonNA+1
    options(width=120); print(head(ranked))
    
    # create a meaningfull output file name if not specified
    if(is.na(output.file)){ 
        output.file="meta_result"
        if(!is.na(keyname)){
            output.file=paste("meta_result",keyname,sep="_")
        }
    }  
    # output will be stored in current directory if not specified
    if(is.na(output.path)){
        output.path=getwd()
    } 
    cat("\n\n * writing to ",output.path)
    write.table(ranked,paste(output.path,output.file,sep="/"))
}

# Performs a meta analysis using the MetaDE.pvalue algorithm from metaDE R library. 
# It ranks the resulting genes in ascending order of Fisher P value.
# 1st parameter: single dataframe merging datasets obtained from prepare_matrix_function
# 2nd parameter: optional, key name to add to logFC column in result table
# 3rd parmeter: path where to store the output file (default to current directory)
# 4rd parameter: optional - output file name 
meta_P_function <- function(allmatrix,keyname=NA,output.path=NA,output.file=NA){
    # P matrix with cbind
    Pcolumns=grep("P.Value", colnames(allmatrix)) # index of columns containing P.value in their name
    Pmatrix=allmatrix[, Pcolumns] # to create and inizialise the matrix to the same number of rows as allmatrix
    rownames(Pmatrix)<-allmatrix$gene
    options(width=120); print(head(Pmatrix))
    
    # Meta P function
    cat("\n\n* doing Meta P\n")
    Pind.res <- list(Pmatrix)
    names(Pind.res) <- c("p")
    metaP<-MetaDE.pvalue(Pind.res,meta.method='Fisher')
    
    #Prepare final dataset
    final <- data.frame(metaP$meta.analysis$pval,metaP$meta.analysis$FDR)
    names(final) <-  c("Fisher.Pvalue","Fisher.FDR")
    final <- final[order(final$Fisher.Pvalue), ]

    rank<-rank(final$Fisher.Pvalue,na.last = TRUE, ties.method = "min") # rank results attending to fisher p value asc
    ranked<-cbind(rank,final)
    NonNAfisher <- which(!is.na(ranked$Fisher.Pvalue))
    lastNonNA <- max(NonNAfisher) #assign same rank if Fisher p value is NA
    ranked$rank[is.na(ranked$Fisher.Pvalue)]<-lastNonNA+1
    options(width=120); print(head(ranked))
    
    # create a meaningfull output file name if not specified
    if(is.na(output.file)){ 
        output.file="meta_result"
        if(!is.na(keyname)){
            output.file=paste("meta_result",keyname,sep="_")
        }
    }  
    # output will be stored in current directory if not specified
    if(is.na(output.path)){
        output.path=getwd()
    } 
    cat("\n\n * writing to ",output.path)
    write.csv(ranked,paste(output.path,output.file,sep="/"))
 }   
