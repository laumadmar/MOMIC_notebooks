##########################################################################################
###### Custom script containing functions required for visualization in Task 2 QC ########
###### This fuctions follow Anderson et al. protocol                              ########
###### Author: Laura Madrid Márquez Aug 2020                                      ########
###### Notes: copy and rename the file to do your own modifications if required   ########
##########################################################################################

###### Missingness ######
# Generates plots to visualize the missingness results. (standard output and pdf)
# 1st parameter: optional - full path to missing data report - with no extension. If not provided it looks for plink.imiss and plink.lmiss in the current directory.
# 2nd parameter: option - boolean indicating whereas plots are sent to standard output 
# 3rd parameter: option - boolean indicating whereas plots are stored as pdf. (In the same directory as the input file)
missingness_function <- function(misspath=NA,std=TRUE,pdf=TRUE)
{
    if (!is.na(misspath)){
        path=dirname(misspath)
        filename=basename(misspath)
    } else{
        path=getwd()
        filename="plink"
    }
    # read data into R
    indmiss<-read.table(file=paste(path,paste0(filename,".imiss"),sep='/'), header=TRUE)
    snpmiss<-read.table(file=paste(path,paste0(filename,".lmiss"),sep='/'), header=TRUE)
   
    if(std){
        # Plot hist to standard output
        hist(indmiss[,6],main="Histogram individual missingness", xlab="Missing call rate") #selects column 6, F_MISS
        hist(snpmiss[,5],main="Histogram SNP missingness", xlab="Missing call rate")
    }
    if(pdf){
    # Create pdf
        pdf(paste(path,paste(filename,"histmiss.pdf",sep="-"),sep="/")) #indicates pdf format and gives title to file
        hist(indmiss[,6],main="Histogram individual missingness", xlab="Missing call rate")
        hist(snpmiss[,5],main="Histogram SNP missingness", xlab="Missing call rate")
        dev.off() # shuts down the current device
    }
}

###### SNPs Missingness ######
# Generates plots to visualize the missingness results. (standard output and pdf)
# 1st parameter: optional - full path to missing data report - with no extension. If not provided it looks for plink.imiss and plink.lmiss in the current directory.
# 2nd parameter: option - boolean indicating whereas plots are sent to standard output 
# 3rd parameter: option - boolean indicating whereas plots are stored as pdf. (In the same directory as the input file)
snp_missingness_function <- function(misspath=NA,std=TRUE,pdf=TRUE)
{
    if (!is.na(misspath)){
        path=dirname(misspath)
        filename=basename(misspath)
    } else{
        path=getwd()
        filename="plink"
    }
    # read data into R
    snpmiss<-read.table(file=paste(path,paste0(filename,".lmiss"),sep='/'), header=TRUE)
   
    if(std){
        # Plot hist to standard output
        hist(snpmiss[,5],main="Histogram SNP missingness", xlab="Missing call rate")
    }
    if(pdf){
    # Create pdf
        pdf(paste(path,paste(filename,"histlmiss.pdf",sep="-"),sep="/")) #indicates pdf format and gives title to file
        hist(snpmiss[,5],main="Histogram SNP missingness", xlab="Missing call rate")
        dev.off() # shuts down the current device
    }
}

###### Individuals Missingness ######
# Generates plots to visualize the missingness results. (standard output and pdf)
# 1st parameter: optional - full path to missing data report - with no extension. If not provided it looks for plink.imiss and plink.lmiss in the current directory.
# 2nd parameter: option - boolean indicating whereas plots are sent to standard output 
# 3rd parameter: option - boolean indicating whereas plots are stored as pdf. (In the same directory as the input file)
ind_missingness_function <- function(misspath=NA,std=TRUE,pdf=TRUE)
{
    if (!is.na(misspath)){
        path=dirname(misspath)
        filename=basename(misspath)
    } else{
        path=getwd()
        filename="plink"
    }
    # read data into R
    indmiss<-read.table(file=paste(path,paste0(filename,".imiss"),sep='/'), header=TRUE)
   
    if(std){
        # Plot hist to standard output
        hist(indmiss[,6],main="Histogram individual missingness", xlab="Missing call rate") #selects column 6, F_MISS
    }
    if(pdf){
    # Create pdf
        pdf(paste(path,paste(filename,"histimiss.pdf",sep="-"),sep="/")) #indicates pdf format and gives title to file
        hist(indmiss[,6],main="Histogram individual missingness", xlab="Missing call rate")
        dev.off() # shuts down the current device
    }
}

###### Sex discrepancy ######
# Generates plots to visualize the sex-check results
# 1st parameter: optional - full path to sexcheck file. If not provided it looks for a plink.sexcheck file in the current directory.
# Output pdf is stored in the same directory as the input file
gender_check_function <- function(filepath=NA){
    if (is.na(filepath)){
        filename="plink.sexcheck"
        filepath=paste(getwd(),filename,sep="/")
    }
    basepath=dirname(filepath)
    gender <- read.table(filepath, header=T,as.is=T)
    male=subset(gender, gender$PEDSEX==1)
    female=subset(gender, gender$PEDSEX==2)
    
    # Plot hist to standard output
    hist(gender[,6],main="Gender", xlab="F (Inbreeding coefficient)\nF<0.2: female calls; F>0.8: male calls")
    hist(male[,6],main="Men",xlab="F (Inbreeding coefficient)\nF<0.2: female calls; F>0.8: male calls")
    hist(female[,6],main="Women",xlab="F (Inbreeding coefficient)\nF<0.2: female calls; F>0.8: male calls")

    pdf(paste(basepath,"Gender_check.pdf",sep="/"))
    hist(gender[,6],main="Gender", xlab="F (Inbreeding coefficient)\nF<0.2: female calls; F>0.8: male calls")
    dev.off()
    
    pdf(paste(basepath,"Men_check.pdf",sep="/"))
    hist(male[,6],main="Men",xlab="F (Inbreeding coefficient)\nF<0.2: female calls; F>0.8: male calls")
    dev.off()
  
    pdf(paste(basepath,"Women_check.pdf",sep="/"))
    hist(female[,6],main="Women",xlab="F (Inbreeding coefficient)\nF<0.2: female calls; F>0.8: male calls")
    dev.off()
}

###### MAF ######
# Generates a plot of the MAF distribution.
# 1st parameter: optional - full path to freq file. If not provided it looks for a MAF_check.frq file in the current directory.
# Output pdf is stored in the same directory as the input file
maf_freq_function <- function(filepath=NA){
    if(is.na(filepath)){
        filepath=paste(getwd(),"MAF_check.frq",sep="/")
    }
    basepath=dirname(filepath)
    maf_freq <- read.table(filepath, header =TRUE, as.is=T)
    
    # Plot hist to standard output
    hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
    
    pdf(paste(basepath,"MAF_distribution.pdf",sep="/"))
    hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
    dev.off()
}

###### HWE ######
# Generates a histogram for HWE results.
# 1st parameter: optional - full path to hwe report. If not provided it looks for a plink.hwe file in the current directory.
# 2nd parameter: optional - full path to hew zoom file (SNPs with HWE p-value below 0.00001)
# Output pdf is stored in the same directory as the input file
hwe_function <- function(filepath=NA, zoomfilepath=NA){
    if(is.na(filepath)){
        filepath=paste(getwd(),"plink.hwe",sep="/")
    }
    if(is.na(zoomfilepath)){
        zoomfilepath=paste(getwd(),"plinkzoomhwe.hwe",sep="/")
    }
    basepath=dirname(filepath)
    
    hwe<-read.table(filepath, header=TRUE)
    # Plot hist to standard output
    hist(hwe[,9],xlab="HWE exact p-value",main="Histogram HWE")
    
    pdf(paste(basepath,"histhwe.pdf",sep="/"))
    hist(hwe[,9],xlab="HWE exact test p-value",main="Histogram HWE")
    dev.off()
    
    hwe_zoom<-read.table(zoomfilepath, header=TRUE)
    hist(hwe_zoom[,9],xlab="HWE exact test p-value",main="Histogram HWE: strongly deviating SNPs only")

    pdf(paste(basepath,"histhwe_below_theshold.pdf",sep="/"))
    hist(hwe_zoom[,9],xlab="HWE exact test p-value",main="Histogram HWE: strongly deviating SNPs only")
    dev.off()
}

###### Heterozygosity ######
# Generates a plot of the heterozygosity rate distribution
# 1st parameter: optional - full path to het report. If not provided it looks for a R_check.het file in the current directory.
# Output pdf is stored in the same directory as the input file
het_function <- function(filepath=NA){
    if(is.na(filepath)){
        filepath=paste(getwd(),"R_check.het",sep="/")
    }
    basepath=dirname(filepath)
    het <- read.table(filepath, head=TRUE)
    het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
    # Plot hist to standard output
    hist(as.numeric(het$HET_RATE), xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")

    pdf(paste(basepath,"heterozygosity.pdf",sep="/"))
    hist(as.numeric(het$HET_RATE), xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
    dev.off()
}

# Generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
# 1st parameter: optional - full path to het report. If not provided it looks for a R_check.het file in the current directory.
# This outputs the file fail-het-qc.txt
het_outliers_list_function <- function(filepath=NA) {
    if(is.na(filepath)){
        filepath=paste(getwd(),"R_check.het",sep="/")
    }
    basepath=dirname(filepath)
    het <- read.table(filepath, head=TRUE)
    het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
    het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
    het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
    write.table(het_fail, paste0(basepath,"/fail-het-qc.txt"), row.names=FALSE)
} 

###### Relatedness ######
# Generate a plot to assess the type of relationship.
# 1st parameter: optional - full path to pihat_min.genome. If not provided it looks this file in the current directory.
# Output pdf is stored in the same directory as the input file
relatedness_function <- function(filepath=NA){
    if(is.na(filepath)){
        filepath=paste(getwd(),"pihat_min.genome",sep="/")
    }
    basepath=dirname(filepath)
    relatedness = read.table(filepath, header=T)
    
    hist(relatedness[,10],breaks=50,main="Histogram relatedness", xlab= "Pihat")
    
    pdf(paste(basepath,"hist_relatedness.pdf",sep="/"))
    relatedness = read.table(filepath, header=T)
    hist(relatedness[,10],breaks=50,main="Histogram relatedness", xlab= "Pihat")
    dev.off()
}

###### Relatedness ######
# Generate a plot to assess the type of relationship.
# 1st parameter: optional - full path to pihat_min.genome. If not provided it looks this file in the current directory.
# 2nd parameter: optional - full path to zoom_pihat.genome. If not provided it looks this file in the current directory.
# Output pdf is stored in the same directory as the input file
relatedness_zoom_function <- function(filepath=NA, zoomfilepath=NA){
    if(is.na(filepath)){
        filepath=paste(getwd(),"pihat_min.genome",sep="/")
        zoomfilepath=paste(getwd(),"zoom_pihat.genome",sep="/")
    }
    basepath=dirname(filepath)
    relatedness = read.table(filepath, header=T)
    relatedness_zoom = read.table(zoomfilepath, header=T)
    
    # Plot to standard output
    par(pch=16, cex=1)
    with(relatedness,plot(Z0,Z1, xlim=c(0,1), ylim=c(0,1), type="n",xlab = "Z0 P(IBD=0)",ylab="Z1 P(IBD=1)", main="Relatedness"))
    with(subset(relatedness,RT=="PO") , points(Z0,Z1,col=4))
    with(subset(relatedness,RT=="UN") , points(Z0,Z1,col=3))
    #legend(1,1, xjust=1, yjust=1, legend=levels(as.factor(relatedness$RT)), pch=16, col=c(4,3))
    legend(1,1, xjust=1, yjust=1, legend=c("Parent-Offspring","Unrelated"), pch=16, col=c(4,3))
         
    # Plot to standard output
    par(pch=16, cex=1)
    with(relatedness_zoom,plot(Z0,Z1, xlim=c(0,0.02), ylim=c(0.98,1), type="n", main="Parent-Offspring relations"))
    with(subset(relatedness_zoom,RT=="PO") , points(Z0,Z1,col=4))
    with(subset(relatedness_zoom,RT=="UN") , points(Z0,Z1,col=3))
    #legend(0.02,1, xjust=1, yjust=1, legend=levels(as.factor(relatedness$RT)), pch=16, col=c(4,3))
    legend(0.02,1, xjust=1, yjust=1, legend=c("Parent-Offspring","Unrelated"), pch=16, col=c(4,3))
    
    hist(relatedness[,10],main="Histogram relatedness", xlab= "Pihat")
    
    # Plot to pdf
    pdf(paste(basepath,"relatedness.pdf",sep="/"))
    par(pch=16, cex=1)
    with(relatedness,plot(Z0,Z1, xlim=c(0,1), ylim=c(0,1), type="n",xlab = "Z0 P(IBD=0)",ylab="Z1 P(IBD=1)", main="Relatedness"))
    with(subset(relatedness,RT=="PO") , points(Z0,Z1,col=4))
    with(subset(relatedness,RT=="UN") , points(Z0,Z1,col=3))
    #legend(1,1, xjust=1, yjust=1, legend=levels(as.factor(relatedness$RT)), pch=16, col=c(4,3))
    legend(1,1, xjust=1, yjust=1, legend=c("Parent-Offspring","Unrelated"), pch=16, col=c(4,3))
    dev.off()
    pdf(paste(basepath,"zoom_relatedness.pdf",sep="/"))
    par(pch=16, cex=1)
    with(relatedness_zoom,plot(Z0,Z1, xlim=c(0,0.02), ylim=c(0.98,1), type="n", main="Parent-Offspring relations"))
    with(subset(relatedness_zoom,RT=="PO") , points(Z0,Z1,col=4))
    with(subset(relatedness_zoom,RT=="UN") , points(Z0,Z1,col=3))
    #legend(0.02,1, xjust=1, yjust=1, legend=levels(as.factor(relatedness$RT)), pch=16, col=c(4,3))
    legend(0.02,1, xjust=1, yjust=1, legend=c("Parent-Offspring","Unrelated"), pch=16, col=c(4,3))
    dev.off()
    pdf(paste(basepath,"hist_relatedness.pdf",sep="/"))
    relatedness = read.table(filepath, header=T)
    hist(relatedness[,10],main="Histogram relatedness", xlab= "Pihat")
    dev.off()
}

###### Population Stratification ######

# Generate population stratification plot.
# 1st parameter: full path to MDS dataset.
# 2nd parameter: full path to racefile.
# Output pdf is stored in the same directory as the input data file
MDS_function <- function(filepath, racefilepath){
    basepath=dirname(filepath)
    
    data<- read.table(filepath,header=TRUE)
    race<- read.table(racefilepath,header=TRUE)
    datafile<- merge(data,race,by=c("IID","FID"))

    
    for (i in 1:nrow(datafile))
    {
    if (datafile[i,14]=="EUR") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.1,0.15),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="blue")}
    par(new=T)
    if (datafile[i,14]=="EAS") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.1,0.15),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="green")}
    par(new=T)
        if (datafile[i,14]=="SAS") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.1,0.15),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="yellow")}
    par(new=T)
    if (datafile[i,14]=="AMR") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.1,0.15),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="orange")}
    par(new=T)
    if (datafile[i,14]=="AFR") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.1,0.15),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="purple")}
    par(new=T)
    if (datafile[i,14]=="OWN") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.1,0.15),xlab="MDS Component 1",ylab="MDS Component 2",pch=3,cex=0.7,col="black")}
    par(new=T)
    }
    
    C1vector<-datafile[which(datafile$SuperPop=="EUR"),]$C1
    C2vector<-datafile[which(datafile$SuperPop=="EUR"),]$C2

    abline(v=mean(C1vector),col="RED",lty=2)
    abline(h=mean(C2vector),col="RED",lty=2)

    abline(v=mean(C1vector)-6*sd(C1vector),col="BLUE",lty=2)
    abline(v=mean(C1vector)+6*sd(C1vector),col="BLUE",lty=2)

    abline(h=mean(C2vector)-6*sd(C2vector),col="BLUE",lty=2)
    abline(h=mean(C2vector)+6*sd(C2vector),col="BLUE",lty=2)
    
    thresholds <-c(mean(C1vector)-6*sd(C1vector),mean(C1vector)+6*sd(C1vector),mean(C2vector)-6*sd(C2vector),mean(C2vector)+6*sd(C2vector))
    write.table(thresholds,file=paste(basepath,"thresholds.txt",sep="/"))
    
    #abline(v=-0.035,lty=3)
    #abline(h=0.035,lty=3)
    
   legend("topright", pch=c(1,1,1,1,3),c("EUR","EAS","SAS","AMR","AFR","OWN"),col=c("blue","green","yellow","orange","purple","black"),bty="o",cex=1)

} 

# Identifies individuals that fails the thresholds. That is individuals which components 1 and 2 are outside +/- 6 SD from the mean
# 1st parameter: full path to MDS dataset.
# 2nd parameter: full path to racefile.
# 3nd parameter: full path to threshold file
# 4rd parameter: output file
get_failures <- function(filepath, racefilepath, thresholds, output){
    data<- read.table(filepath,header=TRUE)
    race<- read.table(racefilepath,header=TRUE)
    mergedata<- merge(data,race,by=c("IID","FID"))
    
    OWNset<-subset(mergedata,SuperPop=='OWN')
    thres = read.table(thresholds,header = TRUE)
    # this filter out individuals outsite the thresholds. 
    failures=subset(OWNset,C1<thres[1,] | C1>thres[2,] | C2<thres[3,] | C2>thres[4,], select=c(IID,FID))
    # write individuals that fails the thresholds
    write.table(failures,output,sep=" ",row.names=FALSE, col.names=FALSE, quote=FALSE)
}