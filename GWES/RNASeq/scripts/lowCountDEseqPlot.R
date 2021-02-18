########################################################################################
###### Custom script that plots the pct of low counts after DESeq                 ######
###### Uses summary(DGE.results)                                                  ######
###### Author: Laura Madrid Márquez Aug 2020                                      ######
###### Notes: copy and rename the file to do your own modifications if required   ######
########################################################################################


pct_low_counts_function <- function(summarypath=NA){    
    # listing files to be read in 
    infiles <- list.files(path=summarypath, # folder with log files
                          pattern=".summary$",
                          full.names = TRUE)

#infiles

results <- lapply(infiles, function(x)
  read.table(x, sep=":", 
             strip.white=TRUE,
             stringsAsFactor=FALSE,
             skip=3, fill = TRUE, header = FALSE) )

#head(results[[1]]) # peek into one of the data.frames within the list


# take pct of low count and remove "%" to keep just the numeric parts
results <- lapply(results, function(x)
  transform(x, V3 = as.numeric(gsub("%", "", sapply(strsplit(V2,","), `[`, 2))) ))

# plot (pct in column V3 row 4)
pct=unlist(lapply(results, function(x)
  x[[4,3]]))
labels=gsub(".summary", "",sapply(strsplit(infiles,"/"), tail, 1))
#barplot(pct,names.arg=labels,las=3,cex.names=0.4,main = "pct low counts")

#pdf(file=paste(summarypath,"plot.pdf",sep="/"), onefile=T, paper='A4r') 
a=barplot(pct,main = "pct low counts")
par(mar = c(9,4,4,2) + 0.1)
text(a[,1], -3.7, srt = 60, adj= 1, xpd = TRUE, labels = labels , cex=0.55)
#dev.off() 
}