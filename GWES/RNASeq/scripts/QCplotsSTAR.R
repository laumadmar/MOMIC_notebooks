########################################################################################
###### Takes Log.final.out from STAR output and generates QC plots covariables.   ######
###### Contrast are also defined by the user as input parameters.                 ######
###### Heatmap, pca and volcano plot are created. [Part of multi-omics pipeline]  ######
###### Author: Laura Madrid Marquez Sep 2018                                      ######
###### Notes: copy and rename the file to do your own modifications               ######
########################################################################################


# library(knitr) # for turning the original script into a nicely formatted pdf
#opts_chunk$set(echo = TRUE, message = FALSE) # knitr options
#opts_knit$set(self.contained = FALSE) # knitr options
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra)) # for composite plotting of ggplots
suppressPackageStartupMessages(library(grid))
options(stringsAsFactors = FALSE)

PlottingAlignmentResultsMain <- function (logs_path,pos=8){
## Reading in data

# listing files to be read in 
infiles <- list.files(path=logs_path, # folder with log files
                           pattern=".Log.final.out$",
                           recursive = TRUE,
                           full.names = TRUE)
#infiles

# iterating over the file list with a function to read in the log files
# generates a list of data frames 
align.results <- lapply(infiles, function(x)
  read.table(x, sep="|", 
             strip.white=TRUE,
             stringsAsFactor=FALSE,
             skip=3, fill = TRUE, header = FALSE) )
typeof(align.results)
head(align.results[[1]]) # peek into one of the data.frames within the list

# removing "%" from some of the values to keep just the numeric parts
align.results <- lapply(align.results, function(x)
  transform(x, V2 = as.numeric(gsub("%", "", x$V2) )))

# some cosmetics of each data frame's name - this is specific for the sample names
# of the files used here!
names(align.results) = sapply(strsplit(as.character(infiles),"/"), `[`, pos)


## Generating a long data frame for ggplot2-based plotting

# concatenating all data frames of align.results
align.results.df <- as.data.frame(do.call(rbind, align.results))

# removing lines without values 
align.results.df <- align.results.df[complete.cases(align.results.df),]

# adding additional columns with information about sample and replicate ID,
# using the information from the row names (which are, in turn, based on the
# names of the individual data frames that were stored in the original 
# list, align.results)
align.results.df$sample <- gsub("(.*)\\_.*", "\\1", row.names(align.results.df))

# If there is no replicates - if not this creates a NA column
#align.results.df$replicate <- as.factor(as.numeric(
#  gsub(".*\\_([0-9]*)\\.[0-9]*", "\\1", row.names(align.results.df))
#))

# check the result - we should have a data frame with 3 columns (or 4 with replicates)
head(align.results.df)

## Making plots
# We don\'t need to visualize every entry from the \texttt{STAR} output since some are redundant, and some even not applicable in our case: 

unique(align.results.df$V1)

#Thus, we first define those QC entries that we are interested in:

filters = c("Number of input reads", "Uniquely mapped reads %", "Number of splices: Total", "Number of splices: Non-canonical")

#{r barPlots, fig.width = 13, fig.height = 10}
# for each entry in "filters", generate a bar chart
plots <- lapply(filters, function(x) 
  PlottingAlignmentResults(x, align.results.df, Legend = FALSE))

# combining plots and legend
grid.arrange(arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow=2)
             #,my.legend$legend, nrow=2,
             #heights= unit.c(unit(1, "npc") - my.legend$lheight, my.legend$lheight)
)
}
       
                
# This function extracts those lines that correspond to the value stored in Filter and generates a bar plot where each replicate is shown with a different color
PlottingAlignmentResults <- function(Filter, DF, Legend=TRUE, PlotMedian = TRUE){
  filtered.df <- DF[which(DF$V1 == Filter),]
  medians <- as.data.frame(aggregate(V2~sample, data=filtered.df, FUN=median))
  filtered.df <- merge(filtered.df, medians, by.x = "sample", by.y = "sample", all.x=TRUE)
  
  p <- ggplot(data=filtered.df, aes(y=V2.x, x=sample)) +
    geom_bar(stat="identity",position=position_dodge()) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=4),
          axis.title.x = element_text(size=6),
          axis.title.y = element_text(size=6),
          legend.position="bottom",
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"),
          legend.title=element_blank()) +
    coord_flip() + ylab("") + ggtitle(Filter)
  
  if(PlotMedian){
    p <- p + geom_errorbar(aes(y=V2.y, ymax=V2.y, ymin=V2.y), linetype="dashed")
  }
  
  if(!Legend){
    p <- p +  theme(legend.position="none")
  }
  
  return(p)
}
