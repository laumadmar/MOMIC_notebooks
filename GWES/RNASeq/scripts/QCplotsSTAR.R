########################################################################################
###### Takes Log.final.out from STAR output and generates QC plots covariables.   ######
###### Contrast are also defined by the user as input parameters.                 ######
###### Heatmap, pca and volcano plot are created. [Part of multi-omics pipeline]  ######
###### Author: Laura Madrid Marquez Sep 2018                                      ######
###### Notes: copy and rename the file to do your own modifications               ######
########################################################################################


# library(knitr) # for turning the original script into a nicely formatted pdf
opts_chunk$set(echo = TRUE, message = FALSE) # knitr options
opts_knit$set(self.contained = FALSE) # knitr options
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra)) # for composite plotting of ggplots
suppressPackageStartupMessages(library(grid))
options(stringsAsFactors = FALSE)


# ################################################################################
# The functions of this script are not very generalized.
# Most of them expect a long data frame with exactly the following format:
# head(align.results.df)
# V1          V2                                                  replicate sample
# SNF2_1.1 Mapping speed, Million of reads per hour      844.70         1   SNF2
# SNF2_1.2                    Number of input reads 11966526.00         1   SNF2
# SNF2_1.3                Average input read length       51.00         1   SNF2
# SNF2_1.4             Uniquely mapped reads number 10289736.00         1   SNF2
# SNF2_1.5                  Uniquely mapped reads %       85.99         1   SNF2
# SNF2_1.6                    Average mapped length       50.73         1   SNF2
##################################################################################
PlottingCorrelation <- function(DF, Var1, Var2, Var1.Label, Var2.Label){
  # Convenience function for the simple plot() function that allows for separate
  # definition of labels and columns that should be compared against each other
  # usage: PlottingCorrelation(DF=aligned.reads.df,
  #               Var1="Number of input reads", Var2="Uniquely mapped reads %",
  #               Var1.Label = "InputReads", Var2.Label="UniquelyMappedFraction")
  m <- matrix(data = c(DF$V2[which(DF$V1 == Var1)],
                       DF$V2[which(DF$V1 == Var2)]),
              ncol = 2)
  colnames(m) <- c(Var1.Label,Var2.Label)
  plot(m)
}

PlottingAlignmentResults <- function(Filter, DF, Legend=TRUE, PlotMedian = TRUE){
  # this function extracts those lines that correspond to the value stored in Filter and generates a bar plot where each replicate is shown with a different color

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

# there is no guide-box. This fails
ExtractLegend <- function(Plot){
  G <- ggplotGrob(Plot)$grobs
  Legend <- G[[which(sapply(G, function(x) x$name) == "guide-box")]]
  Lheight <- sum(Legend$height)
  return(list(legend = Legend, lheight=Lheight))
}

Extract.Histo.Info <- function(InList, # output from hist(...,plot=FALSE)
                               Percentage = TRUE
){
  # this function uses the information from hist() that are stored in lists
  # to make a data frame that's suitable for bar plots
  ll <- length(InList$breaks)
  out.df <- data.frame(breaks1 = InList$breaks[c(1:ll-1)],
                       breaks2 = InList$breaks[c(2:ll)],
                       counts = InList$counts)
  if(Percentage){
    out.df <- transform(out.df,
                        breaks = paste(out.df$breaks1*100, "-", out.df$breaks2*100, sep = ""),
                        breaks1 = NULL, breaks2 = NULL)
  }else{
    out.df <- transform(out.df,
                        breaks = paste(out.df$breaks1, "-", out.df$breaks2, sep = ""),
                        breaks1 = NULL, breaks2 = NULL)
  }
  out.df$breaks <- factor(out.df$breaks, levels = unique(as.character(out.df$breaks)), ordered = TRUE)
  return(out.df)
}

## Reading in data

# listing files to be read in 
infiles <- list.files(path="/mnt/Internal_Datasets/IN-07-MICE/IN-07-MICE-RNA/analysis/alignmentAndCountSTAR", # folder with log files
                           pattern=".Log.final.out$",
                           recursive = TRUE,
                           full.names = TRUE)
#infiles

# iterating over the file list with a function to read in the log files
# generates a __list of data frames__
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
names(align.results) =sapply(strsplit(as.character(infiles),"/"), `[`, 8)


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
# There is no replicates - this creates a NA column
#align.results.df$replicate <- as.factor(as.numeric(
#  gsub(".*\\_([0-9]*)\\.[0-9]*", "\\1", row.names(align.results.df))
#))

# check the result - we should have a data frame with 4 columns 
#( or 3 without replicates - fastq file from the same sample, different lane were merge with STAR)
head(align.results.df)


## Making plots

#Numbers are usually more easily digestible and comparable if they are represented visually, especially if one has many samples (although you will probably rarely encounter 48 samples per replicate as for this example data used here).
#However, we don\'t need to visualize every entry from the \texttt{STAR} output since some are redundant, and some even not applicable in our case: 

unique(align.results.df$V1)


#Thus, we first define those QC entries that we are interested in:

filters = c("Number of input reads", "Uniquely mapped reads %", "Number of splices: Total", "Number of splices: Non-canonical")


# Now, let\'s generate some bar plots.
# The following chunk of code is optimized for returning a composite image with four plots and one legend.
# Note that the functions \texttt{PlottingAligmentResults} and \texttt{ExtractLegend} are part of the script that we sourced in the beginning.

#{r barPlots, fig.width = 13, fig.height = 10}
# for each entry in "filters", generate a bar chart
plots <- lapply(filters, function(x) 
  PlottingAlignmentResults(x, align.results.df, Legend = FALSE))

# getting the legend is a bit tricky (kudos to Luce for figuring it out) - it does not work
# my.legend <- ExtractLegend(PlottingAlignmentResults(align.results.df, 
#                                                     Filter = filters[1], 
#                                                     Legend=TRUE))

# combining plots and legend
grid.arrange(arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow=2)
             #,my.legend$legend, nrow=2,
             #heights= unit.c(unit(1, "npc") - my.legend$lheight, my.legend$lheight)
)
