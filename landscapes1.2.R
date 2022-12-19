#!/usr/bin/env Rscript
# Usage: Rscript --vanilla landscapes1.2.R mode style filename indexname plotname

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 4) {
  stop("At least four arguments must be supplied", call.=FALSE)
} else {
  
  # check mode argument
  if(args[1] %in% c("genome", "chromosome")){
    mode = args[1]
  } else {
    stop("The argument 'mode' can only be 'genome' or 'chromosome'", call.=FALSE)
  }

  # check style of the Y axis argument
  if(args[2] %in% c("bp", "percentage")){
    style = args[2]
  } else {
    stop("The argument 'style' can only be 'bp' or 'percentage'", call.=FALSE)
  }

  filename = args[3]
  indexname = args[4]
  plotname = args[5]

  # collect the remaining arguments that must be the list of chromosomes to take into consideration
  if(mode == "chromosome"){
    chromosomes = c(args[5:length(args)])
  }
}

# load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

### FUNCTION
# read input files
read.landscape = function(file){
  
  require(data.table)
  data = fread(input = file)
  data = data[,c(1:17)]
  
  return(data)
}
###

data = read.landscape(file = filename)

### FUNCTION
# clean dataset
clean.landscape = function(data, mode){
  
  # keep only the columns with: name, size, divergence, ID number
  data = data[,c(5,11,8,17)]
  names(data) = c("Scaffold", "Repeat", "Size", "Divergence")

  # get "class" names from column Repeat
  # get the part of the name after # and then the part before the /
  elements = sapply(strsplit(x = data$Repeat, split = "#"), "[[", 2)
  data$Class = sapply(strsplit(x = elements, split = "/"), "[[", 1)

  # list of categories to keep
  # !!! CHANGE THIS PART TO REMOVE/KEEP CATEGORIES TO PLOT !!!
  # now I kept only the ones I thought were important
  # consider that in this present way, elements labelled as SINE? LINE? etc are not considered
  cat2keep = c("LINE", "LTR", "DNA", "Unknown", "SINE")
  data = data[data$Class %in% cat2keep,]

  # subset for the chromosomes specified in the arguments if mode = chromosome
  # otherwise delete the scaffold information and leave just "Genome" as scaffold name
  if(mode == "genome"){

    data$Scaffold = "Genome"

  } else {

    data = data[data$Scaffold %in% chromosomes,]

  }

  return(data)
}
###

data = clean.landscape(data = data, mode)

### FUNCTION
# create table with bps for each level of divergence per scaffold
abundance = function(data){

  # convert divergence and round it
  data$Divergence = data$Divergence * 100
  data$RoundDiv = floor(data$Divergence)

  ### FACTOR AND CLASS ARE THE IDENTIFIERS YOU CAN WORK ON TO USE GENERAL OR SPECIFIC CLASSES OF REPEATS
  data$Factor = paste(data$Class, data$RoundDiv, sep = "$")

  # convert table into data.frame and sum the bps at each divergence per chromosome
  # if chromosomes are present, otherwise only "genome" is taken into consideration
  data = as.data.frame(data)
  data_bps = aggregate(x = data$Size, by = data[c("Scaffold", "Factor")], FUN = sum, na.rm = TRUE)

  # restore classes (categories) and divergence from Factor column
  data_bps$Class = sapply(strsplit(data_bps$Factor, "\\$"), "[[", 1)
  data_bps$Divergence = sapply(strsplit(data_bps$Factor, "\\$"), "[[", 2)

  # discard divergence higher than 50%
  data_bps$Divergence = as.integer(data_bps$Divergence)
  boo = data_bps$Divergence <= 50
  data_bps = data_bps[boo,]

  # convert size into megabases
  data_bps$Mb = data_bps$x / 1000000

  return(data_bps)
}
###

### FUNCTION
# create table with percentage of repeats for each level of divergence per scaffold
percentage = function(data, indexname, mode){

  # read index file
  index = read.table(indexname, sep = "\t", stringsAsFactors = FALSE)

  # convert divergence and round it
  data$Divergence = data$Divergence * 100
  data$RoundDiv = floor(data$Divergence)

  ### FACTOR AND CLASS ARE THE IDENTIFIERS YOU CAN WORK ON TO USE GENERAL OR SPECIFIC CLASSES OF REPEATS
  data$Factor = paste(data$Class, data$RoundDiv, sep = "$")

  # convert table into data.frame and sum the bps at each divergence per chromosome
  # if chromosomes are present, otherwise only "genome" is taken into consideration
  data = as.data.frame(data)
  data_perc = aggregate(x = data$Size, by = data[c("Scaffold", "Factor")], FUN = sum, na.rm = TRUE)

  # restore classes (categories) and divergence from Factor column
  data_perc$Class = sapply(strsplit(data_perc$Factor, "\\$"), "[[", 1)
  data_perc$Divergence = sapply(strsplit(data_perc$Factor, "\\$"), "[[", 2)

  # discard divergence higher than 50%
  data_perc$Divergence = as.integer(data_perc$Divergence)
  boo = data_perc$Divergence <= 50
  data_perc = data_perc[boo,]

  if(mode == "chromosome"){

    # convert size into percentages for each scaffold based on index file
    o = match(x = data_perc$Scaffold, table = index$V1)
    data_perc$Size_scaffold = index$V2[o]
    data_perc$Percentage = (data_perc$x / data_perc$Size_scaffold) * 100
  
  } else {

    total_size = sum(index$V2)
    data_perc$Size_scaffold = total_size
    data_perc$Percentage = (data_perc$x / total_size) * 100

  }

  return(data_perc)
}
###


###FUNCTION
# get colors
# this function is a bit chaotic, can be simplified but I'm lazy, it works :P
get.colors = function(data_bps){
  # assign a color for each category
  # generate n colors from the Spectral palette
  # the palette can be modified!!
  # you can also do it manually with the colors you prefer of course
  palette = colorRampPalette(brewer.pal(length(unique(data_bps$Class)), "Spectral"))(length(unique(data_bps$Class)))

  # here I order the categories by alphabetic order but another custom order can be done!!!
  levels = unique(data_bps$Class)
  # create a table that contains a column with the categories of repeats and the color assigned to it in the second column
  data_col = data.frame(Repeat = levels, Colors = palette)
  data_col$Colors = as.character(data_col$Colors)

  # this is a fix because I always have problems with ggplot
  o = match(table = data_col$Repeat, x = unique(data_bps$Class))
  colors = data_col$Colors[o]
 
  return(colors)
}
###

if(style == "bp"){

  data_plot = abundance(data)
  
  colors = get.colors(data_plot)

  # save the table with abundance, just in case
  # comment it away if you do not want this output
  write.table(x = data_plot, file = sub(pattern = ".k2p.noCpG.size", replacement = "_divergence_abundance.txt", x = filename), sep ="\t", quote = F, row.names = F, col.names = T)

  # create and save plot
  plot = ggplot(data = data_plot, aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity", colour = "gray90", size = 0.15) + xlab("Divergence (%)") + ylab("Genome occupied by repeats (Mb)") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA)) + xlim(NA, 51) + scale_fill_manual(values = colors) + facet_wrap(vars(Scaffold))


} else {

  data_plot = percentage(data, indexname, mode)

  colors = get.colors(data_plot)

  # save the table with percentages, just in case
  # comment it away if you do not want this output
  write.table(x = data_plot, file = sub(pattern = ".k2p.noCpG.size", replacement = "_divergence_percentage.txt", x = filename), sep ="\t", quote = F, row.names = F, col.names = T)

  # create and save plot
  plot = ggplot(data = data_plot, aes(x = as.integer(Divergence), y = Percentage, fill = Class)) + geom_bar(stat = "identity", colour = "gray90", size = 0.15) + xlab("Divergence (%)") + ylab("Percentage of genome occupied by repeats (%)") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA)) + xlim(NA, 51) + scale_fill_manual(values = colors) + facet_wrap(vars(Scaffold))
  
}


if(mode == "genome"){

  # size of an horizontal A4
  ggsave(filename = paste0(plotname, ".pdf"), plot = plot, device = "pdf", width = 30, height = 21, units = "cm", scale = .5)

} else {

  # size of an horizontal A4
  ggsave(filename = paste0(plotname, "_chromosomes.pdf"), plot = plot, device = "pdf", width = 30, height = 21, units = "cm", scale = .5)

}