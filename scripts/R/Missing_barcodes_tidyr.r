#  A  TIDYVERSE VERSION OF THE WHOLE THING

banzai.output="/Users/moncho/banzai_outputs/banzai_out_20171128_2303/"#Your directory here

#And we need to load some R functions, and libraries
library(tidyverse)
library(gridExtra)
library(stringr)
library (ggplot2)
library (reshape2)
library(Biostrings)
source('remap_banzai_output.r')
source('string.diff.ex.r')
source('per_lib_with_only_green.r')
source('fv.r')
source('distance_to_known.r')
source('Mon_rev_com.r')
source("multiplot.r")

#banzai_output_function<-function(banzai.output){
#If we want to ake this a shiny app, then we could use the out folder as teh input variable

### specify otu mapfile
otu_mapfile <- paste0(banzai.output,"/all_lib/derep.map")
# read otu map
otu_map <- read.table(otu_mapfile, stringsAsFactors = FALSE, header = FALSE)
otu_map_tidy<- read_table2(file = otu_mapfile,col_names = c("DUP","Sample","Count")) %>%
  separate(col = "Sample",into = c("library", "Fwd", "Rev"), sep=";") 
otu_map_tidy[,c("Fwd","Rev")]<-apply ()


otu_map[,4:6]<-with (otu_map, colsplit(V2,";",c("library","Fwd","Rev")))
otu_map[,5:6]<-apply(otu_map[,5:6], 2, function (x) {gsub("ID[0-9][A-Z]=","",x)}) #remove ID1...

otu_map_no_singletons<-subset(otu_map, V3>1)