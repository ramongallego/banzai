#This is an R script, that given a folder with the results of banzai,
#will create some plots and dfs with the missing and unexpected sampe barcodes found
#As it is now, you have to input here the directory but working on it to launch from banzai.sh, given a OK in banzai.params.sh
banzai.output="/Users/moncho/Current_project/Kelly_HoodCanalTidal_20170413_done_sub/output/banzai_out_20170508_1607"#Your directory here

#And we need to load some R functions, and libraries
library(gridExtra)
library(stringr)
library (ggplot2)
library (reshape2)
library(Biostrings)
source('scripts/R/remap_banzai_output.r')
source('scripts/R/string.diff.ex.r')
source('scripts/R/per_lib_with_only_green.r')
source('scripts/R/fv.r')
source('scripts/R/distance_to_known.r')
source('scripts/R/Mon_rev_com.r')
source("scripts/R/multiplot.r")
#banzai_output_function<-function(banzai.output){
#If we want to ake this a shiny app, then we could use the out folder as teh input variable

### specify otu mapfile
otu_mapfile <- paste0(banzai.output,"/all_lib/derep.map")
# read otu map
otu_map <- read.table(otu_mapfile, stringsAsFactors = FALSE, header = FALSE)
otu_map[,4:6]<-with (otu_map, colsplit(V2,";",c("library","Fwd","Rev")))
otu_map[,5:6]<-apply(otu_map[,5:6], 2, gsubseq)
otu_map_no_singletons<-subset(otu_map, V3>1)
#return(write.csv(head(otu_map_no_singletons, file=paste0(banzai.output,"trial_output.csv"))))
#}
### specify sample mapfile
sample_mapfile <- paste0(banzai.output, "/sample_trans.tmp")
# read sample map
sample_map <- read.table(sample_mapfile, stringsAsFactors = FALSE, header = FALSE)

sample_map[,4:6]<-with (sample_map, colsplit(V1,";",c("library","Fwd","Rev")))
sample_map[,5:6]<-apply(sample_map[,5:6], 2, gsubseq)
#Create a df with all tags and its Fwd and Rev
Tags_16S<-data.frame(Fwd=levels(as.factor(sample_map[,"Fwd"])),
                     Tag=1:nlevels(as.factor(sample_map[,"Fwd"])))
Tags_16S[,"Rev"]<-sample_map[,"Rev"][match(Tags_16S[,"Fwd"], sample_map[,"Fwd"])]


#Now run all functions

#otu_map <- rename_samples(otu_map, sample_map)

otu_map_no_singletons<-subset(otu_map, V3>1)
otu_map_no_singletons$library<-as.factor(otu_map_no_singletons$library)
  a<-levels(as.factor(otu_map_no_singletons$library))
   test<-lapply(a,per_lib,otu=otu_map_no_singletons, map=sample_map )
  names(test)<-a
    #The output is in a list, so we have to unlist it and work on each output
  tett2<-unlist(test, recursive = F);rm(test)
  names(tett2)
  #Plots are .plot, missing and unexpectedF and R, and a df to create a second plot
  #Plots: showd them on screen and make a pdf
    b<-grep(pattern = ".plot", x = names(tett2))
    library("gridExtra")
    heatmaps<-marrangeGrob(tett2[b],ncol = 1, nrow = 2)
    heatmaps
   t<-multiplot(tett2[b], cols=1)

   png(paste0(banzai.output,"/A.png"))
   tett2[b][1]
   dev.off()
   png(paste0(banzai.output,"/all.png"))
   heatmaps
   dev.off()
#Seems like there is something really wrong. Do a similar plot only including the intended barcodes
   #heat map with only the expected barcodes - to see chances of crosscontamination
   c<-grep(pattern = ".data_to_plot", x = names(tett2))
   c1<-gsub(pattern = ".data_to_plot",replacement = "", x = names(tett2))
    greenF<-lapply(a,function(x){sample_map[which(sample_map[,"library"]==x),"Fwd"]})
    greenR<-lapply(a,function(x){sample_map[which(sample_map[,"library"]==x),"Rev"]})
    names(greenF)=names(greenR)=a
  #Create a df for all reads, and do a unique plot for all libraries on a project
    full_df_to_plot<-NULL
    for (i in 1:length(c)){
      tett2[[c[i]]]$lib=a[i]
      to_plotA<-subset(tett2[[c[i]]],subset= Fwd %in% greenF[[i]] & Rev %in% greenR[[i]])
      full_df_to_plot<-rbind(full_df_to_plot,tett2[[c[i]]])
      #Heatmaps
      p<- ggplot(data = to_plotA, aes(x=Rev, y=Fwd, fill=V3))   
      p<-p+ geom_raster()+theme_bw()+ggtitle(paste0("Library ", LETTERS[i]))+
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.y = element_text(face="bold"),
              axis.text.x = element_text(face="bold", angle = 90))
      p+scale_fill_distiller(palette = "Spectral" , name="log # Reads", trans='log')+ggsave(filename =paste0(banzai.output,"/Library_",LETTERS[i],"heatmap_log.png") )
      p+scale_fill_distiller(palette = "Spectral" , name="# Reads")+ggsave(filename =paste0(banzai.output,"/Library_",LETTERS[i],"heatmap.png") )
    }
 #Now assign to which Tag each dup could be assigned by its Fwd read
    
    full_df_to_plot$TagF<-Tags_16S[,"Tag"][match(full_df_to_plot$Fwd,Tags_16S$Fwd)]
    
    full_df_to_plot$TagR<-Tags_16S[,"Tag"][match(full_df_to_plot$Rev,Tags_16S$Rev)]
    #Now  determine if paired or unPaired
    full_df_to_plot$Paired<-ifelse(full_df_to_plot$TagF==full_df_to_plot$TagR,"Paired","UnPaired")
    full_df_to_plot[is.na(full_df_to_plot[,"Paired"]),"Paired"]<-"Unmapped"
    full_df_to_plot$Paired<-as.factor(full_df_to_plot$Paired)
    bb<-aggregate(full_df_to_plot$V3,by=(list(full_df_to_plot$lib,full_df_to_plot$TagF,full_df_to_plot$Paired)),sum)
    names(bb)<-c("lib","TagF", "Paired", "V3")
    bb$V3.log<-log2(bb$V3)
    #And finally plot it beautifully, the best plot you've ever seen, yuge
    ggplot(data = subset(bb, !is.na(TagF)), aes(x= Paired, y=V3, fill=Paired))+geom_col()+facet_grid(lib~TagF, scales = "fixed")+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+labs(y="# Reads")+
      ggsave(filename = paste0(banzai.output,"/Pairing_by_library_FwdTags.png"))
    
    ggplot(data = subset(bb, !is.na(TagF)), aes(x= Paired, y=V3.log, fill=Paired))+geom_col()+facet_grid(lib~TagF, scales = "fixed")+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+labs(y="# Reads")+
      ggsave(filename = paste0(banzai.output,"/Pairing_by_library_FwdTags_log.png"))
    
    #And Now with the reverse 
    bb<-aggregate(full_df_to_plot$V3,by=(list(full_df_to_plot$lib,full_df_to_plot$TagR,full_df_to_plot$Paired)),sum)
    names(bb)<-c("lib","TagR", "Paired", "V3")
    bb$V3.log<-log2(bb$V3) 
    ggplot(data = subset(full_df_to_plot, !is.na(TagR)), aes(x= Paired, y=V3, fill=Paired))+geom_col()+
      facet_grid(lib~TagR, scales = "fixed")+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+labs(y="# Reads")+
      ggsave(filename = paste0(banzai.output,"/Pairing_by_library_RevTags.png"))
    
    ggplot(data = subset(full_df_to_plot, !is.na(TagR)), aes(x= Paired, y=V3.log, fill=Paired))+geom_col()+
      facet_grid(lib~TagR, scales = "fixed")+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+labs(y="# Reads")+
      ggsave(filename = paste0(banzai.output,"/Pairing_by_library_RevTags_log.png"))
    
  #Missing samples
    c<-grep(pattern = ".missing", x = names(tett2))
    c1<-do.call("rbind", tett2[c])
    c1$lib_fwd<-with(c1,interaction(library,Fwd, drop = T))
    c1$sample<-sample_map[,3][match(c1$lib_fwd,sample_map$lib_fwd)]
    c1
    write.csv(c1, file=paste0(banzai.output,"/all_missing_samples.csv"))
  #Unexpected samples - We analize each set of unexpected barcodes, on both the Fwd and rev end of the amplicon
    d<-grep(pattern = ".unexpectedF", x= names(tett2))
    d1<-do.call("rbind", tett2[d])
    e<-grep(pattern = ".unexpectedR", x= names(tett2))
    e1<-do.call("rbind", tett2[e])
    write.csv(rbind(d1,e1), file=paste0(banzai.output,"/all_unexpected_barcodes.csv"))
    p<-ggplot(data=d1, aes(x=Mapped,y=Closest_Tag, alpha=N_Reads, fill=as.factor(Diffs)))
    p<- p+geom_raster()+labs(x="")+theme_bw()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
                                                    axis.ticks.y = element_blank(),axis.text.y = element_blank())
    fwd<-p+facet_grid(Closest_Tag~lib, scales="free_y", switch = "y")+theme(strip.text.y = element_text(angle = 180))
    ggsave(fwd,filename = paste0(banzai.output,"/unexpected_Fwd_barcodes.png"), dpi=600)
    q<- ggplot(data=e1, aes(x=lib,y=Closest_Tag, alpha=N_Reads, fill=as.factor(Diffs)))
    q<- q+geom_raster()+labs(x="")+theme_bw()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
                                                    axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                                                    strip.text.y = element_text(angle = 180))
    rev<-q +facet_grid(Closest_Tag~lib, scales="free", switch = "y")
    x_l<-nlevels(as.factor(e1$Barcode))
    y_l<-nlevels(as.factor(e1$lib))
    ggsave(rev,filename = paste0(banzai.output,"/unexpected_Rev_barcodes.png"), dpi=600)
