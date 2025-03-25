#!/usr/bin/env Rscript
#Script to generate CNVkit plot
#Brando Poggiali 05/12/2022

#EXAMPLE OF COMMAND LINE: 
#Rscript --vanilla CNVKIT_graph.R ASL_ONC_ONC00003TTD01R.cnr ASL_ONC_ONC00003TTD01R.cns ASL_ONC_ONC00003TTD01R_modified.cnr ASL_ONC_ONC00004BLD01R_filtered_log2_genes.txt centromere_position.tsv
args=commandArgs(trailingOnly=TRUE)
file_cnr=args[1]
file_cns=args[2]
brando_cnr=args[3]
genes_file=args[4]
centromere_position=args[5]

###Download/Upload packages--------------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(stringr)

###Upload data---------------------------------------------------------------------------------------------------
#genes <- as.character(read.csv("genes.txt",header=FALSE))

cnr_data <- fread(file_cnr,sep="\t")
cns_data<- fread(file_cns,sep="\t")
brando_cnr <- fread(brando_cnr,sep="\t")
dup_del_genes<- fread(genes_file,sep="\t")
centromere_positions <- fread(centromere_position,sep="\t")
colnames(centromere_positions) <-c("CHROM", "START", "END") 

sample_name=str_remove(file_cns, ".cns")

###Filter dataset-------------------------------------------------------------------------------------------------
#set variables
plot_point=(cnr_data$end[1]-cnr_data$start[1])/2
autosomal <- paste0("chr",c(1:22,"X","Y"))

#For loop to plot all chromosomes
for (chr in autosomal){
  cnr_chr <- cnr_data[cnr_data$chromosome==chr]
  cns_chr <- cns_data[cns_data$chromosome==chr]
  brando_cnr_chr <- brando_cnr[brando_cnr$chromosome==chr]
  dup_del_genes_chr <- dup_del_genes[dup_del_genes$chromosome==chr]
  centromere_position_chr <- centromere_positions[centromere_positions$CHROM==chr]
  genes=paste(dup_del_genes_chr$gene, sep=",", collapse=",")
  genes=str_split(genes, ",")[[1]]
  ###Plots------------------------------------------------------------------------------------------------------------


  graph <- ggplot() +
    geom_point(data=cnr_chr, aes(x=(start+plot_point)/1000000, y=log2), color="#616a6b",size=0.4) + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.minor.y = element_line(size=0.1, colour="grey50", linetype="dashed"),
    panel.grid.major.y = element_line(size=0.1, colour = "grey50",linetype="dashed")) +
    ggtitle(paste("Chromosome", str_remove(chr,"chr"))) +xlab("Position (Mb)") + ylim(-5,2) +
    theme(plot.title = element_text(hjust = 0.5,size = 22)) +
    geom_hline(yintercept = 0,size=1.08 ) +  #set black line at y=0
    scale_x_continuous( expand = c(0.01, 0.01)) +  #set side layout graph
    geom_rect(data=as.data.frame(centromere_position_chr), inherit.aes=FALSE, aes(xmin=START/1000000, xmax=END/1000000, ymin=-Inf,ymax=Inf ), color="transparent", fill="orange", alpha=0.27) +  #Plot centromere position
    geom_segment(data=cns_chr, aes(x=start/1000000, xend=end/1000000, y=log2, yend=log2),lineend="round"  ,color="deepskyblue",size=1.45) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 11), axis.title = element_text(size = 15)) 
  

  if (!genes[1]=="") {
  position_genes <- subset(brando_cnr_chr, gene %in% genes)
  unique_position_genes <- position_genes[!duplicated(position_genes$gene)]
  
  a <- seq(1.4,3,by=0.20)
  b <- seq(3.2,5,by=0.20)
  x <- b
  value_to_change <- "b"
  
  for (row in 1:nrow(unique_position_genes)) {
    gene_location=unique_position_genes[row,start]/1000000
    gene_name= unique_position_genes[row,gene]
    
    
    #if vectors are empty fill in again
    if(length(x)==1 && value_to_change=="a"){
      a <- seq(1.4,3,by=0.20)
      a <- a[! a %in% last_value_a]
      x <- a
    }
    
    if(length(x)==1 && value_to_change=="b"){graph
      b <- seq(3.2,5,by=0.20)
      b <- b[! b %in% last_value_b]
      x <- b
    }
    
    y_value=sample(x,1)
    
    #remove value from vectors
    if (value_to_change=="a"){
      a=a[! a %in% y_value]
      last_value_a <- y_value
    }
    if (value_to_change=="b"){
      b=b[! b %in% y_value]
      last_value_b <- y_value
    }
    
    
    
    #change variable to choose the vector to use 
    if (value_to_change=="a") {
      value_to_change <- "b"
      x <- b
    } else {
      value_to_change <- "a"
      x <- a
    }
    
    
    graph <- graph + geom_vline(xintercept=gene_location, color="#f9a4a4", size=0.15) +
      annotate("text",x=gene_location, y=-y_value,label = gene_name, col="red", size=2.80)
    
  }}
  
  ggsave(filename=paste0(sample_name,"_", chr,".png"), width = 15, height = 7)   
  #ggsave(filename=paste0(sample_name,"_", chr,".svg"), width = 15, height = 7) #if you want to save in svg
}



#, panel.grid.major.y = element_line(size=0.15,color = "grey50")











