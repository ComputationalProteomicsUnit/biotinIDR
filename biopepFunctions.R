#---------------------------------------------------------------------------
# Author 	      : Manasa Ramakrishna, mr325@le.ac.uk
# Date started 	: 17th May, 2018
# Last modified : 18th July, 2018
# Aim 		      : Functions used for analysing disordered regions
#--------------------------------------------------------------------------- 
library(mygene)
library(data.table)
library(plyr)
library(stringr)
library(biomaRt)
library(gplots)
library(limma)
library(VennDiagram)
library(gridBase)
library(grid)
library(lattice)
library(gtools)
library(AnnotationDbi)


# -----------------------------------------------------------------------------------------------------
# Function  : 'plotTukeyHSD' 
# Aim       : To draw effect size plots with some modifications to the original 'plotTukeyHSD'
# Input     :
#   tukey.out : Object from a Tukey's HSD test
# Output    : A Tukey's effect size plot with red dotted line indicating cut-off for significance
# -----------------------------------------------------------------------------------------------------

plotTukeyHSD <- plotTukeysHSD <- function(tukey.out,
                                          x.axis.label = "Comparison",
                                          y.axis.label = "Effect Size",
                                          axis.adjust = 0,
                                          adjust.x.spacing = 5){
  
  tukey.out <- as.data.frame(tukey.out[[1]])
  means <- tukey.out$diff
  categories <- row.names(tukey.out)
  groups <- length(categories)
  ci.low <- tukey.out$lwr
  ci.up  <- tukey.out$upr                         
  
  n.means <- length(means)
  
  #determine where to plot points along x-axis
  x.values <- 1:n.means
  x.values <- x.values/adjust.x.spacing
  
  
  # calculate values for plotting limits            
  y.max <- max(ci.up) +                    
    max(ci.up)*axis.adjust
  y.min <- min(ci.low) - 
    max(ci.low)*axis.adjust
  
  if(groups == 2){ x.values <- c(0.25, 0.5)}
  if(groups == 3){ x.values <- c(0.25, 0.5,0.75)}
  
  x.axis.min <- min(x.values)-0.05
  x.axis.max <- max(x.values)+0.05
  
  x.limits <- c(x.axis.min,x.axis.max)
  
  #Plot means
  plot(means ~ x.values,
       xlim = x.limits,
       ylim = c(y.min,y.max),
       xaxt = "n",
       xlab = "",
       ylab = "",
       cex = 1.25,
       pch = 16)
  
  axis(side = 1, 
       at = x.values,
       labels = categories,
  )
  
  #Plot upper error bar 
  lwd. <- 2
  arrows(y0 = means,
         x0 = x.values,
         y1 = ci.up,
         x1 = x.values,
         length = 0,
         lwd = lwd.)
  
  #Plot lower error bar
  arrows(y0 = means,
         x0 = x.values,
         y1 = ci.low,
         x1 = x.values,
         length = 0,
         lwd = lwd.) 
  
  #add reference line at 0
  abline(h = 0, col = 2, lwd = 2, lty =2)
  
  #mtext(text = x.axis.label,side = 1,line = 1.75)
  #mtext(text = y.axis.label,side = 2,line = 1.95)
  #mtext(text = "Error bars = 95% CI",side = 3,line = 0,adj = 0)
}


#------------------------------------------------------------------------------
# Function  : mapPathwayToName (thanks to https://biobeat.wordpress.com/tag/kegg/)
# Aim       : Function to extract pathway description from KEGG given a kegg ID
# Input     : 
#   organism: name of the organism for which you want KEGG pathways
# Output    : object containing KEGG descriptions with IDs as rownames
#------------------------------------------------------------------------------
mapPathwayToName <- function(organism) {
  KEGG_PATHWAY_LIST_BASE <- "http://rest.kegg.jp/list/pathway/"
  pathway_list_REST_url <- paste(KEGG_PATHWAY_LIST_BASE, organism, sep="")
  
  pathway_id_name <- data.frame()
  
  for (line in readLines(pathway_list_REST_url)) {
    tmp <- strsplit(line, "\t")[[1]]
    pathway_id <- strsplit(tmp[1], organism)[[1]][2]
    pathway_name <- tmp[2]
    pathway_name <- strsplit(pathway_name, "\\s+-\\s+")[[1]][1]
    pathway_id_name[pathway_id, 1] = pathway_name
    
  }
  
  names(pathway_id_name) <- "pathway_name"
  return(pathway_id_name)
}

# -----------------------------------------------------------------------------------------------------------
# Function  : myProtMapper 
# Aim		: To use the function 'queryMany' from Bioconductor package mygene as fast and most up-to-date
# Input 
# 	    : ids = a character list of ids which can be uniprot, ensembl gene, gene symbol,etc
#       : id.type = what type of ids have you provided in the 'ids' list. Default = "uniprot"
#       : outlist = list of ids you want as an output. Default = c("interpro","ensembl.gene","go")
#       : modify = Logical, Default = T; Would you like to modify fields such as interpro, enseml, go to make #         them more human readable.
# Output: A dataframe with required ids and input ids 
# -----------------------------------------------------------------------------------------------------------

myProtMapper <- function(ids,id.type="uniprot",out.fields=c("interpro.short_desc","ensembl.gene","go.MF.id","go.CC.id","go.BP.id","pathway.kegg.id"),species=9606,modify=T){
  
  # Get the mapping
  qm = queryMany(ids,scopes=id.type,fields=out.fields,species=species)
  
  # Returning variable
  ret.qm = NULL
  
  # Resolve the mappings to make them human readable
  if(modify == T){
    qm$go.all = NULL
    
    # Interpro mappings
    if(is.element("interpro",colnames(qm))){
      qm$domains = sapply(qm$interpro,function(x) paste(unlist(x),collapse=";"))
    }
    else{
      print("No Interpro domains")
    }
    
    # GO mappings
    if(!is.na(grep("go",colnames(qm)))){
      
      # Grep all the go columns 'go.CC','go.MF','go.BP'
      f = grep("go",colnames(qm), value=T)
      
      qm$go.all = apply(qm[,f], MARGIN=1, FUN = function(x) paste0(as.character(unique(unlist(x))), collapse=";"))
      qm$go.all = gsub("^;","",gsub(";;",";",qm$go.all))
      qm$go.count = lengths(strsplit(qm$go.all,";"))
    }
    else{
      print("No GO terms")
    }
    
    # KEGG mappings
    if(is.element("pathway.kegg",colnames(qm))){
      qm$kegg = sapply(qm$pathway.kegg,function(x) paste(unlist(x),collapse=";"))
    }
    else{
      print("No KEGG pathways")
    }
    
    # Ensembl.gene mappings
    if(is.element("ensembl",colnames(qm))){
      qm$ens = sapply(qm$ensembl,function(x) paste(unlist(x),collapse=";"))
    }
    else{
      print("No Ensembl Ids")
    }
    
    # Return mapped structure with tidy columns
    ret.qm = qm
  }
  else{
    ret.qm = qm
  }
  
  return(data.frame(ret.qm))
}

# -----------------------------------------------------------------------------------------------------
# Function  : 'makeGene2Cat' to produce a 1:1 mapping of uniprot/ensembl/symbols to GO/Interpro terms. 
#              Will be used as input into the 'goseq' function in the gene2cat slot
# Input 
#           : dat = dataframe with ids and go/interpro terms (obtained from myProtMapper)
#           : from.id =  ids you want to map 'from'. Default = "uniprot"
#           : to.id =  ids you want to map to c("interpro","ensembl.gene","go")
#           : splt = symbol you want to split by if there are multiple ids
# Output  : A two column dataframe with Uniprot ids in the first and Go/Interpro in the second
# ------------------------------------------------------------------------------------------------------

makeGene2Cat <- function(dat,from.id,to.id,splt){
  
  cat.frame = dat[,c(from.id,to.id)]
  d.dt = data.table(cat.frame,key=colnames(cat.frame[,from.id]))
  cat.out = data.frame(d.dt[, list(to.id = unlist(strsplit(get(to.id), splt))), by=from.id])
  cat.out = unique(cat.out)
  
  return(cat.out)
}


#--------------------------------------------------------------------------------------------
# Function: rungoseq
# Aim:  goseq analysis
# Input : genes = genelist of interest; g2c = gene to category mapping i.e univ.go or univ.pro
# b = bias data ; bh = Bonferroni p-value cutoff. 
#--------------------------------------------------------------------------------------------

rungoseq<-function(genes,g2c, univ, b, bh){
  
  all.genes = rep(0,length(univ))
  names(all.genes) = univ
  all.genes[which(names(all.genes) %in% genes)] = 1
  t = table(all.genes)
  
  if(is.null(b)){
    pwf = nullp(all.genes,bias.data = NULL,genome="hg19",id="ensGene",plot.fit = T)
  }else{
    pwf = nullp(all.genes,bias.data = b,plot.fit = T)
  }
  GO.wall = goseq(pwf,gene2cat = g2c)
  GO.wall$BH = p.adjust(GO.wall$over_represented_pvalue,method = "BH")
  GO.wall$Foreground = length(intersect(genes,g2c$query))
  GO.wall$Background = length(unique(g2c$query))
  GO.wall$obsRatio = paste(GO.wall$numDEInCat,GO.wall$Foreground, sep="/")
  GO.wall$bgRatio = paste(GO.wall$numInCat,GO.wall$Background, sep="/")
  GO.wall$expectDE = ceiling(GO.wall$Foreground*(GO.wall$numInCat/GO.wall$Background))
  GO.wall$expRatio = paste(GO.wall$expectDE,GO.wall$Foreground, sep = "/")
  GO.wall$foldEnrich = round(GO.wall$numDEInCat/GO.wall$expectDE,2)
  
  GO.enriched = GO.wall[which(GO.wall$BH <= bh),]
  GO.enriched$geneID = sapply(GO.enriched$category,function(x) paste(intersect(genes,g2c$query[grep(x,g2c$to.id)]),collapse="/"))
  GO.enriched$Count = sapply(GO.enriched$category,function(x) length(intersect(genes,g2c$query[grep(x,g2c$to.id)])))
  
  GO.enriched$neg.log10.BH = -log10(GO.enriched$BH)
  GO.enriched$neg.log10.BH[which(GO.enriched$BH == 0)] = max(GO.enriched$neg.log10.BH[is.finite(GO.enriched$neg.log10.BH)])+1
  
  return(list(GO.wall,GO.enriched))
}

#--------------------
# reverseMapping
#--------------------
reversemapping=function(map){
  tmp=unlist(map,use.names=FALSE)
  names(tmp)=rep(names(map),times=as.numeric(summary(map)[,1]))
  return(split(names(tmp),as.vector(tmp)))
}

#----------------------------------------------------------------------------------------------
# Function: enricherPlot
# Aim : Modify DOSE::plot to use colours I like for plotting results of compareCluster
# Default : Will only plot the dot size to show GeneRatio and colour to show adjusted.p.val in grey and gold
# Input : 
#   data    : object from compareCluster function
#   N       : Number of top terms to show on the plot Default = 5
#   colorBy : What numeric value do you want the colour scale based on ? 
#             Default = BH or Benjamini-Hochberg adjusted p-value
#   sizeBy  : What numeric value do you want the size of the dots based on ? Default = obsRatio
#   low.col : What colour would you like your low 'colorBy' values to be ? Default = grey
#   high.col: What colour would you like your high 'colorBy' values to be ? Default = gold
#   trunc.len: At what length do you want your GO/Interpro/KEGG terms truncated ? Default = 40
#   suf     : Suffix for output file
#   all.size: What is the size that you want your legend and label text to be ? Default = 10
#   y.size  : What is the size that you want for your y-axis labels ?
#   x.size  : What is the size that you want for your x-axis labels ?
#----------------------------------------------------------------------------------------------

enricherPlot<-function(data,suf,N=5,colorBy = "BH",sizeBy = "obsRatio",low.col="#E69F00", high.col="#999999",trunc.len=40,all.size=10,y.size=12,x.size=14){
  
  #--------------------------------------------------------------------------
  # Function : topN 
  # Aim : Picks the top "N" terms in an enrichment analysis for each cluster
  #--------------------------------------------------------------------------
  
  topN <- function(res, showCategory){
    ddply(.data = res,
          .variables = .(Cluster),
          .fun = function(df, N) {
            if (length(df$obsRatio) > N) {
              if (any(colnames(df) == "pValue")) {
                idx <- order(df$pValue, decreasing=FALSE)[1:N]
              } else {
                ## for groupGO
                idx <- order(df$Count, decreasing=T)[1:N]
              }
              return(df[idx,])
            } else {
              return(df)
            }
          },
          N=showCategory)
  }
  
  # Convert 'compareCluster' result to a data.frame
  df = data.frame(data)
  
  # 'gcsize' is the number of proteins in each dataset that could be mapped to GO/Interpro/KEGG. It is the denominator in 'GeneRatio'
  # 'size' = GeneRatio is a text field - split its elements and calculate the actual GeneRatio or proportion of genes contributing to term enrichment
  # 'tot.size' = Modified x-axis labels to contain count for each cluster
  # 'mod.desc' = Modify the length of the description of terms to be 40 characters long. Anything longer will be truncated and followed by "..."
  
  gcsize = sapply(df$obsRatio,function(x) strsplit(x,"/")[[1]][2])
  #df$size =  sapply(df$obsRatio,function(x) as.numeric(strsplit(x,"/")[[1]][1])/as.numeric(strsplit(x,"/")[[1]][2]))
  df$tot.size <- paste(as.character(df$Cluster),"\n", "(", gcsize, ")", sep="")
  df$mod.desc = as.character(df$Description)
  df$mod.desc[which(nchar(df$mod.desc)>trunc.len)] = paste(substring(df$mod.desc[which(nchar(df$mod.desc)>trunc.len)],1,trunc.len),"...",sep="")
  
  # Once you've modified the main data frame, subset it to only include the top 'N' terms for each cluster
  # Order this data.frame such that the most enriched terms are at the top of the figure
  df.sub.org = topN(df,N)
  df.sub = df[which(df$mod.desc %in% unique(df.sub.org$mod.desc)),]
  
  idx <- order(df.sub[,colorBy], decreasing = F)
  df.sub$mod.desc <- factor(df.sub$mod.desc, levels=unique(df.sub$mod.desc[idx]))
  
  # Draw the plot
  pdf(paste(outdir,paste(suf,N,"enricher-dotplot.pdf",sep="_"),sep="/"),paper="a4r",width=12,height=8)
  gp = ggplot(df.sub, aes_string(x="tot.size", y="mod.desc", size=sizeBy, color=colorBy)) + geom_point() + scale_size(breaks = c(0,2,4,8,16), limits=c(0,20))+ scale_color_gradient2(low=low.col,mid=high.col, high = "#56B4E9",midpoint = quantile(df.sub[,colorBy],0.98))+xlab("")+ylab("")+guides(size=guide_legend("Fold enrichment",order=1),color=guide_colorbar(title = "-log10(adj.p.value)", title.vjust = 1.0, label.position = "bottom", reverse = F,label.theme = element_text(angle = -90)))+theme_bw()+theme(text = element_text(size=all.size),axis.text.x=element_text(size=x.size),axis.text.y=element_text(size=y.size),legend.direction = "horizontal", legend.position = "top",legend.box = "vertical")
  print(gp)
  dev.off()
  
  return(gp)
}

