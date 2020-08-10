# Visualize the results of SW's break point model with circlize
library(circlize)
library(tidyverse)


get_SV <- function(fname){
  svs = read_tsv(fname)
  svs %>% rowwise() %>%  mutate(SVTYPE= strsplit(extra, "=")[[1]][2]) %>% as.data.frame() -> svs

  # sv.first and sv.second are dataframes with columns: "chr","start","end"
  sv.first <- svs %>% select(chr=chr1, start=coord1) %>% as.data.frame()
  sv.first$chr = paste0("chr", sv.first$chr)
  sv.first$end = sv.first$start

  sv.second <- svs %>% select(chr=chr2, start=coord2) %>% as.data.frame()
  sv.second$chr = paste0("chr", sv.second$chr)
  sv.second$end = sv.second$start
  
  return(list(svs = svs, sv.first=sv.first, sv.second=sv.second))
}


plot_SV <- function(fname){
  res = get_SV(fname)
  svs = res$svs
  sv.first = res$sv.first
  sv.second = res$sv.second

  pa <- c("#377eb8","#d95f02","#33a02c")
  names(pa) <- c("DEL","INS","INV")
  col <- pa[svs$SVTYPE]
  
  
  for(i in 0:(tail(svs,n=1)[1,"cycleID"])) {
    
    fname <- paste0("sv_data_", i, ".tsv")
    outfile <- str_replace(fname, ".tsv", ".png")
    name <- ""
    png(outfile, height=600, width=600)
    
    print(sv.first)
    circos.initializeWithIdeogram(species="hg38")
    circos.genomicLink(sv.first[svs$cycleID==i,], sv.second[svs$cycleID==i,], col=col[svs$cycleID==i])
    
    circos.clear()
    
    title(name)
    dev.off()
  }
}


# Too many segments to show in the circos plots
get_CN <- function(fname, cutoff = 5*10e4){
  cns = read_tsv(fname)

  cns$chr = paste0("chr", cns$chr)
  
  cns %>% select(chr, start, end, cn) -> cnbed
  unique(cns$cn)
  summary(cnbed$end-cnbed$start)
  
  #cnbed %>% filter(value1 + value2 != 1)  %>% filter(end-start > cutoff) -> cnsel
  cnbed %>% filter(end-start > cutoff) -> cnsel
  str(cnsel)
  
  unique(cnsel$cn)
  
  return(cnsel)
}


plot_CN_SV <- function(fcnv, fsv, fout, cutoff = 5*10e4){
  cnsel = get_CN(fcnv)
  
# diploid genome:
#  cols = c("#6283A9","#bdd7e7","#f0f0f0","#FCAE91", "#B9574E", "#76000D", "#3b0107")
#  col_fun = colorRamp2(breaks = seq(0:6)-1, colors = cols)
 
# haploid genome: 
  cols = c("#6283A9","#f0f0f0", "#B9574E", "#3b0107")
  col_fun = colorRamp2(breaks = seq(0:3)-1, colors = cols)
  
  res = get_SV(fsv)
  svs = res$svs
  sv.first = res$sv.first
  sv.second = res$sv.second
  
  pa <- c("#377eb8","#d95f02","#33a02c")
  names(pa) <- c("DEL","INS","INV")
  col <- pa[svs$SVTYPE]
  
  name <- ""
  png(fout, height=600, width=600)
  
  circos.initializeWithIdeogram(species="hg38")
  circos.genomicHeatmap(cnsel, col = col_fun, side = "inside", border = "white")
  circos.genomicLink(sv.first, sv.second, col=col)
  circos.clear()
  
  title(name)
  dev.off()
}

setwd("/Users/samuelwinnall/Documents/UCL/Thesis/Modelling/Current/output/00")
fsv = "sv_data.tsv"
plot_SV(fsv)

fcnv = "cn_data.tsv"
cutoff = 5*10e6
fout = "sv_cnv.png"
plot_CN_SV(fcnv, fsv, fout, cutoff)

