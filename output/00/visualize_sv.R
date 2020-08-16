# Visualize the results of SVGen with circlize
library(circlize)
library(tidyverse)

setwd("/Users/samuelwinnall/Documents/Github/sv-model/output/00")

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

get_CN <- function(fname, cutoff = 5*10e4){
  cns = read_tsv(fname)

  cns$chr = paste0("chr", cns$chr)

  cns %>% select(chr, start, end, cn, haplotype, cycleNum) -> cnbed
  unique(cns$cn)
  summary(cnbed$end-cnbed$start)


  cnbed %>%  filter(end-start > cutoff) -> cnsel
  str(cnsel)

  unique(cnsel$cn)

  return(cnsel)
}


plot_SV <- function(fname){
  res = get_SV(fname)
  svs = res$svs
  sv.first = res$sv.first
  sv.second = res$sv.second

  pa <- c("#377eb8","#d95f02","#33a02c")
  names(pa) <- c("DEL","INS","INV")
  col <- pa[svs$SVTYPE]


  for(i in 0:(tail(svs,n=1)[1,"cycleNum"])) {

    fname <- paste0("sv_data_", i, ".tsv")
    outfile <- str_replace(fname, ".tsv", ".png")
    name <- ""
    png(outfile, height=600, width=600)

    print(sv.first)
    circos.initializeWithIdeogram(species="hg38")
    circos.genomicLink(sv.first[svs$cycleNum==i,], sv.second[svs$cycleNum==i,], col=col[svs$cycleNum==i])
    circos.clear()

    title(name)
    dev.off()
  }
}


plot_CN_SV <- function(fsv, fcnv){

  # get CN data
  cnsel = get_CN(fcnv)

  # CN colors
  # diploid genome - 2 cell cycles means max cn = 2
  cols = c("#6283A9","#f0f0f0", "#B9574E", "#3b0107")
  col_fun = colorRamp2(breaks = seq(0:3)-1, colors = cols)

  # get SV data
  res = get_SV(fsv)
  svs = res$svs
  sv.first = res$sv.first
  sv.second = res$sv.second

  # SV colors
  pa <- c("#377eb8","#d95f02","#33a02c")
  names(pa) <- c("DEL","INS","INV")
  col <- pa[svs$SVTYPE]


  for(i in 0:(tail(cnsel$cycleNum,1))) {

    fname <- paste0("cn_sv_data_", i, ".tsv")
    outfile <- str_replace(fname, ".tsv", ".png")
    name <- ""
    png(outfile, height=600, width=600)

    circos.initializeWithIdeogram(species="hg38")

    #print(cnsel[cnsel$cycleNum==i,])
    circos.genomicHeatmap(cnsel[cnsel$cycleNum==i,], col = col_fun, side = "inside", border = "white")
    circos.genomicLink(sv.first[svs$cycleNum==i,], sv.second[svs$cycleNum==i,], col=col[svs$cycleNum==i])
    circos.clear()

    title(name)
    dev.off()
  }
}


fsv = "sv_data.tsv"
plot_SV(fsv)

fcnv = "cn_data.tsv"
plot_CN_SV(fsv, fcnv)
