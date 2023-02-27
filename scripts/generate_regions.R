library(stringr)
library(tidyverse)
library(GenomicRanges)

#mapping.file = '~/sda1/projects/fragmentomics/data/mapping_b37.interval'
#mapping.file = '~/sda1/projects/fragmentomics/data/mapping_hg19.interval'
#mapping.file = '~/sda1/projects/fragmentomics/data/mapping_hg38.interval'

unmappable.file = '~/sda1/projects/fragmentomics/data/unmappable_b37.interval'
unmappable.file = '~/sda1/projects/fragmentomics/data/unmappable_hg19.interval'
unmappable.file = '~/sda1/projects/fragmentomics/data/unmappable_hg38.interval'


gap.file = '~/sda1/projects/fragmentomics/data/gap_b37.tsv'
gap.file = '~/sda1/projects/fragmentomics/data/gap_hg19.tsv'
gap.file = '~/sda1/projects/fragmentomics/data/gap_hg38.tsv'


cytoband.file = '~/sda1/projects/fragmentomics/data/cytoBandIdeo_b37_clean.csv'
cytoband.file = '~/sda1/projects/fragmentomics/data/cytoBandIdeo_hg19_clean.csv'
cytoband.file = '~/sda1/projects/fragmentomics/data/cytoBandIdeo_hg38_clean.csv'

hic.compartments.file = '~/sda1/projects/fragmentomics/data/hic_compartments_100kb_ebv_2014_b37.txt'
hic.compartments.file = '~/sda1/projects/fragmentomics/data/hic_compartments_100kb_ebv_2014.txt'
hic.compartments.file = '~/sda1/projects/fragmentomics/data/hic_compartments_100kb_ebv_2014_hg38_LiftOver.bed'


chromosomes = as.character(1:22)
chromosomes = paste0('chr', 1:22)


# read mapping regions
# mapping_file is a interval file generated from 'gatk ScatterIntervalsByNs --OUTPUT_TYPE ACGT'
read_interval_list <- function(interval_file, chromosomes = paste0('chr', 1:22), min_width = 300) {
   mapping.df = read_tsv(interval_file, comment = '@', col_names = F, show_col_types = F)
   mapping.df %<>% dplyr::filter(X1 %in% chromosomes) %>% mutate(X3 = X3)
   mapping.gr = makeGRangesFromDataFrame(mapping.df, seqnames.field = 'X1', start.field = 'X2', end.field = 'X3', starts.in.df.are.0based = F)
   mapping.gr = reduce(mapping.gr)
   mapping.gr = mapping.gr[width(mapping.gr) >= min_width]
   print(mapping.gr)
   return(mapping.gr)
}


# read the gap file
read_gap_regions <- function(gap_file, types = c('centromere', 'telomere'), chromosomes = paste0('chr', 1:22)) {
   gap.df = read_tsv(gap_file, comment = '#', col_names = T, show_col_types = F)
   gap.df %<>% dplyr::filter(chrom %in% chromosomes, type %in% types) %>% mutate(chromEnd = chromEnd)
   gap.gr = makeGRangesFromDataFrame(gap.df, seqnames.field = 'chrom', start.field = 'chromStart', end.field = 'chromEnd', starts.in.df.are.0based = T)
   gap.gr = sort(gap.gr)
   print(gap.gr)
   return(gap.gr)
}


# read cytoband file
read_cytoband <- function(cytoband_file, chromosomes = paste0('chr', 1:22)) {
   arms.df = read_tsv(cytoband_file, col_names = F, show_col_types = F)
   arms.df %<>% dplyr::filter(X1 %in% chromosomes)
   arms.df %<>% mutate(X6 = str_split(X4, '\\d+', simplify = T)[, 1])
   temp.list = split(arms.df, f = ~ arms.df$X1 + arms.df$X6)

   f <- function(tb) {  # tb is a tibble, reduce arm region
      g = reduce(makeGRangesFromDataFrame(tb, ignore.strand = T, seqnames.field = 'X1', start.field = 'X2', end.field = 'X3', starts.in.df.are.0based = T, keep.extra.columns = T))
      if (length(unique(tb$X6)) == 1) {  # if all regions in one arm
         g$arm = unique(tb$X6) # reduce function will remove arm column, so have to add it manually
      }
      return(g)
   }
   temp.list = lapply(temp.list, FUN = f)
   names(temp.list) = NULL  # do.call won't work if have list name
   arms.gr = do.call('c', temp.list)
   seqlevels(arms.gr) = chromosomes
   arms.gr = arms.gr[order(seqnames(arms.gr), arms.gr$arm)]
   print(arms.gr)
   return(arms.gr)
}


unmappable.gr = read_interval_list(unmappable.file) # 260
gap.gr = read_gap_regions(gap.file)   # 64
discard.gr = reduce(c(unmappable.gr, gap.gr)) # 260
arms.gr = read_cytoband(cytoband.file)  # 44


regions = list()
ab.df = read_tsv(hic.compartments.file, col_names = T, show_col_types = F)
ab.df %<>% dplyr::filter(chr %in% chromosomes) %>% mutate(end = end + 1)
ab.gr = makeGRangesFromDataFrame(ab.df, starts.in.df.are.0based = T)  # 26421
ab.gr <- ab.gr[-queryHits(findOverlaps(ab.gr, discard.gr))]   # 26237, remove the 100K hic compartments regions that overlaps the discard regions

# # # # # # # # # # # # # # # generate the 100K and 5M hic compartments regions for DELFI, 500K and 1M regios for CNA/aneuploidy detection # # # # # # # # # # # # # # #
# assign chromosome arms to each hic compartments region
temp = findOverlaps(ab.gr, arms.gr)  # find the overlaps with chromsome arms regions
ab.gr <- ab.gr[temp@from]   # sort the regios by arms
ab.gr$arm <- arms.gr$arm[temp@to]
ab.100K.list = split(ab.gr, as.factor(ab.gr))

temp = as_tibble(ab.gr) %>% group_by(seqnames, arm) %>% mutate(seqnames = as.character(seqnames), name = str_c(seqnames, '.', arm, '.', 1:length(arm)))
names(ab.100K.list) = temp$name
regions$AB.100K = ab.100K.list  # 26237

# try to combine each 5 100K regions into a 500K region
ab.df = data.frame(ab.gr) %>% select(-strand)
ab.df %<>% group_by(seqnames, arm) %>% mutate(combine = ifelse(arm == 'p', yes = ceiling((1:length(arm))/5), no = ceiling(rev((1:length(arm))/5))))
temp.list = split(ab.df, f = data.frame(ab.df$seqnames, ab.df$arm, ab.df$combine), drop = T)
temp.c = sapply(temp.list, FUN = nrow)
temp.c = names(temp.c)[temp.c == 5]
ab.500K.list = temp.list[temp.c] # only keep the big regions with 50 100K regions so that the big region will be 5M after concated
ab.500K.list = lapply(ab.500K.list, FUN = function(tb) {tb %>% makeGRangesFromDataFrame() %>% reduce()})
regions$AB.500K = ab.500K.list  # 5233

# try to combine each 10 100K regions into a 1M region
ab.df = data.frame(ab.gr) %>% select(-strand)
ab.df %<>% group_by(seqnames, arm) %>% mutate(combine = ifelse(arm == 'p', yes = ceiling((1:length(arm))/10), no = ceiling(rev((1:length(arm))/10))))
temp.list = split(ab.df, f = data.frame(ab.df$seqnames, ab.df$arm, ab.df$combine), drop = T)
temp.c = sapply(temp.list, FUN = nrow)
temp.c = names(temp.c)[temp.c == 10]
ab.1M.list = temp.list[temp.c] # only keep the big regions with 50 100K regions so that the big region will be 5M after concated
ab.1M.list = lapply(ab.1M.list, FUN = function(tb) {tb %>% makeGRangesFromDataFrame() %>% reduce()})
regions$AB.1M = ab.1M.list  # 2609

# try to combine each 50 100K regions into a 5M region
ab.df = data.frame(ab.gr) %>% select(-strand)
ab.df %<>% group_by(seqnames, arm) %>% mutate(combine = ifelse(arm == 'p', yes = ceiling((1:length(arm))/50), no = ceiling(rev((1:length(arm))/50))))
temp.list = split(ab.df, f = data.frame(ab.df$seqnames, ab.df$arm, ab.df$combine), drop = T)
temp.c = sapply(temp.list, FUN = nrow)
temp.c = names(temp.c)[temp.c == 50]
ab.5M.list = temp.list[temp.c] # only keep the big regions with 50 100K regions so that the big region will be 5M after concated
ab.5M.list = lapply(ab.5M.list, FUN = function(tb) {tb %>% makeGRangesFromDataFrame() %>% reduce()})
regions$AB.5M = ab.5M.list  # 504

# # # # # # # # # # # # # # # generate the arms regions for count fragment size distributions # # # # # # # # # # # # # # #
# combine according to chromosome arms
ab.df = data.frame(ab.gr) %>% select(-strand)
temp.list = split(ab.df, f = data.frame(ab.df$seqnames, ab.df$arm), drop = T)
f <- function(df) {  # tb is a tibble
   g = reduce(makeGRangesFromDataFrame(df))
   if (length(unique(df$arm)) == 1) {
      g$arm = unique(df$arm) # reduce function will remove arm column, so have to add it manually
   }
   return(g)
}
ab.arms.list = lapply(temp.list, FUN = f)
regions$AB.arm = ab.arms.list  # 40


save(regions, file = 'data/fragmentomics_regions.RData')


