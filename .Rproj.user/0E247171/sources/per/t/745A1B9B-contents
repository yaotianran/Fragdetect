source('scripts/collect_fragment_motif.R')
source('scripts/collect_fragment_GC.R')
source('scripts/collect_bin_size_ratio.R')
source('scripts/GC_correct.R')
source('scripts/collect_bin_aneuploidy.R')
source('scripts/model.R')
source('scripts/utils.R')

#source('scripts/pseudo_code.R')

library(tidyverse)
library(magrittr)
library(GenomicAlignments)
library(digest)
library(Rsamtools)
library(bettermc)
library(gtools)

SAMPLE = 'samples.txt'
PARAMETERS = 'parameter.txt'

samples.df = read_tsv(SAMPLE, show_col_types = F, col_names = F, comment = '#')
parameters.list = read_parameters(PARAMETERS) # utils.R

SIZE.BREAKS = as.integer(str_extract_all(parameters.list$size.breaks, '\\d+')[[1]])
DIST.BREAKS = as.integer(str_extract_all(parameters.list$arms.distribution.breaks, '\\d+')[[1]])
REFERENCE = parameters.list$reference

result.file = sprintf('train_data_%s.RData', paste0(str_extract_all(Sys.time(), '\\d+')[[1]], collapse = ''))

dir.create('results', showWarnings = F)
dir.create('plots', showWarnings = F)
dir.create('temp', showWarnings = F)
dir.create('temp/cnv', showWarnings = F)


X2 = 1 # for editor only, no use
regions = list() # for editor only, no use
r = load(parameters.list$region.file)
cnv.gr = write_cnvkit_bin(regions$AB.500K, file = 'data/cnv.bed')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # How data stored # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sample.existed.df = data.frame(matrix(nrow = 0, ncol = 3))  # which samples are really calculate ?

# rows are samples
motif.ratio.ma = matrix(0, nrow = 0, ncol = 4^(parameters.list$motif_upstream + parameters.list$motif_downstream))  # each column is a motif, the value is the motif ratio
size.ratio.ma = matrix(0, nrow = 0, ncol = (length(SIZE.BREAKS) - 1)*length(regions$AB.5M))  # columns are short/media/long size fragments count ratio for each 5M region
dist.ratio.ma = matrix(0, nrow = 0, ncol = (length(DIST.BREAKS) - 1)*length(regions$AB.arm))  # columns are 65-69bp, 70-74bp, ....., 396-400bp size fragments count ratio count ratio for each chromosome arm
aneuploidy.ma = matrix(0, nrow = 0, ncol = 1*length(cnv.gr)) # columns are CNV/CNA of each 500K region


# # # # # # # # # # # # # # # # # # # # # # # # # # # # task payload # # # # # # # # # # # # # # # # # # # # # # # # # # # #
for (i in 1:nrow(samples.df)) {
   id = samples.df$X1[i]
   bam.file = samples.df$X2[i]
   type = samples.df$X3[i]
   message(sprintf('========%s/%s  %s========', i, nrow(samples.df), id))


   # # # # # # read the bam file # # # # # #
   message(sprintf('Reading %s', bam.file))
   frag = read_bam(bam.file, folder = 'results', min.mapq = parameters.list$MAPQ)  # frag is a GenomicRange object, in which each row is a fragment range, utils.R


   # # # # # # collect motif # # # # # #
   # motif.df is a three-column headless data frame (X1, X2, X3) , each row is a motif statistics
   # column one is motif sequence (char),
   # column two is motif count at 5-end (integer)
   # column three is motif count at 3-end (integer)
   message(sprintf('Collecting motif %s', bam.file))
   data.file = sprintf('results/%s_fragments_motif.RData', id)
   if (file.exists(data.file)) {
      message(sprintf('Load data file %s', data.file))
      r = load(data.file)
   } else {
      motif.df = collect_fragment_motif(frag, genome.reference = REFERENCE, upstream = parameters.list$motif_upstream, downstream = parameters.list$motif_downstream) # return motif.df
      save(motif.df, file = data.file)
   }

   temp.ma = gtools::permutations(4, parameters.list$motif_upstream + parameters.list$motif_downstream, c('A', 'T', 'C', 'G'), repeats.allowed = T)
   temp.df = data.frame(X1 = apply(temp.ma, MARGIN = 1, FUN = paste0, collapse = ''), no.use = 0)
   merged.df = merge(temp.df, motif.df, by = 'X1', all.x = T)
   merged.df$X2[is.na(merged.df$X2)] <- 0
   merged.df %<>% mutate(ratio = X2 / sum(X2))
   motif.ratio.ma = rbind(motif.ratio.ma, merged.df$ratio)


   # # # # # # collect fragment size ratio of 5M bin # # # # # #
   # calculate the GC for every fragment
   message(sprintf('Collecting GC for %s fragments', length(frag)))
   data.file = sprintf('results/%s_fragments_GC.RData', id)
   if (file.exists(data.file)) {
      message(sprintf('Load data file %s', data.file))
      r = load(data.file)
   } else {
      frag = collect_fragment_GC(frag, genome.reference = REFERENCE) # add a column 'GC' (float)
      save(frag, file = data.file)
   }

   # collect 100K/5M bin metrics (length ratio, distribution etc.)
   # The input frag is a GenomicRange object with extra 'GC' column
   # bins: a list, in which each item is a bin represented by a GenomicRange object
   # breaks: a integer vector c(i1, i2, .... in). This is the size intervals you want to count
   #   fragments size in interval [i1, i2), [i2, i3), .... [i(n-1), in) will be counted and included in output
   #   fragments with size < i1 and >= in will be discarded
   #   for esample if breaks is c(100, 200, 300), then output will include the count of fragments with size of 100~199bp and 200~299bp
   #   fragments < 100bp and >= 300bp will not included.
   # core : how many CPU core you want to use

   # In AB.5M.ma/AB.100K.ma each row is a chromosome 5M/100K region
   # the output is N x (n + 1) matrix (AB.5M.ma/AB.100K.ma) , N is bin (chromosome 5M/100K region) length, n is interval count
   # 1st column (float) is the mean GC of total fragments (all fragments included) in this bin
   # 2nd ~ nth columns (integer) are the count of fragment in each size interval in this bin
   # (n + 1)th column (integer) is total fragment (all fragments included) in this bin
   message(sprintf('Collecting size for %s fragments', length(frag)))
   data.file = sprintf('results/%s_fragments_size_ratio.RData', id)
   if (file.exists(data.file)) {
      message(sprintf('Load data file %s', data.file))
      r = load(data.file)
   } else {
      message(sprintf('%s 5M bins ...', length(regions$AB.5M)))
      exception = tryCatch(
         {AB.5M.ma = collect_bin_size_ratio(frag, bins = regions$AB.5M, breaks = SIZE.BREAKS, core = parameters.list$core)},
            error = function(e) e,
            warning = function(w) w
            )
      if (inherits(exception, 'error')) {
         AB.5M.ma = collect_bin_size_ratio(frag, bins = regions$AB.5M, breaks = SIZE.BREAKS, core = 2)
      }

      message(sprintf('%s 100K bins ...', length(regions$AB.100K)))
      exception = tryCatch(
         {AB.100K.ma = collect_bin_size_ratio(frag, bins = regions$AB.100K, breaks = SIZE.BREAKS, core = parameters.list$core)},
            error = function(e) e,
            warning = function(w) w
            )
      if (inherits(exception, 'error')) {
         AB.100K.ma = collect_bin_size_ratio(frag, bins = regions$AB.100K, breaks = SIZE.BREAKS, core = 2)
      }

      save(AB.100K.ma, AB.5M.ma, file = data.file)
      message('Save to ', data.file)
   }
   # size.ratio.df is a N x n-1 float data frame
   # N is length if 5M regions
   # columns are short/media/long size fragments count ratio
   size.ratio.df = GC_correct(AB.100K.ma, AB.5M.ma)
   temp.ma = matrix(t(size.ratio.df), nrow = 1)
   size.ratio.ma = rbind(size.ratio.ma, temp.ma)

   # # # # # # collect fragment size distrbution of each chromosome  # # # # # #
   # dist.df is a 67-column headless data frame , each row is a chromosome arm. The value is distribution ratio
   # column 1st-67th are 65-69bp, 70-74bp, ....., 396-400bp fragments count ratio in a chromosome arm.
   message(sprintf('Collecting chromosome size distribution %s', bam.file))
   data.file = sprintf('results/%s_chromosome_dist.RData', id)
   if (file.exists(data.file)) {
      message(sprintf('Load data file %s', data.file))
      r = load(data.file)
   } else {
      dist.ma = collect_bin_size_ratio(frag, bins = regions$AB.arm, breaks = DIST.BREAKS, core = parameters.list$core) # return dist.df
      dist.ma = dist.ma[sort(rownames(dist.ma)), ]
      total.c = dist.ma[, ncol(dist.ma)]
      dist.ma = dist.ma[, c(-1, -ncol(dist.ma))]
      dist.ma = dist.ma / total.c
      save(dist.ma, file = data.file)
   }

   temp.ma = matrix(t(dist.ma), nrow = 1)
   dist.ratio.ma = rbind(dist.ratio.ma, temp.ma)


   # # # # # # collect CNA  # # # # # #
   # aneuploidy.c is a float vector whose length is equal to length of 500K regions
   message(sprintf('Collecting whole genome aneuploidy %s', bam.file))
   data.file = sprintf('results/%s_chromosome_aneuploidy.RData', id)
   if (file.exists(data.file)) {
      message(sprintf('Load data file %s', data.file))
      r = load(data.file)
   } else {
      aneuploidy.c = collect_bin_aneuploidy(bam.file, genome.reference = REFERENCE, core = parameters.list$core)
      #aneuploidy.c = rep(0, length(cnv.gr))
      save(aneuploidy.c, file = data.file)
   }
   aneuploidy.ma = rbind(aneuploidy.ma, aneuploidy.c)


   sample.existed.df = rbind(sample.existed.df, samples.df[i, ])
}

save(sample.existed.df, motif.ratio.ma, size.ratio.ma, dist.ratio.ma, aneuploidy.ma, file = result.file)








