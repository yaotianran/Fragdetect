# some useful helper functions
library(tools)
library(Rsamtools)
library(GenomicAlignments)
library(stringr)


# read parameter
# input a headless two-column text file, column one is parameter, column two is value
# output is a list
read_parameters <- function(paramter_file = 'parameter.txt') {
   parameters.df = read.csv(paramter_file, sep = '\t', header = F, stringsAsFactors = F, row.names = 1, comment.char = '#', strip.white = T, blank.lines.skip = T)
   parameters.list = split(parameters.df$V2, 1:nrow(parameters.df))
   names(parameters.list) = rownames(parameters.df)
   f <- function(s) {
      if (s %in% c("T", "TRUE", "True", "true", "F", "FALSE", "False", "false")) {
         return(as.logical(s))
      } else {
         s = tryCatch(as.numeric(s),
                      error = function(e) {s},
                      warning = function(w) {s})
      }
      return(s)
   }
   parameters.list = lapply(parameters.list, FUN = f)
   print(parameters.list)
   return(parameters.list)
}


# use ms5sum to get file md5 value
md5 <- function(file) {
   command = sprintf('md5sum %s', file)
   r = system(command, intern = T)
   md5 = str_split(r, '\\s+', simplify = T)[1, 1]
   return(md5)
}

# read a bam file and return a fragment GRange
read_bam <- function(bam.file, folder = './fragments', min.mapq = 1, min.length = 10, max.length = 1000, remove.unstranded = T) {
   r = dir.create(folder, showWarnings = F)
   #basename.str = md5(bam.file)
   basename.str = str_split(basename(bam.file), '\\.', simplify = T)[1,1]
   data.file = sprintf('%s/%s.RData', folder, basename.str)

   if (file.exists(data.file)) {
      message(sprintf('Read from %s', data.file))
      r = load(data.file)
   } else {
      chromosomes = paste0('chr', 1:22)
      param <- Rsamtools::ScanBamParam(what = c("seq"), flag = scanBamFlag(isDuplicate = F, isSecondaryAlignment = F, isUnmappedQuery = F), mapqFilter = min.mapq)
      galp <- GenomicAlignments::readGAlignmentPairs(bam.file, param = param)
      frag <- granges(keepSeqlevels(galp, chromosomes, pruning.mode = "coarse"), on.discordant.seqnames = "drop") # convert aligned segments to cfDNA fragments
      frag = frag[min.length <= width(frag) & width(frag) <= max.length]
      if (remove.unstranded) {
         frag = frag[strand(frag) %in% c('+', '-')]
      }
      save(frag, file = data.file)
      rm(galp)
      r = gc()
   }
   message(sprintf('%s  ----->  %s fragments  ----->  %s', bam.file, length(frag), data.file))
   return(frag)
}


# manually change the chromosome aneuploidy by adding or extracting extra fragments
# frag is a Grange object for fragments
# chr6 35,993,818 -36,108,087

change_aneuploidy <- function(frag, range, prop) {

   message('extracting range ...')
   count.c = countOverlaps(frag, range)
   temp.gr = frag[count.c > 0]

   if (length(temp.gr) == 0) {
      message('No hit in range')
      return(NA)
   }
   rest.gr = frag[count.c == 0]

   message('slicing object ...')
   temp.df = slice_sample(data.frame(temp.gr), prop = prop, replace = T)
   temp.gr = sort(makeGRangesFromDataFrame(temp.df))
   return(c(rest.gr, temp.gr))

}


# then
write_cnvkit_bin <- function(bin, file) {
   temp.c = sapply(bin, FUN = length)
   cnv.list = bin[temp.c == 1]
   names(cnv.list) = NULL  # do.call won't work if have list name
   cnv.gr = do.call('c', cnv.list)
   cnv.gr = sort(cnv.gr)
   data.frame(cnv.gr) %>% mutate(strand = NULL, width = NULL) %>% write_tsv(file, col_names = F)
   return(cnv.gr)
}
