
# This function extract the fragment motif and calculate the count
# frag
# the output is a three-column headless data frame (X1, X2, X3) , each row is a motif statistics
# column one is motif sequence (char),
# column two is motif count at 5-end (integer)
# column three is motif count at 3-end (integer)

# temp/cnv/hic.bed
collect_bin_aneuploidy <- function(bam.file, genome.reference, core, execute = 'cnvkit.py') {
   r = NA

   merge_cnr_cns <- function(cnr.file, cns.file) {
      cnr.df = read_tsv(cnr.file, show_col_types = F)
      cnr.gr = makeGRangesFromDataFrame(cnr.df, keep.extra.columns = T)

      cns.df = read_tsv(cns.file, show_col_types = F)
      cns.gr = makeGRangesFromDataFrame(cns.df, keep.extra.columns = T)

      query.df = data.frame(findOverlaps(cns.gr, cnr.gr))
      log2.df = data.frame(matrix(nrow = 0, ncol = 3))
      for (i in 1:nrow(cnr.df)) {
         #message(i)
         if (i %in% query.df$subjectHits) {
            index = match(i, query.df$subjectHits)
            queryHits = query.df$queryHits[index]
            if ('cn' %in% colnames(cns.df)) {
               temp.df = cns.df[queryHits, c('log2', 'depth', 'weight', 'cn')]
            } else {
               temp.df = cns.df[queryHits, c('log2', 'depth', 'weight')]
            }

            log2.df = rbind(log2.df, temp.df)
         } else {
            log2.df = rbind(log2.df, NA)
         }
      }
      colnames(log2.df) = c('log2.seg', 'depth.seg', 'weight.seg', 'cn')
      result.df = cbind(cnr.df, log2.df)
      return(result.df)
   }


   # count
   command.str = sprintf('%s coverage -f %s -p %s -o temp/coverage.cnn %s data/cnv.bed', execute, genome.reference, core, bam.file)
   message(command.str)
   r = system(command.str, intern = T)

   # correct coverage
   command.str = sprintf('%s fix -o temp/correct.cnr temp/coverage.cnn data/empty.antitargetcoverage.for.WGS.coverage.fix.cnn data/reference_hg19.cnn', execute)
   message(command.str)
   r = system(command.str, intern = T)

   # segment
   command.str = sprintf('%s segment -m hmm-tumor -o temp/segment.cns temp/correct.cnr', execute)
   message(command.str)
   r = system(command.str, intern = T)

   cnr.file = 'temp/correct.cnr'
   cns.file = 'temp/segment.cns'
   cnv.df = merge_cnr_cns(cnr.file, cns.file)
   file.remove('temp/coverage.cnn', 'temp/correct.cnr', 'temp/segment.cns')
   return(cnv.df$log2.seg)

}
