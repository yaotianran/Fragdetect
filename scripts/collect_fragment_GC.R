
# This function add an extra 'GC' column to frag GenomicRange object

collect_fragment_GC <- function(frag, genome.reference) {
   r = NA
   temp.df = data.frame(frag) %>% dplyr::filter(strand %in% c('+', '-'))
   temp.df %>% select(seqnames, start, end, strand) %>% write_tsv('temp/pos.tsv', col_names = F)
   command.str = sprintf('python scripts/get_GC.py %s', genome.reference)
   r = system(command.str, intern = T)

   gc.df = read_tsv('temp/gc.tsv', show_col_types = F, col_names = F)
   frag$GC = gc.df$X1
   return(frag)

}
