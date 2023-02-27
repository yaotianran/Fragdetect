
# This function extract the fragment motif and calculate the count
# frag
# the output is a three-column headless data frame (X1, X2, X3) , each row is a motif statistics
# column one is motif sequence (char),
# column two is motif count at 5-end (integer)
# column three is motif count at 3-end (integer)

collect_fragment_motif <- function(frag, genome.reference, upstream, downstream) {
   r = NA
   temp.df = data.frame(frag) %>% dplyr::filter(strand %in% c('+', '-'))
   temp.df %>% select(seqnames, start, end, strand) %>% write_tsv('temp/pos.tsv', col_names = F)
   command.str = sprintf('python3 scripts/get_motif.py %s %s %s', genome.reference, upstream, downstream)
   r = system(command.str, intern = T)
   #r = system2(command.str)

   motif.df = read_tsv('temp/motif.tsv', show_col_types = F, col_names = F)
   return(motif.df)

}
