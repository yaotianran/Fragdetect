



collect_fragment_motif <- function(frag, genome.reference, upstream, downstream) {
   temp.ma = gtools::permutations(4, (parameters.list$motif_upstream + parameters.list$motif_downstream), c('A', 'T', 'C', 'G'), repeats.allowed = T)
   motif.c = apply(temp.ma, MARGIN = 1, FUN = paste0, collapse = '')

   temp.count = sample(1:99999, size = length(motif.c), replace = T)
   motif.df = data.frame(X1 = motif.c, X2 = sample(temp.count, size = length(motif.c)), X3 = sample(temp.count, size = length(motif.c)))
   return(motif.df)
}

collect_fragment_GC <- function(frag, genome.reference) {
   temp = frag
   gc.c = runif(length(frag))
   temp$GC = gc.c
   return(temp)
}

collect_bin_sizeratio <- function(frag, bins, mc.cores) {

   gc.c = runif(length(bins))

   temp.count = sample(1:99999, size = length(bins), replace = T)

   count.ma = matrix(nrow = length(bins), ncol = 5)
   count.ma[, 1] = gc.c

   temp1.count = sample(1:99999, size = length(bins), replace = T)
   count.ma[, 2] = temp1.count

   temp2.count = sample(1:99999, size = length(bins), replace = T)
   count.ma[, 3] = temp2.count

   temp3.count = sample(1:99999, size = length(bins), replace = T)
   count.ma[, 4] = temp3.count

   count.ma[, 5] = temp1.count + temp2.count + temp3.count

   return(count.ma)

}



GC_correct <- function(train, test) {

   result.df = data.frame(matrix(nrow = nrow(test), ncol = 3))

   temp1.count = sample(1:99999, size = nrow(test), replace = T)
   temp2.count = sample(1:99999, size = nrow(test), replace = T)
   temp3.count = sample(1:99999, size = nrow(test), replace = T)
   temp.c = temp1.count + temp3.count + temp3.count

   result.df = data.frame(X1 = temp1.count/temp.c,
                          X2 = temp2.count/temp.c,
                          X3 = temp3.count/temp.c)

   return(result.df)
}


collect_bin_distribution <- function(frag, bins, mc.cores) {

   result.df = data.frame(matrix(nrow = length(bins), ncol = 67))

   temp.ma = matrix(sample(1:99999, size = length(bins) * 67, replace = T), ncol = 67)
   temp.c = t(t(temp.ma) / colSums(temp.ma))

   return(data.frame(temp.c))
}





