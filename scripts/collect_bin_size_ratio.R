# This function give the count of each fragments of each  intervals in each bin
# frag : a GenomicRange object, each row is a fragment range
# bins: a list, in which each item is a bin represented by a GenomicRange object
# breaks: a integer vector c(i1, i2, .... in). This is the size intervals you want to count
#   fragments size in interval [i1, i2), [i2, i3), .... [i(n-1), in) will be counted and included in output
#   fragments with size < i1 and >= in will be discarded
#   for esample if breaks is c(100, 200, 300), then output will include the count of fragments with size of 100~199bp and 200~299bp
#   fragments < 100bp and >= 300bp will not included.
# core : how many CPU core you want to use

# the output is N x (n + 1) matrix , N is bin length, n is interval count
# 1st column (float) is the mean GC of total fragments (all fragments included) in this bin
# 2nd ~ nth columns (integer) are the count of fragment in each size interval in this bin
# (n + 1)th column (integer) is total fragment (all fragments included) in this bin
collect_bin_size_ratio <- function(frag, bins, breaks, core = 2) {
   if (!'GC' %in% colnames(frag@elementMetadata)) {
      message('fragments object has no GC info')
      quit(save = 'no', status = -1)
   }

   # The Fragment Size Ratio or FSR, similar to the original fragmentation profile reported by DELFI [7], was
   # generated using the short/intermediate/long fragments ratios except using different cutoffs:
   # the short, intermediate and long fragments were defined as 65-150bp, 151-220bp and 221-400bp,
   # according to the overall fragment lengths profile in our cohorts.
   # get the following features for each bin: average fragment GC ratio,
   ..get.size.ratio <- function(bin, fragments, breaks) {
      query = findOverlaps(fragments, bin)
      if (length(queryHits(query)) == 0) {  # if no hits
         result = c(gc = NA, short = 0, medium = 0, long = 0, n = 0)
         return(result)
      }
      fragments.hits = fragments[unique(queryHits(query))]  # collect the fragments in the bin
      gc.total = mean(fragments.hits$GC)
      n = length(fragments.hits)

      # when provided by N breaks, we can make N+2 intervals, divided fragments into N+1 groups by size
      size.breaks = c(0, breaks, max(width(fragments.hits)))
      h = hist(width(fragments.hits), breaks = size.breaks, plot = F, right = F)
      count.c = h$counts[2:(length(h$counts) - 1)]  # we discard the shortest and longest fragments

      result = c(gc.total, count.c, n)
      return(result)
   }

   temp.list = bettermc::mclapply(bins, FUN = ..get.size.ratio, frag, breaks = breaks, mc.cores = core)
   bin.ma = do.call(rbind, temp.list)
   return(bin.ma)
}



