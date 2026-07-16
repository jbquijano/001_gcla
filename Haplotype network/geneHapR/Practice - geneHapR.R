library("geneHapR")
library("ape")

seq <- 
  import_seqs("GCCOISeqs_AlignedAllAddKeep_20251022_clipkitsmartgap.fa")


hap <- 
  seqs2hap(seq)

hapSummary <- hap_summary(hap)

tc <- 
  seq[grepl("^TC[0-9]+$", names(seq))]

tchap <- 
  seqs2hap(tc)

tchapSummary <- hap_summary(tchap)

plotHapTable(tchapSummary)

tchapNet <- get_hapNet(tchapSummary,
                     groupName = "Type")

plotHapNet(tchapNet,
           size = "freq",        # circle size
           scale = "log2",       # scale circle with 'log10(size + 1)'
           cex = 0.8,            # size of hap symbol
           col.link = 2,         # link colors
           link.width = 2,       # link widths
           show.mutation = 2,    # mutation types one of c(0,1,2,3)
           legend = c(-12.5, 7)) # legend position

hapDistribution(tchap,
                database = "world",
                AccINFO = AccINFO,
                LON.col = "longitude",
                LAT.col = "latitude", 
                hapNames = c("H001", "H002", "H003"), 
                legend = TRUE)
