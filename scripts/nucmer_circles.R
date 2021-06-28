

install.packages(c("devtools", "ggplot2"))

library(devtools)
install_github("timflutre/rutilstimflutre")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ggbio", "GenomicRanges"))
getwd()
setwd("../concat_mrl20k_q15l15k_pepper/reports/nucmer/")                     

coords_file = "nucmer.contigs_filter2_coords.txt"                                    
                     
load_coords <- function(coords_file, perc.id) {
  library(rutilstimflutre)
  library(GenomicRanges)
                       
  coords <- loadMummer(coords_file, algo = "nucmer")                    # read in nucmer results as a GRanges obj
  coords <- coords[(elementMetadata(coords)[ , "perc.id"] >= perc.id)]  # filter entries with perc.id != 100
  seqlevels(coords) <- seqlevelsInUse(coords)                           # drop seq levels that are no longer used
  coords$qry.name <- names(coords)                                      # set column with query contig name
  coords$qry <- rep("Query", nrow(values(coords)))

  return(coords)
}





faidx_to_GRanges <- function(faidx_file){
  faidx <- read.table(file = faidx_file, header = F, stringsAsFactors = F,
                      col.names = c("name", "contig.len", "offset", 
                                    "linebases", "linewidth"))
  
  # create a GRanges object for reference sequences
  grange <- GRanges(seqnames = Rle(faidx$name), 
                    ranges = IRanges(start = rep(1, nrow(faidx)), 
                                     end = faidx$contig.len))
  
  # add seqlengths to the reference GRanges object
  seqlengths(grange) <- faidx$contig.len
  genome(grange) <- "Reference"
  grange <- sortSeqlevels(grange)
  grange <- sort(grange)
  return(grange)
}



circular_plot_w_ref <- function(reference_GRange, NUCmer_coords){
  # reference_GRange: reference sequence GRanges obj
  # NUCmer_coords: a GRanges object produced by reading in a show-coords processed NUCmer object.
  
  library(ggplot2)
  library(ggbio)
  
  p <- ggbio() + 
    circle(NUCmer_coords, geom = "rect", 
           aes(color = "steelblue", fill = "steelblue")) +  # NUCmer obj
    circle(reference_GRange, geom = "ideo", 
           aes(color = "gray70", fill = "gray70")) +        # Ideogram of ref genome
    circle(reference_GRange, geom = "scale", 
           scale.type = "M", size = 1.5) +                  # Scale from seqlen of ref genome
    # circle(reference_GRange, geom = "text", 
    #       aes(label = seqnames), size = 2) +              # Uncomment for sequence label
    scale_color_manual(name = "Sequence Origin", 
                       labels = c("Reference", "Query"), 
                       values = c("gray70", "steelblue")) +
    scale_fill_manual(name = "Sequence Origin", 
                      labels = c("Reference", "Query"), 
                      values = c("gray70", "steelblue")) +
    ggtitle("Reference vs. Query")
  return(p)
}



load_and_plot_nucmer_w_ref <- function(NUCmer_coords_file, ref_faidx_file, perc.id) {
  # NUCmer_coords_file: string for file path of NUCmer output produced by show coords 
  # ref_faidx_file: reference sequence GRanges obj
  # perc.id: percent id cutoff to show on plot
  
  library(GenomicRanges)
  
  NUCmer_coords <- load_coords(NUCmer_coords_file, perc.id = perc.id)   # Make GRanges obj of nucmer output file
  referenceGR <- faidx_to_GRanges(faidx_file = ref_faidx_file)
  plot <- circular_plot_w_ref(reference_GRange = referenceGR, NUCmer_coords = NUCmer_coords)
  return(plot)
}



load_and_plot_nucmer_w_ref(NUCmer_coords_file = "nucmer.contigs_filter_coords.txt", 
                           ref_faidx_file = "scaffolds.reduced.fasta.fai", 
                           perc.id = 87)
