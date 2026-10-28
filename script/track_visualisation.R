
# Author : Damien Dufour, University of Oslo, 2026 

# Load Packages and Sources -----------------------------------------------


install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    # Try installing with BiocManager first
    tryCatch({
      BiocManager::install(pkg, dependencies = TRUE, update = FALSE, lib="/storage/mathelierarea/processed/local_homes/damiend/R/x86_64-pc-linux-gnu-library/4.5/")
    }, error = function(e) {
      message(paste("BiocManager could not install", pkg, "- trying install.packages instead."))
      install.packages(pkg, dependencies = TRUE, lib="/storage/mathelierarea/processed/local_homes/damiend/R/x86_64-pc-linux-gnu-library/4.5/")
    })
    
    # Load the package after installation
    library(pkg, character.only = TRUE)
  }
}
packages <- c(
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "rtracklayer",
    "Gviz",
    "GenomicRanges",
    "dplyr",
    "stringr",
    "GenomicInteractions",
    "argparse", # to fetch argument from shell
    "conflicted", # filter and select are ambiguous
    "viridisLite"
 )

suppressWarnings(
  sapply(packages, install_and_load)
  )


conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
# Load Data ---------------------------------------------------------------


# parse data
parser <- ArgumentParser(description= 'This program plot genomic track based on bed files')
parser$add_argument('--input', '-i', help= 'Input file')
parser$add_argument('--output_dir', '-o', help= 'Output directory')
parser$add_argument('--coordinates', '-co', help= 'Coordinates of the region to visualise')
parser$add_argument('--modality', '-mo', help= 'What modality to plot')
parser$add_argument('--max', '-ma', help= 'Maximum value on the plot for each modality')
parser$add_argument('--ttt',  help= 'Which treatment to include')
parser$add_argument('--day', "-d", help= 'Which day to include')
parser$add_argument('--title', '-t', help= 'Name of the gene or region to plot')
parser$add_argument('--interaction', help= 'Should the interaction file be plotted', type= 'logical')
parser$add_argument('--interaction_file', help= 'Location of bedpe')
parser$add_argument('--genome', default='hg38', help= 'Reference genome')

xargs<- parser$parse_args()

if (!file.exists(xargs$input)) 
  {
    stop(paste("Error1: The file", xargs$input, "does not exist."))
  } 

sample_sheet <- read.csv(xargs$input, 
                         strip.white = TRUE,
                         header = T)

# Process arguments (always read as strings) 
gene_name  <- xargs$title
genes_loc <- xargs$coordinates

max <- xargs$max%>%
  str_split(",") %>%
  unlist() %>%
  as.numeric() 

day <- xargs$day  %>%
  str_split(",") %>%
  unlist() 

ttt <- xargs$ttt  %>%
  str_split(",") %>%
  unlist() 

modality <- xargs$modality %>%
  str_split(",") %>%
  unlist() 

Interaction <- xargs$interaction
interaction_file <- xargs$interaction_file

# load genome 

 genome <- xargs$genome

if (genome == "hg38") {
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
} else if (genome == "hg19") {
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
} else if (genome == "mm10") {
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
} else if (genome == "mm9") {
  library(TxDb.Mmusculus.UCSC.mm9.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
} else {
  stop("Unsupported genome version provided!")
}



# subset sample_sheet based on args ---------------------------------------

# # Import files -------------------------
# # https://ivanek.github.io/Gviz/reference/settings.html

 filtered_samples <- sample_sheet %>%
      filter(input_modality %in% modality) %>%
      filter(input_day %in% day) %>%
      filter(input_ttt %in% ttt) 
      
bw_files <- filtered_samples %>%
      pull(path)

 sample_names <- paste(filtered_samples$input_day,
                         filtered_samples$input_ttt,
                         sep = "_")

for(i in 1:length(bw_files))
{
if (!file.exists(bw_files[i])) 
  {
    stop(paste("Error2: The file", bw_gr[i], "does not exist."))
  }

print("bigWig files exist!")

tryCatch({
  bw_info <- seqinfo(BigWigFile(bw_files[i]))
  print("BigWig file is correct and readable!")
  print(bw_info)}, error = function(error2) {
  message("Error3: The BigWig file is corrupted or invalid: ", error2$message)
})

}

bw_gr <- lapply(bw_files, function(f) import(f,
                                              format = "bigWig"))

cat("bigWig snippet")
head(bw_gr[1]) 

 chrom  <- str_split_i(genes_loc,":", 1)

 start  <- str_split_i(genes_loc,":", 2) %>%
   str_split_i("-",1) %>%
   str_replace_all(",","") %>%
   as.numeric()

 end    <- str_split_i(genes_loc,":", 2) %>%
   str_split_i("-",2) %>%
   str_replace_all(",","")%>%
   as.numeric()


# define colours

unique_combos <- unique(sample_names)
num_combos    <- length(unique_combos)

base_palette <- viridis(num_combos)

colour_map <- setNames(base_palette, unique_combos)

assigned_colours <- colour_map[sample_names]


# # Create one DataTrack per sample

data_tracks <- Map(
  f = function(gr, nm, colours_plot,max_val) {
    DataTrack(
      range      = gr,
      genome     = genome,
      chromosome = chrom,
      type       = "histogram",      # "l" = line, "histogram" = bars
      name       = nm,
      fill = colours_plot,
      col.frame = "transparent",
      ylim = c(0,max_val),
      col.histogram = "transparent" # change per sample if you like
    )
  },
  max_val = max,
  gr = bw_gr,
  colours_plot = assigned_colours,
  nm = sample_names
)


cat("data_tracks built test")
cat(data_tracks[[1]]@data[1,1])


 itrack <- IdeogramTrack(genome = genome,
                         chromosome = chrom)

 gtrack <- GenomeAxisTrack()



genes  <- GeneRegionTrack(txdb,
                          genome = genome,
                          chromosome = chrom,
                          start = start,
                          end = end,
                          name = "gene",
                          showId=TRUE,
                          geneSymbol=TRUE,
                          col = "black",
                          col.line = "black",
                          col.symbol = "black",
                          fontsize.group = 5)




# check if tracks are ok 

tryCatch(
  {
  plotTracks(genes, chromosome = chrom, from = start, to = end) 
  }, 
  error = function(error4) {
  
  message("error4 : the 'genes' track is not plotting: ", error4$message)
    
})

tryCatch(
  {
  plotTracks(gtrack, chromosome = chrom, from = start, to = end) 
  }, 
  error = function(error4) {
  
  message("error4 : the 'gtrack' track is not plotting: ", error4$message)
    
})

tryCatch(
  {
  plotTracks(data_tracks, chromosome = chrom, from = start, to = end) 
  }, 
  error = function(error4) {
  
  message("error4 : the 'data_tracks' track is not plotting: ", error4$message)
    
})


cat("Tracks ok, plotting...")


if(Interaction == F) {

  
  pdf(paste0(xargs$output,"/figures/",gene_name,"_track.pdf"),
      height = length(bw_files)*.4+0.4,
      width = 3.89,
      pointsize = .5)
  
  plotTracks(
    c(data_tracks, 
      genes, 
      gtrack),
    chromosome = chrom,
    from       = start,
    to         = end,
    background.title = "transparent",
    col.title = "black",
    fontcolor = "black",
    showAxis = F,
    fontsize = 3,
    stacking="squish",
    collapseTranscripts  = "meta",   # or "meta", or "gene" or F
    transcriptAnnotation = "symbol",     # show gene symbols if available
    sizes      = c(rep(.4, length(data_tracks)), 
                   1,
                   .3),
    lwd = .5
  )
  
  dev.off()
}


if(Interaction == T) {
if (!file.exists(interaction_file)) 
  {
    stop(paste("Error5: The file", interaction_file, "does not exist."))
  } 
  
   bedpe <- read.table(interaction_file,
                      header = F,
                      stringsAsFactors = FALSE)

  colnames(bedpe)[1:7] <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "count")

  anchor1 <- GRanges(
    seqnames = bedpe$chr1,
    ranges   = IRanges(bedpe$start1 + 1, bedpe$end1)  
  )
  
  anchor2 <- GRanges(
    seqnames = bedpe$chr2,
    ranges   = IRanges(bedpe$start2 + 1, bedpe$end2)
  )

  gi <- GenomicInteractions(anchor1, anchor2)
  
  
  mcols(gi)$score <- bedpe$count
  
  loopTrack <- InteractionTrack(gi)
  
  # Optional styling
  displayPars(loopTrack) <- list(
    col.outside = "#000000",
    col.interactions = "#000000",
    plot.outside = TRUE  # arcs drawn above panel if needed
  )

pdf(paste0(xargs$output,"/figures/",gene_name,"_track.pdf"),
    height = length(bw_files)*.4+.8,
    width = 3.89,
    pointsize = .5)

plotTracks(
  c(loopTrack, 
    data_tracks, 
    genes, 
    gtrack),
  chromosome = chrom,
  from       = start,
  to         = end,
  background.title = "transparent",
  col.title = "black",
  fontcolor = "black",
  showAxis = F,
  fontsize = 3,
  stacking="squish",
  collapseTranscripts  = "meta",   # or "meta", or "gene" or F
  transcriptAnnotation = "symbol",     # show gene symbols if available
  sizes      = c(.5,
                 rep(.4, length(data_tracks)), 
                 .4,
                 .4),
  lwd = .5
)

dev.off()
}

 write.csv(filtered_samples,
           paste0(xargs$output,"/data/",gene_name,"_input_data.csv"))

cat(paste0(gene_name," done"))












