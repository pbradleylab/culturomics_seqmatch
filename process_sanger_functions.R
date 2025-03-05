
library(sangeranalyseR)
# library(ape)

load_and_trim <- function(data_dir="./seqs",
                          fn_remove=c("Christian_",
                                      "515F_718447_",
                                      ".ab1"),
                          readFeature="Forward Read",
                          readFileName = p,
                          TrimmingMethod = "M2",
                          M2CutoffQualityScore = 25,
                          M2SlidingWindowSize = 15,
                          ...) {
  # list files and extract sample names
  
  seq_files <- list.files(data_dir, ".*.ab1")
  seq_names <- vapply(seq_files, \(x) {
    for (r in fn_remove) {
      x <- gsub(r, "", x)
    }
    x
  }, "")
  names(seq_files) <- seq_names
  
  # Load and trim Sanger reads
  trimmed_reads <- lapply(seq_files, \(sf) {
    p <- file.path(data_dir, sf)
    sangeranalyseR::SangerRead(
      readFeature = readFeature,
      readFileName = p,
      TrimmingMethod = TrimmingMethod,
      M2CutoffQualityScore = M2CutoffQualityScore,
      M2SlidingWindowSize = M2SlidingWindowSize,
      ...)
  })
  
  return(trimmed_reads)
  
}

separate_taxonomy_with_s <- function(inpt, taxa_col, remove=FALSE) {
  inpt <- inpt %>%
    tidyr::separate_wider_delim({{ taxa_col }},
                                names = c("d", "p", "c", "o", "f", "g", "s"),
                                delim = ";",
                                cols_remove=remove) %>%
    mutate(across(c("d", "p", "c", "o", "f", "g", "s"), ~ gsub("[dpcofgs]__","", .))) %>%
    rename_with(~ case_when(
      . == "d" ~ "Domain",
      . == "p" ~ "Phylum",
      . == "c" ~ "Class",
      . == "o" ~ "Order",
      . == "f" ~ "Family",
      . == "g" ~ "Genus",
      . == "s" ~ "Species",
      TRUE ~ .
    ))
  return(inpt)
}

manually_trim <- function(sanger_read, trim5=NULL, trim3=NULL) {
  sanger_read@QualityReport@TrimmingMethod <- "Manual"
  sanger_read@QualityReport@trimmedFinishPos <- sanger_read@QualityReport@rawSeqLength - trim3
  sanger_read@QualityReport@trimmedStartPos <- trim5
  sanger_read
}



interactive_alignment <- function(name="F_F01", trimmed_list=trimmed, vsearch=vsearch_processed) {
  # maximum three tied alignments
  new_cols <- c("rgba(80,120,255,1)", "rgba(80,255,120,1)", "rgba(255,120,80)")
  l <- plotly::layout(
    sangeranalyseR::qualityBasePlot(trimmed_list[[name]]),
    title = list(text=name, x=.05)
  )
  components <- map_chr(l$x$data, ~ .x$name)
  which_trimmed <- which(components=="Trimmed Read")
  
  aln_read_coords <- dplyr::filter(vsearch, qseqid==name)
  for (i in 1:(min(3, nrow(aln_read_coords)))) {
    aln_read <- l$x$data[[which_trimmed]]
    aln_read$name <- paste0("Aligned (", aln_read_coords[i, "Species"][[1]], ")")
    aln_read$line$color <- new_cols[i]
    aln_read$x <- aln_read$x[aln_read_coords[i, "qstart"][[1]] : aln_read_coords[i, "qend"][[1]]]
    aln_read$y <- aln_read$y[aln_read_coords[i, "qstart"][[1]] : aln_read_coords[i, "qend"][[1]]]+(i*3)
    aln_read$text <- rep(paste0(aln_read_coords[["pident"]][i], "% ID<br>", aln_read_coords[["gtdb_taxonomy"]][i]), length(aln_read$x))
    n_comps <- length(l$x$data)
    l$x$data[[n_comps + i]] <- aln_read
  }
  l
}
