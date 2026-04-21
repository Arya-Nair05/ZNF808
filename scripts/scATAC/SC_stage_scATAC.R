#-----------------------------------------------------------------------------
#          SC stage: scATAC-seq analysis using Signac 
##-----------------------------------------------------------------------------

## loading the libraries

  library(Seurat)
  library(Signac)
  library(AnnotationHub)
  library(tidyverse)
  library(GenomicRanges)
  library(viridis)
  library(ggvenn)
  library(ggrepel)
  library(patchwork)
  library(presto)
  
  set.seed(1234)


## Loading the master metadata file
  metadata <- read.csv(file = "scATAC_all_metadata.csv",
                       header = TRUE)


# Filtering for the SC stage samples:
  metadata <- metadata[grep("^SC", metadata$Sample),]

# adding sample name to the barcode since a few barcodes are shared between different samples
  metadata$barcode <- paste0(metadata$Sample,"_", metadata$barcode)
  metadata <- metadata %>% 
    column_to_rownames("barcode")


## Creating Fragments objects for each sample

# Saving the fragment file path
  f1_path <- "/projects/data/znf808/scatacseq/sc/SC_WT_1/outs/fragments.tsv.gz"
  f2_path <- "/projects/data/znf808/scatacseq/sc/SC_WT_2/outs/fragments.tsv.gz"
  f3_path <- "/projects/data/znf808/scatacseq/sc/SC_KO_1/outs/fragments.tsv.gz"
  f4_path <- "/projects/data/znf808/scatacseq/sc/SC_KO_2/outs/fragments.tsv.gz"

# For each sample's fragment file creating a fragments object and modifying the cell barcode by adding the sample name as prefix

  path_list <- list("SC_WT_1"= f1_path, "SC_WT_2"=f2_path, "SC_KO_1"=f3_path, "SC_KO_2"=f4_path)
  frag_objects_list <- list()
  
  for (i in names(path_list)){
    
    fragpath <- path_list[[i]]
    print(i)
    print(fragpath)
    total_counts <-  CountFragments(fragpath)
    cutoff <- 500 
    barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB
    
    # Create a fragment object
    frags <- CreateFragmentObject(path = fragpath, cells = barcodes)
    
    
    # adding the sample name as a prefix to the cell barcodes of the fragments object
    cell_names <- Cells(frags)
    new_cell_names <- paste0(i, "_", cell_names)
    names(new_cell_names) <- cell_names
    frags <- RenameCells(frags, new_cell_names)
    
    
    frag_objects_list[[i]] <- frags
    
  }

## Calling peaks on the list of 4 fragments file:
  
  peaks <- CallPeaks(path_list, 
                     macs2.path = "/penrose/projects/ana207/ls/envs/signac_macs/bin/macs3",
                     format="BEDPE")

## Loading the hg38 blacklist regions bed file
## Filtering the peaks overlapping with the blacklist regions
  # 1355 peaks were filtered out
  
  peaks <- subsetByOverlaps(peaks, blacklist_hg38_unified, invert=TRUE) 

## Creating feature barcode matrix
  # Merging the new cell barcode names of each of the fragments objects

  all_new_barcodes <- unlist(lapply(frag_objects_list,Cells))
  
  counts_BL <- FeatureMatrix(fragments = frag_objects_list,
                          features = peaks_copy,
                          cells = all_new_barcodes)

## Creating chromatin assay object
  chrom_assay_BL <- CreateChromatinAssay(counts = counts_BL,
    sep = c("-", "-"),
    fragments = frag_objects_list,
    min.cells = 10,
    min.features = 200)
  
  chrom_SC_rep_BL <- CreateSeuratObject(counts = chrom_assay_BL,
                                     assay = "peaks",
                                     meta.data = metadata)

# dim(chrom_SC_rep): 557901  27320

## Checking if there are peaks mapped to the chromosome scaffolds
  seqnames(granges(chrom_SC_rep_BL)) %in% standardChromosomes(granges(chrom_SC_rep_BL))
  # only 24 levels present

### Annotation
  ah <- AnnotationHub()
  query(ah, "EnsDb.Hsapiens.v98") ## latest version 113 available-- AH119325
  ensdb_v98 <- ah[["AH75011"]]
  annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
  genome(annotations) <- "hg38"
  
  Annotation(chrom_SC_rep_BL) <- annotations

# The peaks have the chromosomes in the UCSC format (chr1,chr2) while the ah has it in NCBI format (1,2)
# Converting seqlevels style of the annotation object from NCBI to UCSC

  ncbi_levels <- seqlevels(annotations) 
  ucsc_levels <- paste0("chr", ncbi_levels)
  
  Annotation(chrom_SC_rep_BL) <- renameSeqlevels(Annotation(chrom_SC_rep_BL), 
                                              setNames(ucsc_levels, ncbi_levels))

### QC metrics
  chrom_SC_rep_BL <- NucleosomeSignal(object = chrom_SC_rep_BL)
  chrom_SC_rep_BL <- TSSEnrichment(object = chrom_SC_rep_BL)
  chrom_SC_rep_BL$pct_reads_in_peaks <- chrom_SC_rep_BL$peak_region_fragments / chrom_SC_rep_BL$passed_filters * 100
  chrom_SC_rep_BL$nucleosome_group <- ifelse(chrom_SC_rep_BL$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

  chrom_SC_rep_BL$blacklist_ratio <- FractionCountsInRegion(object = chrom_SC_rep_BL,
                                                       assay = 'peaks',
                                                       regions = blacklist_hg38_unified)

#------------------------------------------------------------------------------------------------
#                        Visualising QC metrics
#------------------------------------------------------------------------------------------------

# TSS enrichment score
  DensityScatter(chrom_SC_rep_BL, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# Fragment length periodicity
  FragmentHistogram(object = chrom_SC_rep_BL) #,group.by = 'nucleosome_group')


  VlnPlot(chrom_SC_rep_BL,
          features = c('nCount_peaks','TSS.enrichment','blacklist_ratio',
                       'nucleosome_signal','pct_reads_in_peaks'), pt.size = 0.1, ncol = 5)

# Violin plots of all other QC metrics with threshold marked:
  nCount_peaks_thresh <- 15000
  tss_thresh <- 3
  blacklist_ratio_thresh <- 1e-3
  nucleosome_thresh <- 4
  reads_in_peak_thresh <- 50


  vlns <- VlnPlot(chrom_SC_rep_BL,
                  features = c('nCount_peaks','TSS.enrichment','blacklist_ratio',
                               'nucleosome_signal','pct_reads_in_peaks'), pt.size = 0.1,combine = FALSE)
  
  vlns <- Map(function(p, y) p + geom_hline(yintercept = y, color = "red"),
              vlns,c(nCount_peaks_thresh,
                     tss_thresh,
                     blacklist_ratio_thresh,
                     nucleosome_thresh,
                     reads_in_peak_thresh))

  wrap_plots(vlns, ncol = 5)


#------------------------------------------------------------------------------------------------
#                        Normalisation and clustering
#------------------------------------------------------------------------------------------------

## Filtering out all the low-quality cells
  chrom_SC_rep_BL_filtered <- subset(chrom_SC_rep_BL,
                           subset = nCount_peaks > nCount_peaks_thresh &
                             TSS.enrichment > tss_thresh &
                             blacklist_ratio < blacklist_ratio_thresh &
                             nucleosome_signal < nucleosome_thresh &
                             pct_reads_in_peaks > reads_in_peak_thresh)

# dim(chrom_SC_rep) after filtering: 557901  23509= 3811 cells removed

## Normalisation and dimensional reduction
  chrom_SC_rep_BL_filtered <- RunTFIDF(chrom_SC_rep_BL_filtered)
  chrom_SC_rep_BL_filtered <- FindTopFeatures(chrom_SC_rep_BL_filtered, min.cutoff = 'q0')
  chrom_SC_rep_BL_filtered <- RunSVD(chrom_SC_rep_BL_filtered, n = 80)

  DepthCor(chrom_SC_rep_BL_filtered, n = 80)
  
## Clustering
  ## when 3:80 dims are considered, 42 clusters are formed, out of which most clusters only have 1-4 cells
  ## so the num of dims are reduced to 40
  
## Clustering
  chrom_SC_rep_BL_filtered <- RunUMAP(object = chrom_SC_rep_BL_filtered,
                            reduction = 'lsi',
                            dims = 3:40,
                            min.dist = 0.01,
                            spread = 3)

  chrom_SC_rep_BL_filtered <- FindNeighbors(object = chrom_SC_rep_BL_filtered,
                                  reduction = 'lsi',
                                  dims = 3:40)

  chrom_SC_rep_BL_filtered <- FindClusters(object = chrom_SC_rep_BL_filtered,
                                 verbose = FALSE, 
                                 algorithm = 3,
                                 resolution = 0.5)

## Splitting the sample name into stage and genotype
  meta_data <- chrom_SC_rep_BL_filtered@meta.data
  meta_data <-  meta_data %>% 
    separate(Sample, into = c("Stage", "Genotype", "Rep"), remove = FALSE)
  
  chrom_SC_rep_BL_filtered@meta.data <- meta_data

## Saving the filtered and clustered object
  saveRDS(chrom_SC_rep_BL_filtered, file = "per_sample_chrom_SC_filtered_clustered_BL_removed.RDS")

## Viusalisation
  DimPlot(object = chrom_SC_rep_BL_filtered, label = TRUE)
  DimPlot(object = chrom_SC_rep_BL_filtered, group.by = "Genotype")
  DimPlot(object = chrom_SC_rep_BL_filtered, group.by = "Rep")
  

#------------------------------------------------------------------------------------------------  
#                       Finding the Marker Peaks
#------------------------------------------------------------------------------------------------  

# All marker peaks:
  da_peaks_bl_removed <- FindAllMarkers(object=chrom_SC_rep_BL_filtered, test.use="wilcox", min.pct=0.1)


## saving the list of marker genes
  write.csv(da_peaks_bl_removed, file = "per_sample_SC_DA_peaks_BL_removed.csv")




