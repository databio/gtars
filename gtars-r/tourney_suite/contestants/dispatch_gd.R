#!/usr/bin/env Rscript
# Shared dispatcher for GenomicDistributions vs gtars comparison.
# Usage: Rscript dispatch_gd.R <params.json> <output_path> <backend>
#
# backend = "gd" (GenomicDistributions R package) or "gtars" (Rust)
#
# Data helpers (getGeneModels, getChromSizes, getGenomeBins) always come from
# GenomicDistributions — gtars doesn't export these. Only the computation
# functions differ between backends.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicDistributions)
  library(gtars)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: dispatch_gd.R <params.json> <output_path> <backend>")
}

params_path  <- args[1]
output_path  <- args[2]
backend      <- args[3]

params <- fromJSON(params_path)
func   <- params[["function"]]

# --- Test data ---
data("vistaEnhancers", package = "episcope")

# --- Dispatch ---

switch(func,

  calcWidth = {
    if (backend == "gd") {
      result <- GenomicDistributions::calcWidth(vistaEnhancers)
    } else {
      result <- gtars::calcWidth(vistaEnhancers)
    }
    writeLines(as.character(result), output_path)
  },

  calcNeighborDist = {
    if (backend == "gd") {
      result <- GenomicDistributions::calcNeighborDist(vistaEnhancers)
    } else {
      result <- gtars::calcNeighborDist(vistaEnhancers)
    }
    writeLines(as.character(result), output_path)
  },

  calcNearestNeighbors = {
    if (backend == "gd") {
      result <- GenomicDistributions::calcNearestNeighbors(vistaEnhancers)
    } else {
      result <- gtars::calcNearestNeighbors(vistaEnhancers)
    }
    writeLines(as.character(result), output_path)
  },

  calcFeatureDist = {
    geneModels <- GenomicDistributions::getGeneModels("hg19")
    tss <- GenomicRanges::promoters(geneModels$genesGR, upstream = 0, downstream = 1)
    if (backend == "gd") {
      result <- GenomicDistributions::calcFeatureDist(vistaEnhancers, tss)
    } else {
      result <- gtars::calcFeatureDist(vistaEnhancers, tss)
    }
    out <- ifelse(is.na(result), "NA", as.character(result))
    writeLines(out, output_path)
  },

  calcPartitions = {
    geneModels <- GenomicDistributions::getGeneModels("hg19")
    if (backend == "gd") {
      partitionList <- GenomicDistributions::genomePartitionList(
        geneModels$genesGR, geneModels$exonsGR,
        geneModels$threeUTRGR, geneModels$fiveUTRGR
      )
      result <- GenomicDistributions::calcPartitions(vistaEnhancers, partitionList)
    } else {
      partitionList <- gtars::genomePartitionList(
        geneModels$genesGR, geneModels$exonsGR,
        threeUTRGR = geneModels$threeUTRGR,
        fiveUTRGR = geneModels$fiveUTRGR
      )
      result <- gtars::calcPartitions(vistaEnhancers, partitionList)
    }
    fwrite(result, output_path, sep = "\t")
  },

  calcExpectedPartitions = {
    geneModels <- GenomicDistributions::getGeneModels("hg19")
    chromSizes <- GenomicDistributions::getChromSizes("hg19")
    if (backend == "gd") {
      genomeSize <- sum(as.numeric(chromSizes))
      partitionList <- GenomicDistributions::genomePartitionList(
        geneModels$genesGR, geneModels$exonsGR,
        geneModels$threeUTRGR, geneModels$fiveUTRGR
      )
      result <- GenomicDistributions::calcExpectedPartitions(
        vistaEnhancers, partitionList, genomeSize = genomeSize
      )
    } else {
      partitionList <- gtars::genomePartitionList(
        geneModels$genesGR, geneModels$exonsGR,
        threeUTRGR = geneModels$threeUTRGR,
        fiveUTRGR = geneModels$fiveUTRGR,
        chromSizes = chromSizes
      )
      result <- gtars::calcExpectedPartitions(
        vistaEnhancers, partitionList, genomeSize = chromSizes
      )
    }
    fwrite(result, output_path, sep = "\t")
  },

  calcSummarySignal = {
    data("exampleOpenSignalMatrix_hg19", package = "episcope")
    if (backend == "gd") {
      result <- GenomicDistributions::calcSummarySignal(
        vistaEnhancers, exampleOpenSignalMatrix_hg19
      )
    } else {
      result <- gtars::calcSummarySignal(
        vistaEnhancers, exampleOpenSignalMatrix_hg19
      )
    }
    fwrite(result$signalSummaryMatrix, output_path, sep = "\t")
  },

  calcChromBins = {
    chromSizes <- GenomicDistributions::getChromSizes("hg19")
    if (backend == "gd") {
      genomeBins <- GenomicDistributions::getGenomeBins(chromSizes)
      result <- GenomicDistributions::calcChromBins(vistaEnhancers, genomeBins)
    } else {
      result <- gtars::regionDistribution(vistaEnhancers, 250L)
    }
    fwrite(result, output_path, sep = "\t")
  },

  stop(paste("Unknown function:", func))
)
