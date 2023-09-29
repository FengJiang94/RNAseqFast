#' Merge FeatureCount outputs into one table
#'
#' This function merges FeatureCounts outputs from RNAseqPE and summarizes the results
#'    into a text table.
#'
#' The output table can be used for downstream DE analysis.
#'
#' @param filepath Path to the folder storing FeatureCounts outputs. Use the "raw_counts" folder in RNAseqPE outputs.
#' @param Assay File prefix of the output table".
#' @param colDataPath File name of the Sample information text file. This file contain at least three columns:
#' sample_name, short_name, condition_1, condition_2, ... Sample_name should be the same as the "stat/sample,csv"
#' in RNAseqPE outputs.
#' @param write_to_file True of False, write to a file or not. If True, then use assay as file prefix.
#' @return Merged FeatureCounts table.
#' @export


####################### generate count_matrix ############################
# FeatureCounts outputs are multiple .txt files (each one for one sample)
MergeFeatureCounts<-function(filepath, Assay, colDataPath, write_to_file = T)
{
  #get the file names
  currentpath<-getwd()
  setwd(filepath)
  file_list<-list.files(pattern = "*.txt$")
  cat(paste("load file:\n", file_list, collapse = "\n", sep = ""))
  #read all the files into one list of matirx
  counts_list<- lapply(file_list, function(x) read.csv(x, skip = 1, header = TRUE,row.names= 1, sep= "\t"))
  setwd(currentpath)
  cat("\n")
  cat("load file done")
  # clean up the matrixs in the list.(only keep c=the col coataining raw counts)
  for (i in 1:length(counts_list)){
    counts_list[[i]][,1:5]<-NULL
    names<-as.character(unlist(strsplit(colnames(counts_list[[i]]), split = ".",fixed = TRUE))[1])
  }
  head(counts_list[[1]])
  # convert the list of matrix into one matrix (rownames= gene id, colnames= sample id)
  counts_matrix<-counts_list[[1]]
  for (i in 2:length(counts_list)){
    counts_matrix<-merge(counts_matrix, counts_list[[i]], by ='row.names')
    rownames(counts_matrix)<-counts_matrix$Row.names
    counts_matrix$Row.names<-NULL
  }
  # change the sample IDs into sample names (V300069051_ to sample ID)
  sample_info<-read.csv(colDataPath,header=T)
  sample_info$sample_name <- paste(sample_info$sample_name, "_mapped.sort.bam", sep = "")
  # reorder columns in counts_matrix to match the order of sample_info$sample.ID
  counts_matrix <- counts_matrix[,sample_info$sample_name]
  head(counts_matrix)
  # change the colnames of counts_matrix to sample.ID
  colnames(counts_matrix)<-sample_info$short_name
  if ( write_to_file == T){
    write.csv(counts_matrix,file = paste(Assay, "counts_matrix.csv", sep = "_"))
  }
  return (counts_matrix)
}



