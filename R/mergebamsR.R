#' Merge BAM files
#'
#' This function merges multiple BAM files into a single output file. It checks for the existence of each input BAM file 
#' and the output directory. The function allows optional customization of output names and prefixes.
#'
#' @param bams A vector of file paths for the BAM files to be merged.
#' @param out_path The directory path where the merged BAM file will be saved. The function will stop if the specified output path does not exist.
#' @param names Optional; a vector of names to assign to the merged BAM files. If not provided, the names will be set to empty list.
#' @param prefixes Optional; a vector of prefixes to prepend to the BAM file names during merging. If not provided, no prefixes are used.
#'
#' @return Does not return a value; it generates a merged BAM file at the specified output path.
#'
#' @examples
#' # Assuming you have valid paths to BAM files:
#' bam_files <- c("path/to/bam1.bam", "path/to/bam2.bam")
#' output_path <- "path/to/output/directory"
#' mergebams(bam_files, output_path)
#'
#'@references This documentation was written by ChatGPT v4 - OpenAI, conversation with the author, 5-1-2024.

mergebams<-function(bams, out_path, names=NULL, prefixes=NULL){
  exists<-sapply(bams, file.exists)
  if(!file.exists(out_path)){stop(paste0("Provided out_path not found: ", out_path))}
  if(is.null(prefixes)){
    prefixes = rep("", length(bams))
  }
  if(is.null(names)){
    names<-vector(mode = "list", length = length(bams))
  }
  if(all(exists)){
    mergebams_rust_helper(bams, out_path, names, prefixes)
  } else {
    message(paste0("Files not found:\n", paste(bams[!exists], collapse="\n")))
  }
  
}

subsetbam<-function(inputbam, tags, outputbams, prefixes = NULL, TAG="CB"){
  # if(length(inputbams)!=length(tags)) {stop("Input number of bam files is not equal to number of tags")}
  if(length(tabs)!=length(outputbams)) {stop("Input number of output bam files is not equal to number of output bams")}
  exists<-file.exists(bam)
  if(is.null(prefixes)){
    prefixes<-rep("", length(outputbams))
  }
  if(length(prefixes)!=length(outputbams)) {stop("Input number of prefixes is not equal to number of output bams")}
  if(exists){
    subsetbam_rust_helper(inputbam, tags, outputbams, prefixes)
  } else {
    message(paste0("File not found:\n\t", inputbam))
  }
}


#' Peek into a BAM File with Tag Filtering
#'
#' This function provides a quick look into the contents of a BAM file, filtered by a specific tag.
#' It validates the input parameters, checks for the file's existence, and utilizes a helper function to
#' handle the BAM processing if the file exists. An error is raised if the input conditions are not met.
#'
#' @param bam A character string specifying the path to a single BAM file.
#'            Only a single file should be specified.
#' @param n An integer, defaulting to 100, indicating the number of entries to peek.
#'          `n` must be greater than 0.
#' @param TAG A character string specifying the tag to filter by within the BAM file.
#'            The default is "CB" (Cell Barcode).
#'
#' @return The function itself does not return a value; it operates through side effects
#'         such as invoking a Rust helper function or printing messages to the console.
#'
#' @details If the file does not exist, a message will be displayed. If `n` is less than 1
#'          or if more than one BAM file is specified, the function will stop with an error.
#'
#' @examples
#' # Assuming 'example.bam' is a valid BAM file path:
#' peekbam("example.bam", n = 10, TAG = "CB")
#' 
#'@references This documentation was written by ChatGPT v4 - OpenAI, conversation with the author, 5-2-2024.
#'
peekbam <- function(bam, n=100, TAG="CB"){
  if(as.integer(n)<1){stop("n must be more than 1")}
  if(length(bam)>1){stop("More than one bam file supplied")}
  exists<-file.exists(bam)
  if(exists){
    peekbam_rust_helper(bam, n, TAG)
  } else {
    message(paste0("File not found:\n", paste(bam, collapse="\n")))
  }
  
}
