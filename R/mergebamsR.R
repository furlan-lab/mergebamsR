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
#'@export

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


#' Subset a BAM file based on specific tags
#'
#' This function subselects reads from a BAM file based on the provided tags and 
#' writes the subselected reads into separate output BAM files. The subselection 
#' is done using a specified tag (default "CB"). This function can utilize multiple cores
#' to perform the operation more efficiently.
#'
#' @param inputbam A character string specifying the path to the input BAM file.
#' @param tags A character vector of tags to filter the reads by.
#' @param outputbams A character vector specifying the paths to the output BAM files.
#' @param prefixes Optional; a character vector of prefixes to add to the read names in the output BAM files.
#'        If NULL, no prefixes are added. Defaults to NULL.
#' @param TAG Optional; the tag used for subselection, default is "CB".
#' @param cores Optional; the number of cores to use for the process, default is 1.
#' @param verbose Option to increase verbosity of output
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' subsetbam(inputbam = "path/to/input.bam", 
#'           tags = c("TAG1", "TAG2"), 
#'           outputbams = c("path/to/output1.bam", "path/to/output2.bam"),
#'           prefixes = c("prefix1", "prefix2"),
#'           TAG = "CB",
#'           cores = 2)
#'
#' @details
#' It's important that the length of `tags` is equal to the length of `outputbams`.
#' If `prefixes` is provided, its length must also match the length of `outputbams`.
#' If any of these conditions is not met, the function will stop and throw an error.
#' @export

subsetbam<-function(inputbam, tags, outputbams, prefixes = NULL, TAG="CB", cores=1, verbose=F){
  # if(length(inputbams)!=length(tags)) {stop("Input number of bam files is not equal to number of tags")}
  if(length(tags)!=length(outputbams)) {stop("Input number of output bam files is not equal to number of output bams")}
  exists<-file.exists(inputbam)
  if(verbose){
    message(paste0("Found file: ", inputbam, "\n"))
  }
  if(is.null(prefixes)){
    prefixes<-rep("", length(outputbams))
    if(verbose){
      message(paste0("No prefixes supplied"))
    }
  } else {
    if(verbose){
      message(paste0("Found ", length(prefixes), " prefixes"))
    }
  }
  if(length(prefixes)!=length(outputbams)) {stop("Input number of prefixes is not equal to number of output bams")}
  if(any(sapply(outputbams, file.exists))) {stop("One of the outputbam file exists.  Remove it and rerun subsetbam")}
  if(exists){
    if(verbose){
      message(paste0("Running subset_bam using TAG = ", TAG, "with ", cores, " core(s)"))
    }
    nc<-pbmcapply::pbmclapply(1:length(tags), function(i){
      subsetbam_rust_helper(inputbam = inputbam, tags = tags[i], outputbams = outputbams[i], prefixes = prefixes[i], tag = TAG, cores = 1)
    }, mc.cores = cores)
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
#'@export
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
