#' MergebamsR
#' @export
#' 

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

