#' Parse Quast output.
#'
#' This function parses individual quast output tables and imports them to R.
#'
#' @param file_name required file path; report.tsv output file from quast
#' @param accession optional string; accession or unique ID for the output
#' @return A dataframe with the quast statistics
#' @export
parse_quast <- function(file_name,accession=NA) {
  # read in file
  temp <- read.delim(file_name,sep="\t", comment.char = "")
  file_basename <- basename(file_name)
  # if accession is not provided, attempt to strip it from the filename
  # first get basename of the file
  if (is.na(accession)) {
    accession <- gsub("_report.tsv","",file_basename)
  }

  # check if there are two or three columns in the output
  # this will depend on the way quast was run
  if (ncol(temp)==2) {
    colnames(temp) <- c("Assembly","Scaffolds")
  } else if (ncol(temp)==3) {
    colnames(temp) <- c("Assembly","Scaffolds","Scaffolds_broken")
  }
  # add file name and accession as a new column
  temp$file_name <- file_basename
  temp$Accession <- accession

  # reorder columns
  if (ncol(temp)==4) {
    output <- temp[,c("Accession","file_name","Assembly","Scaffolds")]
  } else if (ncol(temp)==5) {
    output <- temp[,c("Accession","file_name","Assembly","Scaffolds","Scaffolds_broken")]
  }

  # return the data

  return(output)

}
