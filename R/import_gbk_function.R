#' Import antismash gbk files.
#'
#' This function is used within antismash_parseR to import antismash gbk files and extract ORF info.
#' @param gbk_file path to gbk_file
#' @return genes and associated info extracted from the gbk
#' @export

import_gbk <- function(gbk_file) {
  feature_list=c("gene")
  gb <- geneviewer::read_gbk(gbk_file,features=feature_list)
  # get feature details for each gene
  # first: name of cluster
  name <- names(gb)
  names(gb)
  # generalize it
  genelist <- gb[[1]]$FEATURES$gene

  # extract gene info
  df <- data.table::rbindlist(lapply(genelist, function(x) as.list(x)), use.names = TRUE, fill = TRUE)

  # repopulate Name column if necessary
  if (!"Name" %in% colnames(df)) {
    df$Name <- df$ID
  }

  genes_out <- df |>
    dplyr::select(Name, ID, region) |>
    dplyr::mutate(Name = dplyr::if_else(is.na(Name), ID, Name))

  # clean up the column names
  colnames(genes_out) <- c("locus_tag","gene","start_end_strand")

  # split name and add columns that will match all_orfs and all_products files
  chrm <- unlist(strsplit(name,".region"))[1]
  genes_out$chrm <- chrm
  region <- as.numeric(unlist(strsplit(name,".region"))[2])
  genes_out$region <- region
  # reorder columns
  genes_out <- as.data.frame(genes_out)
  order <- c("chrm","region","locus_tag","gene","start_end_strand")
  genes_out <- genes_out[,order]
  # print the BGC name in the console as the function runs
  # print(paste0("Parsed the .gbk file for: ",name))
  return(genes_out)
}
