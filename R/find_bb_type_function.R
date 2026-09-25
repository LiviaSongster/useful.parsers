#' Find backbone enzyme types.
#'
#' This function is used within antismash_parseR to find backbone enzyme types.
#' for specific genes, parse the description column and extract the enzyme type.
#' split by "<br>" and keep biosynthetic regions.
#' @param index gene index in the table
#' @return vector list flattened by ;
#' @export

find_bb_type <- function(index,all_orfs) {
  temp <- unlist(strsplit(all_orfs$description[index],"<br>"))
  temp2 <- temp[grepl(".*biosynthetic \\(rule-based-clusters\\) ",temp)]
  head(temp2)
  # remove everything including and before "(rule-based-clusters) "
  temp3 <- sub(".*\\(rule-based-clusters\\) ", "", temp2)
  head(temp3)
  # split at the : and keep unique
  temp4 <- strsplit(temp3,": ")

  if (length(temp4) > 1) {
    # take the first element of each split and keep only unique ones
    output <- unique(sapply(temp4, `[`, 1))

  } else {
    # only one entry → take its first part
    output <- temp4[[1]][1]
  }

  # if output has more than one thing, flatten with a semicolon
  if (length(output) > 1) {
    output <- paste(output, collapse = ";")
  }
  return(output)
}
