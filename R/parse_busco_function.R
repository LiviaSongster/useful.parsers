#' Parse BUSCO output.
#'
#' This function parses individual busco summary tables and imports them to R.
#'
#' @param file_name required file path; short_summary.txt output file from BUSCO
#' @param accession optional string; accession or unique ID for the output
#' @return A dataframe with the busco results
#' @export
parse_busco <- function(file_name,accession=NA) {

  # read in results
  temp <- readLines(file_name)
  file_basename <- basename(file_name)
  # if accession is not provided, attempt to strip it from the filename
  # first get basename of the file
  if (is.na(accession)) {
    accession <- gsub("_short_summary.txt","",file_basename)
  }

  # get busco information
  version <- sub("^# BUSCO version is:\\s*", "", temp[1])
  lineage <- sub("^# The lineage dataset is:\\s*", "", temp[2])
  input <- sub("^# Summarized benchmarking in BUSCO notation for file\\s*", "", temp[3])
  mode <- sub("^# BUSCO was run in mode:\\s*", "", temp[4])

  # extract values
  c <- as.numeric(sapply(temp[grep("Complete BUSCOs", temp)], function(x) strsplit(x, "\t")[[1]][2]))
  s <- as.numeric(sapply(temp[grep("Complete and single-copy BUSCOs", temp)], function(x) strsplit(x, "\t")[[1]][2]))
  d <- as.numeric(sapply(temp[grep("Complete and duplicated BUSCOs", temp)], function(x) strsplit(x, "\t")[[1]][2]))
  f <- as.numeric(sapply(temp[grep("Fragmented BUSCOs", temp)], function(x) strsplit(x, "\t")[[1]][2]))
  m <- as.numeric(sapply(temp[grep("Missing BUSCOs", temp)], function(x) strsplit(x, "\t")[[1]][2]))
  t <- as.numeric(sapply(temp[grep("Total BUSCO groups searched", temp)], function(x) strsplit(x, "\t")[[1]][2]))

  # if busco was run in genome mode, parse assembly statistics:
  if (mode == "genome") {
    scaf <- as.numeric(sapply(temp[grep("Number of scaffolds", temp)], function(x) strsplit(x, "\t")[[1]][2]))
    con <- as.numeric(sapply(temp[grep("Number of contigs", temp)], function(x) strsplit(x, "\t")[[1]][2]))
    totlen <- as.numeric(sapply(temp[grep("Total length", temp)], function(x) strsplit(x, "\t")[[1]][2]))
    gap <- as.numeric(sapply(temp[grep("Percent gaps", temp)], function(x) strsplit(x, "\t")[[1]][2]))
    n50s <- as.numeric(sapply(temp[grep("Scaffold N50", temp)], function(x) strsplit(x, "\t")[[1]][2]))
    n50c <- as.numeric(sapply(temp[grep("Contigs N50", temp)], function(x) strsplit(x, "\t")[[1]][2]))

    # prepare output table
    output <- data.frame(
      Accession = accession,
      file_name = file_basename,
      Version = version,
      Lineage = lineage,
      Input = input,
      Mode = mode,
      Percent_single_copy = s / t * 100,
      Complete_buscos = c,
      Complete_single = s,
      Complete_duplicated = d,
      Fragmented = f,
      Missing = m,
      Total_buscos = t,
      Scaffolds = scaf,
      Contigs = con,
      Total_length = totlen,
      Perc_gaps = gap,
      Scaffold_N50 = n50s,
      Contig_N50 = n50c
    )
  } else {
    output <- data.frame(
      Accession = accession,
      file_name = file_basename,
      Version = version,
      Lineage = lineage,
      Input = input,
      Mode = mode,
      Percent_single_copy = s / t * 100,
      Complete_buscos = c,
      Complete_single = s,
      Complete_duplicated = d,
      Fragmented = f,
      Missing = m,
      Total_buscos = t
    )
  }

  return(output)
}
