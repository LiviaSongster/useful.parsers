# useful.parsers

Author: Livia Oster, songs005\@umn.edu

This R package contains useful parsers for the outputs of busco, quast, and antismash.

## Installation

```{r}
pak::pak("LiviaSongster/useful.parsers")

# library and test it
library(useful.parsers)

# see functions in the package
getNamespaceExports("useful.parsers")
# there should be 5:
# [1] "antismash_parseR" "import_gbk"       "find_bb_type"     "parse_busco"     
# [5] "parse_quast"   
# if not, Session > Restart and try install again

# check the help page for each function by typing:
?antismash_parseR()
?parse_busco()
?parse_quast()

# these two are functions that are used within antismash_parseR
?import_gbk()
?find_bb_type()
```

## Example usage

## parse_busco

```{r}
# --------------------------
# 1. parse busco output
# --------------------------

# --------------------------
# 1A. parse a single file:

dir="/Users/songs005/Documents/Scripts/R_packages/useful.parsers/data"
file_name=paste0(dir,"/busco/NRRL3249_short_summary.txt")
busco_output <- parse_busco(file_name)

# check the result
head(busco_output)

# --------------------------
# 2B. batch process a dir full of busco summary tables:
# first get list of files in the dir
busco_list <- list.files(paste0(dir,"/busco"),full.names = TRUE)

# then run parse_busco on everything and bind it together
# this will also return an error message on any files that did not parse properly
busco_batch_output <- do.call(rbind, lapply(busco_list, function(x) {
  tryCatch({
    parse_busco(x)  # Attempt to run parse_busco
  }, error = function(e) {
    message("Error in file: ", x)  # Print the file name causing the error
    return(NULL)  # Return NULL so it can be filtered out later
  })
}))

# check the result
head(busco_batch_output)
unique(busco_batch_output$Accession)
```

## parse_quast

NOTE - the values parsed here will all be character, NOT numeric, so make sure to fix that downstream if you need to do any averaging or plotting.

```{r}
# --------------------------
# 2. parse quast output
# --------------------------

# --------------------------
# 2A. parse a single file:

dir="/Users/songs005/Documents/Scripts/R_packages/useful.parsers/data"
file_name=paste0(dir,"/quast/MTD517_report.tsv")
quast_output <- parse_quast(file_name)

# check the result
head(quast_output)

# --------------------------
# 2B. batch process a dir full of quast outputs:

# first get list of files in the dir
quast_list <- list.files(paste0(dir,"/quast"),full.names = TRUE)

# then run parse_quast on everything and bind it together
# this will also return an error message on any files that did not parse properly
quast_batch_output <- do.call(rbind, lapply(quast_list, function(x) {
  tryCatch({
    parse_quast(x)  # Attempt to run parse_quast
  }, error = function(e) {
    message("Error in file: ", x)  # Print the file name causing the error
    return(NULL)  # Return NULL so it can be filtered out later
  })
}))

# check the result
head(quast_batch_output)
unique(quast_batch_output$Accession)

```

## antismash_parseR

```{r}
# --------------------------
# 3. parse antismash output
# --------------------------
dir="/Users/songs005/Documents/Scripts/R_packages/useful.parsers/data"
dir_name=paste0(dir,"/antismash/NRRL3510_antismash_7.1.0")

antismash_parseR(dir_name,
                 output_dir=paste0(dir,"/antismash/parseR_output"),
                 output_gbks=TRUE)

# process a list of dirs
anti_dir_list <- list.dirs(paste0(dir,"/antismash"),
                           recursive=FALSE,
                           full.names = TRUE)

# remove the parseR output dir
anti_dir_list <- anti_dir_list[
  basename(anti_dir_list) != "parseR_output"
]

lapply(anti_dir_list, function(x) {
  tryCatch({
    antismash_parseR(x,
                     output_dir=paste0(dir,"/antismash/parseR_output"),
                     output_gbks=TRUE)
  }, error = function(e) {
    message("Error in file: ", x)  # Print the file name causing the error
    return(NULL)  # Return NULL so it can be filtered out later
  })
})

```
