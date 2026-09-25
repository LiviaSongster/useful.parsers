#' Parse antismash output.
#'
#' This function parses antismash output directories to find lists of genes, products, and backbone enzymes and imports them to R.
#' This function assumes you are running it on the antismash results directory.
#' Specifically, it parses the regions.js file and all of the region .gbk files.
#' These files must all be present and not compressed in the output directory.
#' Other files are not necessary.
#'
#' @param dir_name required dir path; antismash output directory
#' @param accession optional string; accession or unique ID for the output, for naming the outputs
#' @param output_dir optional dir path; directory for saving output files; defaults to directory specified in dir_name
#' @param output_product TRUE/FALSE; output products for each BGC region. defaults TRUE.
#' @param output_orfs TRUE/FALSE; output all orfs in each BGC region as specified in the regions.js file. Sometimes antismash changes locus tags for some genes. defaults TRUE.
#' @param output_gbks TRUE/FALSE; output all orfs in the GBK files for each BGC region. Contains extra border genes. Also contains unmodified locus tags for each gene. defaults FALSE.
#' @param include_seqs TRUE/FALSE; include protein translation in the output_orfs csv for each orf in the BGC region. defaults FALSE.
#' @return csv files in the output directory that contain gene IDs in each BGC region or backbone enzyme types.
#' @return Three possible outputs, depending on function arguments above:
#' @return antismash_parseR_all_orfs.csv; c("Accession", "cluster", "chrm", "region", "locus_tag", "gene", "type", "bb_type", "start", "end", "strand", "translation")
#' @return antismash_parseR_all_products.csv; c("Accession", "cluster", "chrm", "region", "product", "product_category", "backbone_locus_tag", "backbone_gene", "bb_type")
#' @return antismash_parseR_all_gbk_genes.csv; c("Accession", "chrm", "region", "locus_tag", "gene", "start_end_strand", "excluded_from_cluster", "cluster")
#' @export

antismash_parseR <- function(dir_name,
                             accession=NA,
                             output_dir=NA,
                             output_product=TRUE,
                             output_orfs=TRUE,
                             output_gbks=FALSE,
                             include_seqs=FALSE) {
  file_basename <- basename(dir_name)

  # if accession is not provided, attempt to strip it from the directory name
  # first get basename of the file
  if (is.na(accession)) {
    accession <- gsub("_antismash.*","",file_basename)
  }

  # define output dir
  if (is.na(output_dir)) {
    output_dir <- dir_name
  } else {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }

  # find json file in the specified dir
  file_path <- list.files(dir_name,pattern="regions.js",full.names = TRUE)

  # Read the file
  cx <- V8::v8()
  cx$source(file_path) # now the variable 'data' is defined in V8
  data_from_js <- cx$get("recordData")  # fetch the JS variable 'data' into R

  #parse and flatten the json record

  # list of rows where data_from_js$regions is populated
  non_empty_regions <- purrr::keep(data_from_js$regions, ~ nrow(.x) > 0)
  # Get indices of non-empty regions
  non_empty_indices <- which(purrr::map_int(data_from_js$regions, nrow) > 0)
  # Name the list with the corresponding seq_id
  names(non_empty_regions) <- data_from_js$seq_id[non_empty_indices]

  # now for each entry in non_empty_regions,
  # concatenate the tables from "orfs" together
  # add new columns showing chromosome and region, corresponding to the list index

  all_orfs <- purrr::imap_dfr(non_empty_regions, function(region_df, chrm) {
    purrr::imap_dfr(region_df$orfs, function(orf_df, orf_index) {

      if (is.null(orf_df)) {
        return(NULL)
      }

      dplyr::mutate(
        orf_df,
        chrm = chrm,
        region = orf_index
      )
    })
  })
  # add user-specified accession column
  all_orfs$Accession <- accession

  # add backbone type column
  all_orfs$bb_type <- "N/A"
  # find indexes where type=="biosynthetic"
  biosyn_idx <- which(all_orfs$type == "biosynthetic")

  # run find_bb_type on each index and save the output as the corresponding column in all_orfs
  all_orfs$bb_type[biosyn_idx] <- sapply(
    biosyn_idx,
    find_bb_type,
    all_orfs = all_orfs
  )


  ############
  # NEXT - PARSE GBK FILES TO FILL IN MISSING GENE NAMES
  gb_filenames <- list.files(dir_name,pattern=".gbk",full.names = TRUE)
  gb_filenames <- gb_filenames[grepl("region",gb_filenames)]

  ##### import_gbk for all

  all_gbk <- do.call(
    rbind,
    lapply(gb_filenames, function(f) import_gbk(gbk_file = f))
  )

  # some of the entries in the gbk are not included in the antismash cluster
  # mark them here in case you ever want to look back at them:
  all_gbk$excluded_from_cluster <- ifelse(
    all_gbk$locus_tag %in% all_orfs$locus_tag,
    "no",
    "yes"
  )

  # add unique cluster column
  all_gbk$cluster <- paste0(all_gbk$chrm,".region",sprintf("%03s", all_gbk$region))
  all_gbk <- all_gbk[order(all_gbk$cluster), ]
  all_gbk$Accession <- accession
  # order the columns
  order=c("Accession", "chrm", "region", "locus_tag", "gene", "start_end_strand", "excluded_from_cluster", "cluster")
  all_gbk2 <- all_gbk[,order]
  # save
  if (output_gbks==TRUE) {
    write.csv(all_gbk2,
              paste0(output_dir,"/",accession,"_antismash_parseR_all_gbk_genes.csv"),
              row.names=FALSE)
    print(paste0("All GBK ORFs have been parsed for: ",dir_name))
  }

  # add gene column to all_orfs
  simple <- all_gbk[,c("locus_tag","gene")]
  all_orfs_2 <- merge(all_orfs,simple,by="locus_tag")

  # add cluster column

  all_orfs_2$cluster <- paste0(all_orfs_2$chrm,".region",sprintf("%03s", all_orfs_2$region))

  # reorder columns
  if (include_seqs ==TRUE) {
    order <- c("Accession","cluster","chrm","region","locus_tag","gene","type","bb_type","start","end","strand","translation")
  } else {
    order <- c("Accession","cluster","chrm","region","locus_tag","gene","type","bb_type","start","end","strand")
  }
  all_orfs_2 <- all_orfs_2[,order]

  all_orfs_2 <- all_orfs_2[order(all_orfs_2$cluster), ]

  # ONE WEIRD THING THAT SOMETIMES HAPPENS: the dna column will parse across two separate rows and mess up the csv formatting. this is because excel has a char limit for displaying csv. so I decided to remove the dna column.

  # now save output:
  if (output_orfs==TRUE) {
    write.csv(all_orfs_2,
              paste0(output_dir,"/",accession,"_antismash_parseR_all_orfs.csv"),
              row.names=FALSE)
    print(paste0("Cluster ORFs have been parsed for: ",dir_name))
  }

  ##########
  # NEXT - parse each cluster by backbone type
  # concatenate region_df$products and region_df$product_categories
  # then merge them by orf_index
  # Create product table
  all_products <- purrr::imap_dfr(non_empty_regions, function(region_df, chrm) {

    # Number of ORFs in this region
    n_orfs <- length(region_df$orfs)

    # Safely get products and product_categories
    products <- region_df$products
    if (is.null(products)) {
      products <- vector("list", n_orfs)
    }

    product_categories <- region_df$product_categories
    if (is.null(product_categories)) {
      product_categories <- vector("list", n_orfs)
    }

    # Combine by ORF index
    purrr::map_dfr(seq_len(n_orfs), function(i) {
      tibble::tibble(
        chrm = chrm,
        region = i,
        product = paste(products[[i]], collapse = ";"),
        product_category = paste(product_categories[[i]], collapse = ";")
      )
    })
  })

  # add user-specified accession column
  all_products$Accession <- accession

  # also add two columns specifying the backbone genes
  backbones <- all_orfs_2[,c("cluster","locus_tag","gene","chrm","region","bb_type")]
  backbones <- subset(backbones,bb_type != "N/A")
  # flatten so there is only one row for each backbones$cluster. for each column, separate using semicolon space
  backbones_flat <- backbones |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      across(c(locus_tag, gene, bb_type), ~ paste(., collapse = "//")),
      .groups = "drop"
    )
  # now merge onto all_products by cluster column
  all_products$cluster <- paste0(all_products$chrm,".region",sprintf("%03s", all_products$region))

  all_products_2 <- merge(all_products,backbones_flat,by="cluster")

  # reorder and rename columns
  order <- c("Accession","cluster","chrm","region","product","product_category","locus_tag","gene","bb_type")
  all_products_2 <- all_products_2[,order]
  colnames(all_products_2)[7:9] <- c("backbone_locus_tag","backbone_gene","bb_type")

  # sort by cluster
  all_products_2 <- all_products_2[order(all_products_2$cluster), ]

  # now save output:
  if (output_product==TRUE) {
    write.csv(all_products_2,
              paste0(output_dir,"/",accession,"_antismash_parseR_all_products.csv"),
              row.names=FALSE)
    print(paste0("Cluster products have been parsed for: ",dir_name))

  }
}
