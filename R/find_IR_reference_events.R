#'Process GFF File and Identify Unique Exons with Wiggle Room
#'
#' This function processes a GFF file to filter the species exon data, joins it with transcript-to-gene mapping data,
#' and identifies the exons to be used for comparison considering a user-defined wiggle room
#' which depends on annotation. Sugessted wiggle room for human data is 10. This takes considerable time to run
#' which depends on the size of your gff file, Should be ran only once per analysis
#'
#' @param gff_file A string specifying the path to the GFF file. The file should match the species genome version used in other steps.
#' @param fa_file A string specifying the path to the transcript-to-gene mapping file in FASTA format should match the verison of the gff..
#' @param wiggle An integer specifying the wiggle room in nucleotides for comparing exon positions. Default is 10.
#'
#' @return A data frame containing unique exons identified in the retained intron transcripts, considering the specified wiggle room.
#' @export
#'
#' @examples
#' gff_file <- "ref/gencode.v38.chr_patch_hapl_scaff.annotation.gff3.gz"
#' fa_file <- "ref/gencode.v38.transcripts.fa.gz"
#' wiggle <- 10
#' unique_exons <- process_gff_and_find_unique_exons(gff_file, fa_file, wiggle)
#' print(unique_exons)
find_IR_reference_events <- function(gff_file, fa_file, wiggle = 10) {
  #Use isoformic maketxtogene
  txtogene_v38 <- isoformic::make_tx_to_gene(file_path = fa_file)
  #Read the gff
  gencode_annotation <- read_delim(gff_file, delim = "\t", escape_double = FALSE,
                                   col_names = FALSE, trim_ws = TRUE, skip = 7)
  #Process the gff
  gencode_annotation_tomod <- gencode_annotation
  gencode_annotation_tomod[1:12] <- gencode_annotation_tomod$X9 %>% str_split_fixed(';', 12)
  gff_data <- cbind(gencode_annotation, gencode_annotation_tomod)
  names(gff_data) <- c("chr", "HAVANA", "type", "left_pos", "right_pos", "exclude_1",
                       "strand", "exclude_2", "etc", "ID", "V1", "V2", "V3", "V4",
                       "V5", "V6", "V7", "V8", "V9", "V10", "V11")
  gff_filtered <- gff_data %>% dplyr::select(chr, type, left_pos, right_pos, strand, ID)
  gff_filtered$ID <- str_sub(gff_filtered$ID, start = 4) %>% str_replace("exon:", "")

  df_split <- tidyr::separate(gff_filtered, ID, into = c("tx_id", "exon_count"), sep = ":")
  txtogene_v38 <- txtogene_v38 %>% dplyr::select(transcript_id, transcript_name, gene_name, transcript_type)
  df_split_notgene_dic <- df_split %>%
    dplyr::left_join(txtogene_v38, by = c("tx_id" = "transcript_id"))

  #Check wiggle room
  within_wiggle_room <- function(pos1, pos2, wiggle) {
    abs(pos1 - pos2) <= wiggle
  }

  #Find unique exons
  find_unique_exons <- function(df, wiggle) {
    retained_exons <- df %>% filter(transcript_type == "retained_intron", type == "exon")
    other_exons <- df %>% filter(transcript_type != "retained_intron", type == "exon")

    unique_exons <- retained_exons %>%
      rowwise() %>%
      filter(!any(
        within_wiggle_room(left_pos, other_exons$left_pos, wiggle) &
          within_wiggle_room(right_pos, other_exons$right_pos, wiggle)
      )) %>%
      ungroup()

    return(unique_exons)
  }

  # Apply function to find unique exons
  unique_exons <- df_split_notgene_dic %>%
    group_by(gene_name) %>%
    group_modify(~ find_unique_exons(.x, wiggle = wiggle)) %>%
    ungroup()

  return(unique_exons)
}



