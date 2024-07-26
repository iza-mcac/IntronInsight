get_sequence <- function(chromosome, start, end, strand) {
  if (strand == "-") {
    seq <- getSeq(genome, chromosome, start, end)
    seq <- reverseComplement(DNAString(seq))
  } else {
    seq <- getSeq(genome, chromosome, start, end)
  }
  return(as.character(seq))
}

# Get possible ESEs and ESSs
get_adjacent_exon_sequences <- function(chromosome, start, end, strand, adjacent_distance = 10) {
  upstream_start <- start - adjacent_distance
  downstream_end <- end + adjacent_distance

  if (strand == "+") {
    upstream_sequence <- get_sequence(chromosome, upstream_start, start - 1, strand)
    downstream_sequence <- get_sequence(chromosome, end + 1, downstream_end, strand)
  } else {
    downstream_sequence <- reverseComplement(DNAString(get_sequence(chromosome, upstream_start, start - 1, strand)))
    upstream_sequence <- reverseComplement(DNAString(get_sequence(chromosome, end + 1, downstream_end, strand)))
  }
  return(list(upstream_sequence = upstream_sequence, downstream_sequence = downstream_sequence))
}

#' Calculate Intron Metrics
#'
#' This function calculates various metrics for a given intron sequence. It can operate in two modes:
#' "coordinate" mode, where the user provides chromosome coordinates to fetch the sequence, and
#' "sequence" mode, where the user directly provides the intron sequence. If it is using the "sequence" mode ESE
#' and ESS extraction do not work. This also requires a genome that the user selected before everything.
#'
#' @param mode Character string specifying the mode of operation. Use "coordinate" to fetch the sequence
#' from the genome using chromosome coordinates, or "sequence" to provide the sequence directly.
#' Default is "coordinate".
#' @param chromosome Character string specifying the chromosome (e.g., "chr1"). Required if mode is "coordinate".
#' @param start Integer specifying the start position of the intron on the chromosome. Required if mode is "coordinate".
#' @param end Integer specifying the end position of the intron on the chromosome. Required if mode is "coordinate".
#' @param strand Character string specifying the strand ("+" or "-"). Required in both modes.
#' @param sequence Character string specifying the intron sequence. Required if mode is "sequence".
#'
#' @return A data frame containing the calculated metrics for the intron:
#' \item{chromosome}{Chromosome name}
#' \item{start}{Start position of the intron}
#' \item{end}{End position of the intron}
#' \item{strand}{Strand information}
#' \item{intron_length}{Length of the intron}
#' \item{percentage_cg}{Percentage of CG content in the intron}
#' \item{count_polyP}{Count of polyP (C and T) tracts}
#' \item{polyp_sequence}{Sequence of the polyP tract}
#' \item{intron_start_motif}{First 10 nucleotides of the intron}
#' \item{intron_end_motif}{Last 10 nucleotides of the intron}
#' \item{ESS_region}{ESS region sequence (only for "coordinate" mode)}
#' \item{ESE_region}{ESE region sequence (only for "coordinate" mode)}
#' \item{sequence}{Full intron sequence}
#' @export
#'
#' @examples
#' #' # Example using chromosome coordinates
#' calculate_intron_metrics(mode = "coordinate", chromosome = "chr1", start = 100000, end = 100100, strand = "-")
#'
#' # Example using direct sequence input
#' calculate_intron_metrics(mode = "sequence", sequence = "AGCTTGCGATCGTAGCTAGCTAGCTGACTGACTGACTGACTGACTGACTGAC", strand = "-")
#'
calculate_intron_metrics <- function(mode = "coordinate", chromosome = NULL, start = NULL, end = NULL, strand = NULL, sequence = NULL) {
  if (mode == "sequence") {
    if (is.null(sequence) || is.null(strand)) {
      stop("For sequence mode, provide the sequence and strand.")
    }
    # User provided the sequence directly
    intron_sequence <- as.character(sequence)
  } else if (mode == "coordinate") {
    if (is.null(chromosome) || is.null(start) || is.null(end) || is.null(strand)) {
      stop("For coordinate mode, provide chromosome, start, end, and strand.")
    }
    # Use chromosome coordinates to fetch the sequence
    intron_sequence <- as.character(DNAString(get_sequence(chromosome, start, end, strand)))
  } else {
    stop("Invalid mode. Use 'coordinate' or 'sequence'.")
  }

  intron_length <- nchar(intron_sequence)

  # Percent CG
  c_count <- str_count(intron_sequence, "C")
  g_count <- str_count(intron_sequence, "G")
  cg_count <- c_count + g_count
  total_bases <- nchar(intron_sequence)
  percentage_cg <- (cg_count / total_bases) * 100

  # Extract maximum polyP tract
  polyp_sequence <- substr(intron_sequence, max(1, total_bases - 59), total_bases)
  # Count the Cs and Ts
  c_count_polyP <- str_count(polyp_sequence, "C")
  t_count_polyP <- str_count(polyp_sequence, "T")
  count_polyP <- c_count_polyP + t_count_polyP

  # Get the first and last 10 nucleotides of the intron to check for GU AG
  intron_start_motif <- substr(intron_sequence, 1, min(10, total_bases))
  intron_end_motif <- substr(intron_sequence, max(1, total_bases - 9), total_bases)

  if (mode == "coordinate") {
    adjacent_sequences <- get_adjacent_exon_sequences(chromosome, start, end, strand)
    ess_region <- as.character(adjacent_sequences$upstream_sequence)
    ese_region <- as.character(adjacent_sequences$downstream_sequence)
  } else {
    ess_region <- NA
    ese_region <- NA
  }

  result_df <- data.frame(
    chromosome = ifelse(mode == "coordinate", chromosome, NA),
    start = ifelse(mode == "coordinate", start, NA),
    end = ifelse(mode == "coordinate", end, NA),
    strand = strand,
    intron_length = intron_length,
    percentage_cg = percentage_cg,
    count_polyP = count_polyP,
    polyp_sequence = as.character(polyp_sequence),
    intron_start_motif = as.character(intron_start_motif),
    intron_end_motif = as.character(intron_end_motif),
    ESS_region = ess_region,
    ESE_region = ese_region,
    sequence = as.character(intron_sequence),
    stringsAsFactors = FALSE
  )

  return(result_df)
}
