#' Function that assigns variant oncogenicity evidence based on ACMG guidelines
#'
#' @param var_calls sample calls with dbnsfp annotations
#' @param pcgr_config pcgr configuration object
#' @param pcgr_data pcgr data object
#'
#' @return var_calls
#'
#' @export
assign_oncogenicity_evidence <- function(var_calls, pcgr_config, pcgr_data) {


  annotation_cols_necessary <-
    c("VAR_ID",
      "CONSEQUENCE",
      "PROTEIN_CHANGE",
      "HGVSp_short",
      "MUTATION_HOTSPOT",
      "MUTATION_HOTSPOT_CANCERTYPE",
      "SYMBOL",
      "ONCOGENE",
      "TUMOR_SUPPRESSOR",
      "TUMOR_SUPPRESSOR_EVIDENCE",
      "ONCOGENE_EVIDENCE",
      "LOSS_OF_FUNCTION",
      "INTRON_POSITION",
      "EXON_POSITION",
      "NULL_VARIANT",
      "EAS_AF_GNOMAD",
      "NFE_AF_GNOMAD",
      "AFR_AF_GNOMAD",
      "AMR_AF_GNOMAD",
      "SAS_AF_GNOMAD")

  cols_present <- assertable::assert_colnames(
    var_calls,
    annotation_cols_necessary,
    only_colnames = F,
    quiet = T)

  acmg_ev_codes <-

    ### Beningn effects of somatic variants
    c("ACMG_SBVS1",
      ## Very high MAF: > 0.05 in gnomAD - any 5 general continental populations
      ## AFR/AMR/EAS/NFE/SAS

      "ACMG_SBS1",
      ## High MAF: > 0.01 in gnomAD - any 5 general continental populations
      ## AFR/AMR/EAS/NFE/SAS

      "ACMG_SBP1",
      ## Multiple lines of computational evidence support a benign
      ## effect on the gene or gene product
      ## (conservation, evolutionary, splicing impact, etc. - from dbNSFP

      "ACMG_SBP2",
      ## Synonymous (silent) variant for which splicing prediction
      ## algorithms predict no effect on the splice consensus sequence
      ## nor the creation of a new splice site and the nucleotide is
      ## not highly conserved

      #"ACMG_SBS2",
      ## Well established in invitro/in vivo functional studies show
      ## no oncogenic effects

      ## Pathogenic effects of somatic variants
      "ACMG_OVS1",
      ## Null variant - predicted as LoF by LOFTEE
      ## - Nonsense, frameshift, canonical splice sites, initiation codon,
      ##   single-exon/multi-exon deletion
      ## - Tumor suppressor gene

      "ACMG_OS1",
      ## Same amino acid change as previously established oncogenic variant
      ## NOTE: Not relevant since current implementation does not consider
      ## any criteria where nucleotide change is affecting annotation

      # "ACMG_OS2",
      ## Well established in invitro/in vivo functional studies show
      ## oncogenic effect of the variant

      "ACMG_OS3",
      ## Located in a mutation hotspot
      ## - >= 50 samples with a somatic variant at the AA position (cancerhotspots.org)
      ## - Same amino acid change in >= 10 samples (cancerhotspots.org)

      "ACMG_OM1",
      ## Located in a critical and well-established part of a functional domain
      ## (active site of an enzyme)

      "ACMG_OM2",
      ## Protein length changes as a result of in-frame deletions/insertions in a
      ## known oncogene/tumor suppressor genes or stop-loss variants in a
      ## tumor suppressor gene

      "ACMG_OM3",
      ## Missense variant at an amino acid residue where a different missense
      ## variant determined to be oncogenic (using this standard) has been
      ## documented. Amino acid difference from reference amino acid should
      ## be greater or at least approximately the same as for missense change
      ## determined to be oncogenic.

      "ACMG_OM4",
      ## Located in a mutation hotspot
      ## - < 50 samples with a somatic variant at the AA position (cancerhotspots.org)
      ## - Same amino acid change in >= 10 samples (cancerhotspots.org)
      ## - Not applicable if OM1 or OM3 is applicable

      "ACMG_OP1",
      ## All used lines of computational support an oncogenic effect
      ## of a variant (conservation, evolutionary, splicing effect)

      "ACMG_OP2",
      ## Somatic variant in a gene in a malignancy with a single genetic
      ## etiology. Example:retinoblastoma is caused by bi-allelic
      ## RB1 inactivation.

      "ACMG_OP3",
      ## Located in a mutation hotspot
      ## - Same amino acid change in < 10 samples (cancerhotspots.org)
      ## - Not applicable if OM1 or OM3 is applicable

      "ACMG_OP4")
      ## Absent from controls (gnomAD)
      ## - Extremely low MAF

  path_columns <-
    c(acmg_ev_codes,
      "N_INSILICO_CALLED",
      "N_INSILICO_DAMAGING",
      "N_INSILICO_TOLERATED",
      "N_INSILICO_SPLICING_NEUTRAL",
      "N_INSILICO_SPLICING_AFFECTED",
      "hotspot_symbol",
      "hotspot_codon",
      "hotspot_alt_aa",
      "hotspot_pvalue",
      "hotspot_ttype",
      "hotspot_ttype_samples_site",
      "hotspot_ttype_samples_aa")
  var_calls <- var_calls[, !(colnames(var_calls) %in% path_columns)]

  ttype_sample <- "Cancer, NOS"
  if(pcgr_config[["t_props"]][["tumor_type"]] != "Cancer, NOS"){
    ttype_sample <- stringr::str_replace(
      stringr::str_replace_all(
        pcgr_config[["t_props"]][["tumor_type"]],
        " ", "_"),
      "/", "@")
  }

  mutation_hotspots <- var_calls %>%
    dplyr::select(
      .data$VAR_ID,
      .data$SYMBOL,
      .data$CONSEQUENCE,
      .data$HGVSp_short,
      .data$MUTATION_HOTSPOT,
      .data$MUTATION_HOTSPOT_CANCERTYPE) %>%
    dplyr::filter(
      !is.na(CONSEQUENCE) &
        !is.na(HGVSp_short) &
        !is.na(.data$MUTATION_HOTSPOT) &
        !is.na(.data$MUTATION_HOTSPOT_CANCERTYPE)
    )

  if(NROW(mutation_hotspots) > 0){
    mutation_hotspots <- mutation_hotspots %>%
      tidyr::separate(
        .data$MUTATION_HOTSPOT,
        c("hotspot_symbol",
          "hotspot_codon",
          "hotspot_alt_aa",
          "hotspot_pvalue"),
        sep = "\\|", remove = T, extra = "drop") %>%
      tidyr::separate_rows(
        .data$MUTATION_HOTSPOT_CANCERTYPE, sep=","
      ) %>%
      tidyr::separate(
        .data$MUTATION_HOTSPOT_CANCERTYPE,
        c("hotspot_ttype",
          "hotspot_ttype_samples_site",
          "hotspot_ttype_samples_aa"),
        sep="\\|", remove = T, extra = "drop"
      ) %>%
      dplyr::mutate(hotspot_ttype_samples_site = as.numeric(
        hotspot_ttype_samples_site
      )) %>%
      dplyr::mutate(hotspot_ttype_samples_aa = as.numeric(
        hotspot_ttype_samples_aa
      ))

    hotspots_any_tt <- mutation_hotspots %>%
      dplyr::group_by(.data$VAR_ID,
                      .data$SYMBOL,
                      .data$CONSEQUENCE,
                      .data$hotspot_symbol,
                      .data$hotspot_codon,
                      .data$hotspot_alt_aa) %>%
      dplyr::summarise(hotspot_ttype_samples_site = sum(
        hotspot_ttype_samples_site),
        hotspot_ttype_samples_aa = sum(
          hotspot_ttype_samples_aa
        ),
        .groups = "drop"
      )




  }



  ## Assign logical ACMG evidence indicators
  #
  #
  # ACMG_OP1 - Multiple lines (>=7) of insilico evidence support a
  #             deleterious effect on the gene or gene product
  ##           (conservation, evolutionary, splicing impact, etc.)
  # ACMG_SBP1 - Multiple lines (>=7) of insilico evidence support a benign effect.
  #
  # Computational evidence for deleterious/benign effect is taken from
  # invidual algorithm predictions in dbNSFP: SIFT,Provean,MutationTaster,
  # MutationAssessor,M_CAP,MutPred,FATHMM,FATHMM-mkl,DBNSFP_RNN,dbscSNV_RF,
  # dbscSNV_AdaBoost
  # 1) Damaging: Among all possible protein variant effect predictions, at
  #              least six algorithms must have made a call,
  #              with at least 5 predicted as damaging/D
  #              (possibly_damaging/PD), and at most two
  #              predicted as tolerated/T (OP1)
  #       - at most 1 prediction for a splicing neutral effect
  #    Exception: if both splice site predictions indicate damaging effects;
  #    ignore other criteria
  # 2) Tolerated: Among all possible protein variant effect predictions, at
  #    least six algorithms must have made a call,
  #    with at least 5 predicted as tolerated, and at most two
  #    predicted as damaging (SBP1)
  #    - 0 predictions of splice site affected

  dbnsfp_min_majority <- 7
  dbnsfp_max_minority <- 2
  dbnsfp_min_called <- dbnsfp_min_majority

  var_calls <- cpsr::get_insilico_prediction_statistics(var_calls)

  ## Assign logical ACMG evidence indicators based on computational
  ## evidence that support damaging/benign effect
  var_calls <- var_calls %>%
    dplyr::mutate(
      ACMG_OP1 =
        dplyr::if_else(
          .data$N_INSILICO_CALLED >= dbnsfp_min_called &
            .data$N_INSILICO_DAMAGING >= dbnsfp_min_majority &
            .data$N_INSILICO_TOLERATED <= dbnsfp_max_minority &
            .data$N_INSILICO_SPLICING_NEUTRAL <= 1, TRUE,
          FALSE, FALSE)) %>%
    dplyr::mutate(
      ACMG_SBP1 = dplyr::if_else(
        .data$N_INSILICO_CALLED >= dbnsfp_min_called &
          .data$N_INSILICO_TOLERATED >= dbnsfp_min_majority &
          .data$N_INSILICO_DAMAGING <= dbnsfp_max_minority &
          .data$N_INSILICO_SPLICING_AFFECTED == 0, TRUE,
        FALSE, FALSE)) %>%
    dplyr::mutate(ACMG_OP1 = dplyr::case_when(
      .data$N_INSILICO_SPLICING_AFFECTED == 2 ~ TRUE,
      TRUE ~ as.logical(.data$ACMG_OP1)))

  ## Assign logical ACMG evidence indicators based on population frequency
  ## data in non-cancer samples from gnomAD (Dominant vs. recessive
  ## modes of inheritance)
  # 'ACMG_SBSV1'   -  Very high MAF (> 0.5% in any gnomAD pop subsets) -
  # 'ACMG_SBS1' -     High MAF (> 0.1% in any gnomAD pop subset) -

  if ("EAS_AF_GNOMAD" %in% colnames(var_calls) &
      "SAS_AF_GNOMAD" %in% colnames(var_calls) &
      "NFE_AF_GNOMAD" %in% colnames(var_calls) &
      "AMR_AF_GNOMAD" %in% colnames(var_calls) &
      "AFR_AF_GNOMAD" %in% colnames(var_calls)) {

    var_calls <- var_calls %>%
      dplyr::mutate(ACMG_SBVS1 = dplyr::if_else(
        (!is.na(EAS_AF_GNOMAD) & EAS_AF_GNOMAD > 0.05) |
          (!is.na(SAS_AF_GNOMAD) & SAS_AF_GNOMAD > 0.05) |
          (!is.na(NFE_AF_GNOMAD) & NFE_AF_GNOMAD > 0.05) |
          (!is.na(AMR_AF_GNOMAD) & AMR_AF_GNOMAD > 0.05) |
          (!is.na(AFR_AF_GNOMAD) & AFR_AF_GNOMAD > 0.05),
        TRUE,
        FALSE
      )) %>%
      dplyr::mutate(ACMG_SBS1 = dplyr::if_else(
        ((!is.na(EAS_AF_GNOMAD) & EAS_AF_GNOMAD > 0.01) |
          (!is.na(SAS_AF_GNOMAD) & SAS_AF_GNOMAD > 0.01) |
          (!is.na(NFE_AF_GNOMAD) & NFE_AF_GNOMAD > 0.01) |
          (!is.na(AMR_AF_GNOMAD) & AMR_AF_GNOMAD > 0.01) |
          (!is.na(AFR_AF_GNOMAD) & AFR_AF_GNOMAD > 0.01)) &
          ACMG_SBVS1 == FALSE,
        TRUE,
        FALSE
      ))
  }


  ## Assign logical ACMG evidence indicators
  #
  # 'ACMG_OVS1' - Null variant (frameshift, nonsense) -
  #               predicted as LoF by LOFTEE -
  #               in tumor suppressor gene
  # 'ACMG_OM2 - Protein length changes as a result of
  #             in-frame deletions/insertions in a known
  #             oncogene/tumor suppressor genes or stop-loss
  #             variants in atumor suppressor gene

  var_calls <- var_calls %>%
    dplyr::mutate(
      ACMG_OVS1 =
        dplyr::if_else(
          .data$NULL_VARIANT == T &
            .data$LOSS_OF_FUNCTION == T &
            .data$TUMOR_SUPPRESSOR == T,
          TRUE, FALSE, FALSE))

  var_calls <- var_calls %>%
    dplyr::mutate(
      ACMG_OM2 =
        dplyr::if_else(
          (stringr::str_detect(
            .data$CONSEQUENCE,
            "^(inframe_deletion|inframe_insertion)") &
            (.data$TUMOR_SUPPRESSOR == T |
               .data$ONCOGENE == T)) |
            (stringr::str_detect(
              .data$CONSEQUENCE,
              "^stop_lost"
            ) &
              .data$TUMOR_SUPPRESSOR == T),
          TRUE, FALSE, FALSE)
    )


  ## Assign a logical ACMG evidence indicator
  # ACMG_SBP2 - Silent/intronic variant outside of the splice site consensus
  var_calls <- var_calls %>%
    dplyr::mutate(
      ACMG_SBP2 =
        dplyr::if_else((
          (as.integer(.data$INTRON_POSITION) < 0 & as.integer(.data$INTRON_POSITION) < -3) |
            (as.integer(.data$INTRON_POSITION) > 0 & as.integer(.data$INTRON_POSITION) > 6) |
            (as.integer(.data$EXON_POSITION) < 0 & as.integer(.data$EXON_POSITION) < -2) |
            (as.integer(.data$EXON_POSITION) > 0 & as.integer(.data$EXON_POSITION) > 1)) &
            stringr::str_detect(
              .data$CONSEQUENCE,
              "^(synonymous_variant|intron_variant|splice_region_variant)"),
          TRUE, FALSE, FALSE))

  ## Assign a logical ACMG evidence indicator
  # ACMG_BP3 - Variants in promoter or untranslated regions
  var_calls <- var_calls %>%
    dplyr::mutate(
      ACMG_BP3 =
        dplyr::if_else(
          stringr::str_detect(
            .data$CONSEQUENCE,
            "^(downstream|upstream|5_prime_UTR_variant|3_prime_UTR_variant)"),
          TRUE, FALSE, FALSE))


  ## Assign logical ACMG evidence indicators
  # ACMG_PS1 - coinciding with known pathogenic missense variants
  # (yet with different nucleotide change)
  # ACMG_PM5 - occurs at the same codon as a known pathogenic missense variant
  # ACMG_BSC1 - coinciding with known benign missense variants
  # ACMG_BMC1 - occurs at the same codon as a known benign missense variant

  var_calls$codon_prefix <- NA
  if (nrow(var_calls[!is.na(var_calls$CONSEQUENCE) &
                     var_calls$CONSEQUENCE == "missense_variant", ]) > 0) {
    var_calls[var_calls$CONSEQUENCE == "missense_variant", ]$codon_prefix <-
      stringr::str_match(var_calls[var_calls$CONSEQUENCE == "missense_variant", ]$HGVSp_short, "p\\.[A-Z]{1}[0-9]{1,}")
  }

  if (nrow(var_calls[!is.na(var_calls$codon_prefix), ]) > 0) {
    var_calls_pathogenic_codon <-
      dplyr::left_join(
        dplyr::filter(
          dplyr::select(
            var_calls, .data$VAR_ID,
            .data$codon_prefix, .data$SYMBOL),
          !is.na(.data$codon_prefix)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["P"]][["codon"]],
        by = c("codon_prefix" = "codon_prefix", "SYMBOL" = "symbol")) %>%
      dplyr::filter(.data$clinvar_pathogenic_codon == T)
    var_calls <- var_calls %>%
      dplyr::left_join(
        dplyr::select(
          var_calls_pathogenic_codon,
          .data$VAR_ID, .data$clinvar_pathogenic_codon),
        by = c("VAR_ID"))

    var_calls_benign_codon <-
      dplyr::left_join(
        dplyr::filter(dplyr::select(var_calls, .data$VAR_ID, .data$codon_prefix, .data$SYMBOL),
                      !is.na(.data$codon_prefix)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["B"]][["codon"]],
        by = c("codon_prefix" = "codon_prefix", "SYMBOL" = "symbol")) %>%
      dplyr::filter(.data$clinvar_benign_codon == T)

    var_calls <- var_calls %>%
      dplyr::left_join(
        dplyr::select(var_calls_benign_codon, .data$VAR_ID, .data$clinvar_benign_codon),
        by = c("VAR_ID")) %>%
      dplyr::mutate(ACMG_PM5 = dplyr::if_else(
        .data$clinvar_pathogenic_codon == TRUE, TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(ACMG_BMC1 = dplyr::if_else(
        .data$clinvar_benign_codon == TRUE, TRUE, FALSE, FALSE))
  }else{
    var_calls$ACMG_PM5 <- FALSE
    var_calls$ACMG_BMC1 <- FALSE
    var_calls$clinvar_pathogenic_codon <- NA
    var_calls$clinvar_benign_codon <- NA
  }


  if (nrow(var_calls[!is.na(var_calls$HGVSp_short), ]) > 0) {
    var_calls_pathogenic_hgvsp <-
      dplyr::left_join(
        dplyr::filter(
          dplyr::select(var_calls, .data$VAR_ID, .data$HGVSp_short, .data$SYMBOL),
          !is.na(.data$HGVSp_short)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["P"]][["peptide_change"]],
        by = c("HGVSp_short" = "hgvs_p", "SYMBOL" = "symbol")) %>%
      dplyr::filter(.data$clinvar_pathogenic == T)

    var_calls <- var_calls %>%
      dplyr::left_join(
        dplyr::select(var_calls_pathogenic_hgvsp,
                      .data$VAR_ID, .data$clinvar_pathogenic), by = c("VAR_ID"))

    var_calls_benign_hgvsp <-
      dplyr::left_join(
        dplyr::filter(
          dplyr::select(var_calls, .data$VAR_ID, .data$HGVSp_short, .data$SYMBOL),
          !is.na(.data$HGVSp_short)),
        pcgr_data[["clinvar"]][["cpg_loci"]][["hi"]][["B"]][["peptide_change"]],
        by = c("HGVSp_short" = "hgvs_p", "SYMBOL" = "symbol")) %>%
      dplyr::filter(.data$clinvar_benign == T)

    var_calls <- var_calls %>%
      dplyr::left_join(
        dplyr::select(var_calls_benign_hgvsp, .data$VAR_ID, .data$clinvar_benign),
        by = c("VAR_ID")) %>%
      dplyr::mutate(ACMG_PS1 =
                      dplyr::if_else(.data$clinvar_pathogenic == TRUE,
                                     TRUE, FALSE, FALSE)) %>%
      dplyr::mutate(ACMG_BSC1 =
                      dplyr::if_else(.data$clinvar_benign == TRUE,
                                     TRUE, FALSE, FALSE))
  }else{
    var_calls$ACMG_PS1 <- FALSE
    var_calls$ACMG_BSC1 <- FALSE
    var_calls$clinvar_pathogenic <- NA
    var_calls$clinvar_benign <- NA
  }

  ## if previously found coinciding with pathogenic variant (ACMG_PS1),
  # set ACMG_PM5 to false
  var_calls <- var_calls %>%
    dplyr::mutate(
      ACMG_PM5 =
        dplyr::case_when(.data$ACMG_PM5 == T & .data$ACMG_PS1 == T ~ FALSE,
                         TRUE ~ as.logical(.data$ACMG_PM5))) %>%
    ## if previously found coinciding with benign variant (ACMG_BSC1),
    ##  set ACMG_BMC1 to false
    dplyr::mutate(
      ACMG_BMC1 =
        dplyr::case_when(.data$ACMG_BMC1 == T & .data$ACMG_BSC1 == T ~ FALSE,
                         TRUE ~ as.logical(.data$ACMG_BMC1))) %>%
    dplyr::select(
      -c(.data$clinvar_pathogenic_codon, .data$clinvar_benign_codon,
         .data$clinvar_pathogenic,
         .data$clinvar_benign, .data$cpsr_gene_moi, .data$gad_af))

  ##Assign logical ACMG level
  # PM1 - missense variant in a somatic mutation hotspot as
  # determined by cancerhotspots.org (v2)
  var_calls <- var_calls %>%
    tidyr::separate(
      .data$MUTATION_HOTSPOT,
      c("hotspot_symbol", "hotspot_codon",
        "hotspot_alt_aa", "hotspot_pvalue"),
      sep = "\\|", remove = F, extra = "drop") %>%
    dplyr::mutate(
      hotspot_codon =
        dplyr::if_else(
          !is.na(.data$hotspot_codon),
          paste0("p.", .data$hotspot_codon),
          as.character(NA))) %>%
    dplyr::mutate(
      ACMG_PM1 =
        dplyr::if_else(!is.na(.data$hotspot_codon) &
                         !is.na(.data$hotspot_symbol) &
                         !is.na(.data$codon_prefix) &
                         .data$SYMBOL == .data$hotspot_symbol &
                         .data$hotspot_codon == .data$codon_prefix,
                       TRUE, FALSE))

  return(var_calls)
}


#' Function that assigns final oncogenicity classification (B, LB, VUS, P, LP)
#' based on accumulated scores from different ACMG criteria and pre-defined
#' cutoffs
#'
#' @param var_calls data frame with variant calls
#'
#' @return var_calls data frame with oncogenicity classification appended
#'
#' @export
determine_oncogenicity_classification <- function(var_calls) {

  evidence_codes <- pcgrr::pcgr_acmg_oncogenicity[["evidence_codes"]]

  path_cols <- c("ONCOGENICITY_CLASSIFICATION", "ONCOGENICITY_CLASSIFICATION_DOC",
                 "ONCOGENICITY_CLASSIFICATION_CODE",
                 "pcgr_score_pathogenic", "pcgr_score_benign")
  var_calls <- var_calls[, !(colnames(var_calls) %in% path_cols)]

  var_calls$ONCOGENICITY_CLASSIFICATION <- "VUS"
  var_calls$ONCOGENICITY_CLASSIFICATION_DOC <- ""
  var_calls$ONCOGENICITY_CLASSIFICATION_CODE <- ""
  var_calls$pcgr_score_pathogenic <- 0
  var_calls$pcgr_score_benign <- 0

  i <- 1
  while (i <= nrow(evidence_codes)) {
    category <- evidence_codes[i, ]$category
    pole <- evidence_codes[i, ]$oncogenicity_pole
    description <- evidence_codes[i, ]$description
    pcgr_evidence_code <- evidence_codes[i, ]$pcgr_evidence_code
    score <- evidence_codes[i, ]$path_score
    if (pcgr_evidence_code %in% colnames(var_calls)) {
      var_calls <- var_calls %>%
        dplyr::mutate(
          pcgr_score_benign = .data$pcgr_score_benign +
            dplyr::if_else(pole == "B" & !!rlang::sym(pcgr_evidence_code) == T,
                           score, 0)) %>%
        dplyr::mutate(
          pcgr_score_pathogenic = .data$pcgr_score_pathogenic +
            dplyr::if_else(pole == "P" & !!rlang::sym(pcgr_evidence_code) == T,
                           score, 0)) %>%
        dplyr::mutate(
          ONCOGENICITY_CLASSIFICATION_DOC =
            paste0(.data$ONCOGENICITY_CLASSIFICATION_DOC,
                   dplyr::if_else(!!rlang::sym(pcgr_evidence_code) == T,
                                  paste0("- ", description), ""),
                   sep = "<br>")) %>%
        dplyr::mutate(
          ONCOGENICITY_CLASSIFICATION_CODE =
            paste0(.data$ONCOGENICITY_CLASSIFICATION_CODE,
                   dplyr::if_else(!!rlang::sym(pcgr_evidence_code) == T,
                                  pcgr_evidence_code, ""), sep = "|"))
    }
    i <- i + 1
  }

  lb_upper_limit <- -1
  lb_lower_limit <- -6
  b_upper_limit <- -7
  vus_lower_limit <- 0
  vus_upper_limit <- 5
  lp_lower_limit <- 6
  lp_upper_limit <- 9
  p_lower_limit <- 10


  var_calls <- var_calls %>%
    dplyr::mutate(
      ONCOGENICITY_CLASSIFICATION_CODE =
        stringr::str_replace_all(
          stringr::str_replace_all(
            .data$ONCOGENICITY_CLASSIFICATION_CODE,
            "(\\|{2,})", "|"),
          "(^\\|)|(\\|$)", "")
    ) %>%
    dplyr::mutate(
      ONCOGENICITY_CLASSIFICATION_DOC =
        stringr::str_replace_all(
          stringr::str_replace_all(
            .data$ONCOGENICITY_CLASSIFICATION_DOC,
            "(<br>){2,}", "<br>"), "(^(<br>))|((<br>)$)", "")) %>%

    ## Adjust scores in cases where critera are acting as a
    ## prerequisite for other criteria
    dplyr::mutate(
      pcgr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(
            .data$ONCOGENICITY_CLASSIFICATION_CODE, "ACMG_PVS") &
            stringr::str_detect(
              .data$ONCOGENICITY_CLASSIFICATION_CODE, "ACMG_PM2_2"),
          .data$pcgr_score_pathogenic - 1, .data$pcgr_score_pathogenic)) %>%
    dplyr::mutate(
      pcgr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(
            .data$ONCOGENICITY_CLASSIFICATION_CODE, "ACMG_PVS") &
            stringr::str_detect(
              .data$ONCOGENICITY_CLASSIFICATION_CODE, "ACMG_PM2_1"),
          .data$pcgr_score_pathogenic - 0.5, .data$pcgr_score_pathogenic)) %>%
    dplyr::mutate(
      pcgr_score_pathogenic =
        dplyr::if_else(
          stringr::str_detect(
            .data$ONCOGENICITY_CLASSIFICATION_CODE, "ACMG_OVS1_10") &
            stringr::str_detect(
              .data$ONCOGENICITY_CLASSIFICATION_CODE, "ACMG_PP3"),
          .data$pcgr_score_pathogenic - 0.5, .data$pcgr_score_pathogenic)) %>%

    ## Add scores accumulated with benign criteria and pathogenic criteria
    dplyr::mutate(
      ONCOGENICITY_SCORE =
        dplyr::if_else(
          .data$pcgr_score_benign == 0,
          .data$pcgr_score_pathogenic,
          .data$pcgr_score_benign)) %>%
    dplyr::mutate(
      ONCOGENICITY_SCORE =
        dplyr::if_else(
          .data$pcgr_score_benign < 0 &
            .data$pcgr_score_pathogenic > 0,
          .data$pcgr_score_benign + .data$pcgr_score_pathogenic,
          .data$ONCOGENICITY_SCORE)) %>%
    dplyr::mutate(
      ONCOGENICITY_CLASSIFICATION =
        dplyr::case_when(
          .data$ONCOGENICITY_SCORE <= lb_upper_limit &
            .data$ONCOGENICITY_SCORE >= lb_lower_limit ~ "Likely_Benign",
          .data$ONCOGENICITY_SCORE <= b_upper_limit ~ "Benign",
          .data$ONCOGENICITY_SCORE <= vus_upper_limit &
            .data$ONCOGENICITY_SCORE >= vus_lower_limit ~ "VUS",
          .data$ONCOGENICITY_SCORE >= p_lower_limit ~ "Pathogenic",
          .data$ONCOGENICITY_SCORE >= lp_lower_limit &
            .data$ONCOGENICITY_SCORE <= lp_upper_limit ~ "Likely_Pathogenic",
          TRUE ~ as.character("VUS"))) %>%
    dplyr::select(-c(.data$pcgr_score_benign, .data$pcgr_score_pathogenic))

  return(var_calls)

}
