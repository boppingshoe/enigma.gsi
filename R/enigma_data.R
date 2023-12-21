
#' Prepping isotope-genetics BMM input data2
#'
#' @description
#' Prep input data for Bayesian mixture model (i.e., Pella-Madusa model) enhanced with isotope information.
#'
#' @param mixture_data Mixture data in GCL or \pkg{rubias} format.
#' @param baseline_data Baseline data in GCL or \pkg{rubias} format.
#' @param pop_info Population information for the baseline. A tibble with columns
#'   collection (collection names), repunit (reporting unit names),
#'   grpvec (group numbers), origin (wild/hatchery).
#' @param isobreaks Indicate where to cut and put isoscape values into bins, or how many bins the isoscape values should be grouped. Recommend manually set where to cut.
#' @param isoscape Option to provide isoscape if there are no isotope values (`sr_val`) for the baseline. Isoscape is a tibble with three columns:
#'  * `collection` - collection names,
#'  * `sr_mean` - mean Sr readings for each population of the genetic baseline,
#'  * `sr_sd` - standard deviation for Sr readings for each population of the genetic baseline.
#' @param file Where you want to save a copy of input data as a RDS file.
#'   Need to type out full path and extension `.Rds`.
#'   Leave it empty if you don't want to save a copy.
#' @param loci Optional. Provide loci for the mixture or baseline as a fail-safe check.
#'
#' @return A list objects as the input data for enigma_mdl()
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' # prep input data
#' enigma_data <- prep_enigma_data(mixture_data = mix_iso, baseline_data = baseline, pop_info = ayk_pops60)
#'
#' @export
prep_enigma_data <-
  function(mixture_data, baseline_data, pop_info, isobreaks = 5, isoscape = NULL, file = NULL, loci = NULL) {

    start_time <- Sys.time()

    # identify loci for each stage
    # make sure no colnames other than marker names have ".1" at the end
    loci_base <-
      dplyr::tibble(locus = names(baseline_data)) %>%
      dplyr::filter(grepl("\\.1$", locus)) %>%
      dplyr::mutate(locus = substr(locus, 1, nchar(locus) - 2)) %>%
      dplyr::pull(locus)

    # input error check

    loci_mix <-
      dplyr::tibble(locus = names(mixture_data)) %>%
      dplyr::filter(grepl("\\.1$", locus)) %>%
      dplyr::mutate(locus = substr(locus, 1, nchar(locus) - 2)) %>%
      dplyr::pull(locus)

    error_message <- check_loci_pops(loci, loci_base, loci_mix)

    if ("all good" %in% error_message) {
      loci <- loci_base
      message("Compiling input data, may take a minute or two...")
    } else {
      stop(error_message)
    }

    # change column name if the data are gcl objects
    # to match rubias input data name convention
    if ("SILLY_CODE" %in% names(baseline_data))
      baseline_data <- dplyr::rename(baseline_data, collection = SILLY_CODE)

    if ("SillySource" %in% names(mixture_data))
      mixture_data <- dplyr::rename(mixture_data, indiv = SillySource)

    if(!any(grepl("sr_val", names(baseline_data)))) {
      if(is.null(isoscape)) stop("Need to provide isoscape or have sr_val in baseline.")

      baseline_data <- baseline_data %>%
        dplyr::mutate(
          iden = {sapply(collection, function(i) which(i == isoscape$collection))},
          sr_val = {sapply(iden, function(i) rnorm(1, isoscape$sr_mean[i], isoscape$sr_sd[i]))}
        )
    }

    # tally allele for each baseline and mixture sample
    base <- allefreq(baseline_data, baseline_data, loci, collect_by = collection) %>%
      dplyr::right_join(pop_info, by = c("collection" = "collection"), keep = FALSE) %>%
      dplyr::relocate(!dplyr::ends_with(as.character(0:9)), .after = collection) %>%
      dplyr::mutate(dplyr::across(dplyr::ends_with(as.character(0:9)), ~tidyr::replace_na(., 0))) %>%
      dplyr::left_join({
        dplyr::select(baseline_data, collection, sr_val) %>%
          dplyr::mutate(sr_bin = cut(sr_val, breaks = isobreaks)) %>%
          with(., table(collection, sr_bin)) %>%
          tibble::as_tibble() %>%
          tidyr::pivot_wider(names_from = "sr_bin", values_from = "n")
      }, by = "collection")

    if (is.null(isoscape)) {
      isoscape <-
        dplyr::select(baseline_data, collection, sr_val) %>%
        dplyr::group_by(collection) %>%
        dplyr::summarise(sr_mean = mean(sr_val),
                         sr_sd = sd(sr_val))
    }

    isoscape_ordered <- dplyr::select(base, collection) %>%
      dplyr::left_join(isoscape, by = "collection")

    mix <- allefreq(mixture_data, baseline_data, loci)

    mix <- mix %>%
      dplyr::bind_cols(dplyr::select(mixture_data, sr_val)) %>%
      dplyr::left_join({
        dplyr::select(mixture_data, indiv, sr_val) %>%
          dplyr::mutate(sr_bin = cut(sr_val, breaks = isobreaks)) %>%
          with(., table(indiv, sr_bin)) %>%
          tibble::as_tibble() %>%
          tidyr::pivot_wider(names_from = "sr_bin", values_from = "n")
      }, by = "indiv")

    # numbers of allele types
    nalleles <- lapply(loci, function(loc) {
      dplyr::tibble(locus = loc,
                    call = baseline_data %>%
                      dplyr::select(dplyr::all_of(loc), paste0(loc, ".1")) %>%
                      unlist() %>% unique() %>% .[!is.na(.)],
                    altyp = seq.int(dplyr::n_distinct(call)) %>% factor())
    }) %>% dplyr::bind_rows() %>%
      dplyr::group_by(locus) %>%
      dplyr::summarise(n_allele = max(as.numeric(altyp)), .groups = "drop")

    n_alleles <- nalleles %>%
      dplyr::pull(n_allele) %>%
      stats::setNames(nalleles$locus)

    n_alleles["isoscape"] <- length(isobreaks) - 1

    # group names
    grp_nms <- base %>%
      dplyr::arrange(grpvec) %>%
      dplyr::pull(repunit) %>%
      unique()

    # wild or hatchery
    if ("origin" %in% names(base)) {
      wildpops <- base %>%
        dplyr::filter(origin == "wild") %>%
        dplyr::pull(collection)
      hatcheries <- base %>%
        dplyr::filter(origin == "hatchery") %>%
        dplyr::pull(collection)
    } else {
      wildpops <- base %>% dplyr::pull(collection)
      hatcheries <- NULL
    }

    # iden if specified in mixture data
    if (any(grepl("known_", names(mixture_data)))) {
      iden <- mixture_data %>%
        dplyr::select(tidyr::contains("known_")) %>%
        dplyr::pull()
      if (!all(stats::na.omit(iden) %in% c(wildpops, hatcheries))) {
        stop(c("Unidentified populations found in 'known_collection': ",
               paste0(unique(stats::na.omit(iden)[which(!stats::na.omit(iden) %in% c(wildpops, hatcheries))]), ", ")))
      }
      iden <- factor(iden, levels = c(wildpops, hatcheries)) %>%
        as.numeric()
    } else {
      iden <- NULL
    }

    # output
    iso_dat = list(
      x = mix,
      y = base,
      iden = iden,
      nalleles = n_alleles,
      groups = base$grpvec,
      group_names = grp_nms,
      wildpops = wildpops,
      hatcheries = hatcheries,
      isoscape = isoscape_ordered
    )

    if (!is.null(file)) saveRDS(iso_dat, file = file)

    print(Sys.time() - start_time)

    return(iso_dat)

  }


#' Allele frequency
#'
#' Calculate allele frequency for each locus
#'   for individual fish or a collection/population.
#'
#' @param gble_in Genotype table.
#' @param gle_ref Reference genetypr table.
#' @param loci loci names.
#' @param collect_by At what level to group by.
#'
#' @noRd
allefreq <- function(gble_in, gble_ref, loci, collect_by = indiv) {

  alleles = lapply(loci, function(loc) {
    dplyr::tibble(locus = loc,
                  call = gble_ref %>%
                    dplyr::select(dplyr::all_of(loc), paste0(loc, ".1")) %>%
                    unlist() %>%
                    unique() %>%
                    .[!is.na(.)],
                  altyp = seq.int(dplyr::n_distinct(call)) %>% factor)
    }) %>% dplyr::bind_rows()

  n_alleles = alleles %>%
    dplyr::group_by(locus) %>%
    dplyr::summarise(n_allele = max(as.numeric(altyp)), .groups = "drop")

  scores_cols = sapply(loci, function(locus) {
    c(locus, paste0(locus, ".1"))
    }) %>%
    as.vector()

  gble_in %>%
    dplyr::select(c({{ collect_by }}, dplyr::all_of(scores_cols))) %>%
    tidyr::pivot_longer(
      cols = -{{ collect_by }},
      names_to = "locus",
      values_to = "allele"
    ) %>%
    dplyr::mutate(
      locus = stringr::str_replace(string = locus, pattern = "\\.1$", replacement = "")
    ) %>%
    dplyr::left_join(alleles,
                     by = c("locus" = "locus", "allele" = "call"),
                     keep = FALSE) %>%
    dplyr::group_by({{ collect_by }}, locus) %>%
    dplyr::count(altyp, .drop = FALSE) %>%
    dplyr::filter(!is.na(altyp)) %>%
    dplyr::left_join(n_alleles,
                     by = c("locus" = "locus"),
                     keep = FALSE) %>%
    dplyr::filter(as.numeric(altyp) <= n_allele) %>%
    dplyr::select(-n_allele) %>%
    tidyr::unite("altyp", c(locus, altyp)) %>%
    tidyr::pivot_wider(names_from = altyp, values_from = n) %>%
    dplyr::ungroup()

}


#' Error check
#'
#' Check loci and population information in input data.
#'
#' @param loci_provided User provided loci 1 info.
#' @param loci_base Loci info from stage 1 baseline.
#' @param loci_mix All loci in mixture data.
#'
#' @noRd
check_loci_pops <- function(loci_provided, loci_base, loci_mix) {

  if (!setequal(loci_base, loci_mix)) {
    return(c("Different loci found in mixture sample and baseline: ",
             paste0(c(setdiff(loci_mix, loci_base), setdiff(loci_base, loci_mix)), ", ")))
  }

  # check loci if provided
  if (!is.null(loci_provided)) {
    if (!setequal(loci_base, loci_provided)) {
      return(
        c("Unidentified loci in baseline or provided list: ",
          paste0(c(setdiff(loci_base, loci_provided), setdiff(loci_provided, loci_base)), ", "))
        )
    }
  }

  return("all good")

}


utils::globalVariables(c(".", "SILLY_CODE", "SillySource", "altyp", "collection",
                         "grpvec", "indiv", "locus", "n", "n_allele", "origin", "repunit"))









