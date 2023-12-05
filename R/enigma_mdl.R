
#' Isotope-genetics BMM
#'
#' @description
#' Run Bayesian mixture model (i.e., Pella-Madusa model) enhanced with isotope information.
#'
#' @param dat_in Input data (as a list object).
#' @param nreps Total number of iterations (includes burn-ins).
#' @param nburn Number of warm-up runs.
#' @param thin Frequency to thin the output.
#' @param nchains Number of independent MCMC processes.
#' @param nadapt Number of adaptation run (default is 0). Only available when
#'   running model in fully Bayesian mode.
#' @param keep_burn Logical (default = FALSE). To save the burn-ins or not.
#' @param cond_gsi Logical (default = TRUE). To run the model in conditional GSI mode.
#' @param file File path to save the output in RDS file. Need to type out the full path and extension `.Rds`. Default = NULL for not saving the output.
#' @param seed Random seed for reproducibility. Default = NULL (no random seed).
#'
#' @return A list containing:
#'  - Summary of the estimates
#'  - Trace of posterior samples
#'  - Individual assignment history
#'
#' @importFrom magrittr %>%
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#'
#' @examples
#' # prep input data
#' enigma_data <- prep_enigma_data(mixture_data = mix_iso, baseline_data = baseline, pop_info = ayk_pops60)
#'
#' # run model
#' enigma_out <- enigma_mdl(enigma_data, 20, 10, 1, 3)
#'
#' @export
enigma_mdl <- function(dat_in, nreps, nburn, thin, nchains, nadapt = 0, keep_burn = FALSE, cond_gsi = TRUE, file = NULL, seed = NULL) {

  ### ballroom categories ### ----
  categories <- c("Live, Werk, Pose", "Bring It Like Royalty", "Face", "Best Mother", "Best Dressed", "High Class In A Fur Coat", "Snow Ball", "Butch Queen Body", "Weather Girl", "Labels", "Mother-Daughter Realness", "Working Girl", "Linen Vs. Silk", "Perfect Tens", "Modele Effet", "Stone Cold Face", "Realness", "Intergalatic Best Dressed", "House Vs. House", "Femme Queen Vogue", "High Fashion In Feathers", "Femme Queen Runway", "Lofting", "Higher Than Heaven", "Once Upon A Time")

  encouragements <- c("I'm proud of ", "keep your ", "be a ", "find strength in ", "worry will never change ", "work for a cause, not for ", "gradtitude turns what we have into ", "good things come to ", "your attitude determines your ", "the only limits in life are ", "find joy in ", "surround yourself with only ", "if oppotunity doesn't knock, build ")

  ### data input ### ----
  x <- dat_in$x %>%
    dplyr::select(ends_with(as.character(0:9))) %>%
    dplyr::select(order(colnames(.))) %>%
    as.matrix() # mixture
  y <- dat_in$y %>%
    dplyr::select(ends_with(as.character(0:9))) %>%
    dplyr::select(order(colnames(.))) %>%
    as.matrix() # base

  sr_val <- dat_in$x$sr_val
  sr_mean <- dat_in$isoscape$sr_mean
  sr_sd <- dat_in$isoscape$sr_sd

  if (is.null(dat_in$iden)) {
    iden <- rep(NA, nrow(x))
  } else iden <- dat_in$iden # iden info

  nalleles <- dat_in$nalleles # number of allele types
  grps <- dat_in$groups # vector id for reporting groups (aka groupvec)
  grp_names <- dat_in$group_names # reporting groups

  wildpops <- dat_in$wildpops
  if (is.null(wildpops)) wildpops <- paste0("yrow_", seq(nrow(y)))
  K <- length(wildpops)

  if (is.null(dat_in$hatcheries)) {
    if (length(grps) > length(wildpops)) {
      stop("There were hatcheries (known pops) in the reporting groups, but no hatchery was found in the baseline data.")
    }
    hatcheries <- NULL
    H <- 0
  } else {
    hatcheries <- dat_in$hatcheries
    H <- length(hatcheries)
  }

  allpops <- c(wildpops, hatcheries)

  if (any(iden > (K + H) | iden > length(grps), na.rm = TRUE)) stop("Unidentified populations in `iden`. Maybe there are hatcheries in the data that are not listed in the reporting groups?")

  na_i <- which(is.na(iden))

  iden <- factor(iden, levels = seq(K + H))

  trait_fac <- factor(rep(names(nalleles), nalleles), levels = names(nalleles))
  states <- paste0(paste0(trait_fac, "_"), unlist(lapply(nalleles, seq.int)))

  # alleles <- states[seq.int(sum(nalleles))] # allele types

  ### specifications ### ----
  rdirich <- function(alpha0) {
    if (sum(alpha0)) {
      vec = stats::rgamma(length(alpha0), alpha0, 1)
      vec = vec / sum(vec)
      vec[vec == 0] = .Machine$double.xmin
      vec
    }
    else{
      rep(0, length(alpha0))
    }
  } # og random dirichlet by jj

  message(paste0("Running model... and ", sample(encouragements, 1), sample(categories, 1), "!"))

  run_time <- Sys.time()

  if (cond_gsi) nadapt = 0
  n_burn <- ifelse(keep_burn, 0, nburn)

  chains <- seq(nchains)
  cl <- parallel::makePSOCKcluster(nchains)
  doParallel::registerDoParallel(cl, cores = nchains)
  if (!is.null(seed)) doRNG::registerDoRNG(seed, once = TRUE)

  ### initial values ### ----
  # hyper-param for relative freq q (allele) and pi (age class)
  beta <- # actually beta and gamma
    matrix(0,
           nrow = nrow(y),
           ncol = ncol(y))#,
           # dimnames = dimnames(y))
  beta[1:K, ] <-
    matrix(
      rep(1 / nalleles, nalleles),
      nrow = K, # number of wildpops (i.e. collection)
      ncol = sum(nalleles),
      byrow = TRUE
    ) # genetic part of prior (beta)

  t_q <- apply(y + beta, 1, function(rw) {
      unlist(tapply(rw, trait_fac, function(betty) {
        if (sum(betty)) {
          betty / sum(betty)
        } else {
          rep(1, length(betty))
        }
      }, simplify = FALSE)[names(nalleles)])
    }) # transposed (allele freq)

  freq <- matrix(
    0,
    nrow = nrow(x),
    ncol = K + H,
    dimnames = list(rownames(x), allpops)
  )

  if (H > 0 & length(na_i) < nrow(x)) {
    freq[-na_i, hatcheries] <-
      t(sapply(as.integer(iden[-na_i]) - K, function(m) {
        ans = rep(0L, H)
        ans[m] = 1L
        ans
      }))
  } # only when both reporting groups and sample include hatcheries

  # genotype freq prod f(x_m|q_k)
  # rows = indiv, cols = pops
  freq[na_i, wildpops] <- exp(x[na_i,] %*% log(t_q[, 1:K]))

  pPrior <- # alpha, hyper-param for p (pop props)
    (1/ table(grps)/ max(grps))[grps]

  iso <-
    sapply(sr_val, function(sr) {
      1 / sqrt((2 * pi * sr_sd^2)) * exp(-1 * (sr - sr_mean)^2 / (2 * sr_sd^2))
    }) %>% t()

  iso <- tidyr::replace_na(iso, replace = 1)

  iden[na_i] <- unlist( lapply(na_i, function(m) {
    sample(K, 1, FALSE, iso[m, ] * (pPrior * freq[m, ])[seq.int(K)])
  }))

  p <- rdirich(table(iden) + pPrior)

  ### parallel chains ### ----
  `%dorng%` <- doRNG::`%dorng%`
  out_list <- foreach::foreach(
    ch = chains, .packages = c("magrittr", "tidyr", "dplyr")
    ) %dorng% {

    p_out <- iden_out <- list()

    ## gibbs loop ##
    for (rep in seq(nreps + nadapt)) {

       if (!cond_gsi & rep > nadapt) { # no cond gsi or passed adapt stage

         x_sum <- matrix(0L, nrow = nrow(y), ncol = ncol(y))

         x_sum[as.integer(sort(unique(iden))),] <-
           rowsum(x, iden) %>%
           tidyr::replace_na(0) # colsums for new assignment

        beta_prm <- y + beta + x_sum # posterior q ~ dirich(b')

        t_q <- apply(beta_prm, 1, function(rw) {
          unlist(tapply(rw, INDEX = trait_fac, FUN = rdirich))
          })

        freq[na_i, wildpops] <- exp(x[na_i,] %*% log(t_q[, 1:K]))

        }

      iden[na_i] <- unlist( lapply(na_i, function(m) {
        sample(K, 1, FALSE, iso[m, ] * (p * freq[m, ])[seq.int(K)])
      }))

      p <- rdirich(table(iden) + pPrior)

      # record output based on keep or not keep burn-ins
      if (rep > nadapt) { # after adaptation stage
        if ((rep-nadapt) > n_burn & (rep - nadapt - n_burn) %% thin == 0) {

          it <- (rep - nadapt - n_burn) / thin
          p_out[[it]] <- c(p, it, ch)

          iden_out[[it]] <- iden

        } # if rep > nburn & (rep-nburn) %% thin == 0
      } # if rep > nadapt

    } # end gibbs loop

    out_items <- list(p_out, iden_out)

    lapply(out_items, function(oi) {
      sapply(oi, rbind) %>%
        t() %>%
        dplyr::as_tibble()
    })

  } # end parallel chains

  parallel::stopCluster(cl)

  ### summary ### ----

  out_list1 <- lapply(out_list, function(ol) ol[[1]])

  p_combo <-
    lapply(out_list1,
           function(ol) ol %>%
             tidyr::pivot_longer(cols = 1:(ncol(.) - 2)) %>%
             dplyr::mutate(
               grpvec = rep(grp_names[grps],
                            (nreps- nburn * isFALSE(keep_burn)) / thin)
             ) %>%
             dplyr::rename(itr = 1, popn = 2) %>%
             dplyr::group_by(itr, grpvec) %>%
             dplyr::summarise(p = sum(value), .groups = "drop") %>%
             tidyr::pivot_wider(names_from = grpvec, values_from = p) %>%
             dplyr::select(-itr))

  keep_list <- ((nburn * keep_burn + 1):(nreps - nburn * isFALSE(keep_burn)))[!((nburn * keep_burn + 1):(nreps - nburn * isFALSE(keep_burn))) %% thin] / thin

  mc_pop <- coda::as.mcmc.list(
    lapply(p_combo,
           function(rlist) coda::mcmc(rlist[keep_list,])) )

  summ_pop <-
    lapply(p_combo, function(rlist) rlist[keep_list,]) %>%
    dplyr::bind_rows() %>%
    tidyr::pivot_longer(cols = 1:ncol(.), names_to = "group") %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      mean = mean(value),
      median = stats::median(value),
      sd = stats::sd(value),
      ci.05 = stats::quantile(value, 0.05),
      ci.95 = stats::quantile(value, 0.95),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      GR = {if (nchains > 1) {
        coda::gelman.diag(mc_pop,
                          transform = FALSE,
                          autoburnin = FALSE,
                          multivariate = FALSE)$psrf[, "Point est."]
      } else {NA}},
      mpsrf = {if (nchains > 1) {
        my.gelman.diag(mc_pop,
                       # transform = FALSE,
                       # autoburnin = FALSE,
                       multivariate = TRUE)$mpsrf
      } else {NA}},
      n_eff = coda::effectiveSize(mc_pop)
    ) %>%
    dplyr::mutate(grp_fac = factor(group, levels = grp_names)) %>%
    dplyr::arrange(grp_fac) %>%
    dplyr::select(-grp_fac)

  print(Sys.time() - run_time)
  message(Sys.time())

  out <- list()

  out$summ <- summ_pop

  out$trace <- p_combo %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      itr = rep(1:((nreps - nburn * isFALSE(keep_burn)) / thin),
                times = nchains),
      chain = rep(1:nchains,
                  each = (nreps - nburn * isFALSE(keep_burn)) / thin)
      )

  out$idens <-
    lapply(out_list, function(ol) ol[[2]]) %>%
    dplyr::bind_rows()

  if(!is.null(file)) saveRDS(out, file = file)

  return(out)

}


#' Gelman-Rubin diagnostics
#'
#' Don't remember where I stole this from
#'
#' @param x MCMC object
#' @param confidence CI level
#' @param multivariate Logical.
#'
#' @noRd
my.gelman.diag <- function(x, confidence = 0.95, multivariate = TRUE) {
  #, transform = FALSE, autoburnin = FALSE) {
  x <- coda::as.mcmc.list(x)
  # if (coda::nchain(x) < 2)
  #   stop("You need at least two chains")
  # if (autoburnin && stats::start(x) < stats::end(x)/2)
  #   x <- stats::window(x, stats::start = stats::end(x)/2 + 1)
  Niter <- coda::niter(x)
  Nchain <- coda::nchain(x)
  Nvar <- coda::nvar(x)
  xnames <- coda::varnames(x)
  # if (transform)
  #   x <- coda::gelman.transform(x)
  x <- lapply(x, as.matrix)
  S2 <- array(sapply(x, stats::var, simplify = TRUE),
              dim = c(Nvar, Nvar, Nchain)
  )
  W <- apply(S2, c(1, 2), mean)
  xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE),
                 nrow = Nvar, ncol = Nchain)
  B <- Niter * stats::var(t(xbar))
  if (Nvar > 1 && multivariate) {  #ph-edits
    # CW <- chol(W)
    #    #This is W^-1*B.
    # emax <- eigen(
    #  backsolve(CW, t(backsolve(CW, B, transpose = TRUE)), transpose = TRUE),
    # symmetric = TRUE, only.values = TRUE)$values[1]
    emax <- 1
    mpsrf <- sqrt((1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter)
  }  else {
    mpsrf <- NULL
  }

  w <- diag(W)
  b <- diag(B)
  s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain)
  muhat <- apply(xbar, 1, mean)
  var.w <- apply(s2, 1, stats::var)/Nchain
  var.b <- (2 * b^2)/(Nchain - 1)
  cov.wb <- (Niter/Nchain) * diag(stats::var(t(s2), t(xbar^2)) - 2 *
                                    muhat * stats::var(t(s2), t(xbar)))
  V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter
  var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b +
              2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
  df.V <- (2 * V^2)/var.V
  df.adj <- (df.V + 3)/(df.V + 1)
  B.df <- Nchain - 1
  W.df <- (2 * w^2)/var.w
  R2.fixed <- (Niter - 1)/Niter
  R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
  R2.estimate <- R2.fixed + R2.random
  R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) *
    R2.random
  psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
  dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
  out <- list(psrf = psrf, mpsrf = mpsrf, B = B, W = W) #added ph
  class(out) <- "gelman.diag"
  return(out)
}


utils::globalVariables(c(".", "ch", "chain", "itr", "name", "name_fac", "value"))













