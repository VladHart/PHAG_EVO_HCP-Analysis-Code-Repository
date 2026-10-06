# ==============================================================================
# PHAG-EVO_CH HCP questionnaire: complete RAW-RESPONSE-ONLY analysis
# Version: 5.0.0 | 2026-10-05
#
# Run in RStudio:
#   source("PHAG_EVO_HCP_complete_raw_only_analysis_v5.R")
# Or at the command line:
#   Rscript PHAG_EVO_HCP_complete_raw_only_analysis_v5.R input.csv results_directory
#
# Statistical status:
# - Original primary specification: nine-item Pearson PAF, regression factor score,
#   five-category proportional-odds outcome, adjustment for dressing and diagnosis.
# - Added diagnostic follow-up is exploratory, not retrospectively prespecified.
# - A failed diagnostic is recorded as failed/inconclusive, NEVER as a pass.
# - A one-factor candidate is not automatically validated by positive loadings.
# - Polychoric repair is explicit, exported, and not evidence of valid measurement.
# - Q-mode is patient-profile PCA, NOT identification of unique HCP viewpoints.
# - Site and user-account clusters are NOT verified clinical-rater clusters.
#
# ==============================================================================

HCP_SCRIPT_VERSION <- "5.0.0_raw_responses_only"
HCP_SCRIPT_PATH <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
                           error = function(e) NA_character_)
if (length(HCP_SCRIPT_PATH) != 1L) HCP_SCRIPT_PATH <- NA_character_

hcp_default_config <- function() {
  list(seed = 20260407L, n_pa = 5000L, n_pa_common = 1000L,
       n_pa_poly = 1000L, n_pa_q = 2000L, n_boot = 1000L,
       n_boot_site = 1000L, n_boot_cor = 2000L,
       quantile = 0.95, min_sim_success = 0.90,
       # Same principal-axis convergence criterion as the previous R analysis.
       # This preserves its primary point estimates; MINRES is an added sensitivity.
       pa_min_err = 0.001, pa_max_iter = 1000L,
       poly_correct = 0.5, poly_global = FALSE,
       poly_eigen_tolerance = 1e-8,
       diagnostic_gradient_tolerance = 1e-4,
       diagnostic_condition_limit = 1e10,
       install_missing = FALSE, make_figures = TRUE,
       save_models = TRUE)
}

hcp_analysis <- function(input_file, output_root = "PHAG_EVO_HCP_raw_only_outputs",
                         config = list()) {
  cfg <- utils::modifyList(hcp_default_config(), config)
  counts <- c("n_pa","n_pa_common","n_pa_poly","n_pa_q","n_boot","n_boot_site","n_boot_cor")
  if(any(!vapply(cfg[counts],function(x) length(x)==1L && is.finite(x) && x>=2 && x==as.integer(x),logical(1))))
    stop("Simulation and bootstrap counts must be integers of at least 2.")
  old_options <- options(stringsAsFactors = FALSE, digits = 17)
  on.exit(options(old_options), add = TRUE)
  deps <- c("psych", "GPArotation", "MASS", "ordinal", "sandwich", "survey",
            "brglm2", "digest")
  missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && isTRUE(cfg$install_missing)) {
    utils::install.packages(missing, repos = "https://cloud.r-project.org",
                            dependencies = c("Depends", "Imports", "LinkingTo"))
    missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
  }
  if (length(missing)) stop("Install the required packages, then rerun:\n",
    "install.packages(c(", paste(sprintf('"%s"', missing), collapse = ", "),
    "), repos = 'https://cloud.r-project.org')", call. = FALSE)
  if (!file.exists(input_file)) stop("Input CSV not found: ", input_file,
    "\nPlace the raw analytical CSV beside the script or set options(hcp.input = 'full/path.csv').")
  input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
  run_label <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  out <- file.path(output_root, run_label)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) stop("Cannot create output directory: ", out)
  dir.create(file.path(out, "tables")); dir.create(file.path(out, "figures"))
  dir.create(file.path(out, "private")); dir.create(file.path(out, "models"))
  audit <- new.env(parent = emptyenv())
  audit$tables <- list(); audit$descriptions <- list(); audit$privacy <- list()
  audit$events <- list(); audit$run_status <- "incomplete"; audit$fatal <- ""
  audit$start <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  std <- hcp_default_config()
  audit$run_mode <- if(any(vapply(counts,function(nm) cfg[[nm]]<std[[nm]],logical(1))))
    "reduced_iteration_run_NOT_for_manuscript" else "standard_or_higher_iteration_run"
  audit$step <- 0L
  message("Raw-only analysis output: ", normalizePath(out, winslash = "/"))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(cfg$seed)

  # ---- Audit/export utilities -------------------------------------------------
  bind_rows_base <- function(xs) {
    xs <- Filter(function(x) !is.null(x) && is.data.frame(x), xs)
    if (!length(xs)) return(data.frame())
    nm <- unique(unlist(lapply(xs, names), use.names = FALSE))
    xs <- lapply(xs, function(x) {
      for (n in setdiff(nm, names(x))) x[[n]] <- rep(NA, nrow(x))
      x[, nm, drop = FALSE]
    })
    ans <- do.call(rbind, xs); rownames(ans) <- NULL; ans
  }
  event <- function(section, type, detail) {
    audit$events[[length(audit$events) + 1L]] <- data.frame(
      section = section, type = type, detail = paste(detail, collapse = " | "))
  }
  capture <- function(label, expr) {
    ww <- character(); mm <- character(); err <- NULL
    val <- tryCatch(withCallingHandlers(force(expr),
      warning = function(w) {ww <<- c(ww, conditionMessage(w)); invokeRestart("muffleWarning")},
      message = function(m) {mm <<- c(mm, conditionMessage(m)); invokeRestart("muffleMessage")}),
      error = function(e) {err <<- conditionMessage(e); NULL})
    if (length(ww)) event(label, "warning", unique(ww))
    if (length(mm)) event(label, "message", unique(mm))
    if (!is.null(err)) event(label, "error", err)
    list(value = val, error = err, warnings = unique(ww), messages = unique(mm))
  }
  emit <- function(name, x, description = "", private = FALSE) {
    if (is.matrix(x)) x <- data.frame(row_id = rownames(x), x, check.names = FALSE)
    if (!is.data.frame(x)) x <- as.data.frame(x, check.names = FALSE)
    rownames(x) <- NULL
    if (any(vapply(x, is.list, logical(1)))) stop("Non-tabular list column in ", name)
    audit$tables[[name]] <- x; audit$descriptions[[name]] <- description
    audit$privacy[[name]] <- private
    utils::write.csv(x, file.path(out, if (private) "private" else "tables",
                                  paste0(name, ".csv")), row.names = FALSE, na = "NA",
                     fileEncoding = "UTF-8")
    invisible(x)
  }
  emit_failure <- function(name, result, description = "") {
    emit(name, data.frame(status = "not_estimable", reason =
      if (!is.null(result$error)) result$error else "No result returned"), description)
  }
  matrix_df <- function(M, id = "item") {
    ans <- data.frame(row = rownames(M), M, check.names = FALSE)
    names(ans)[1] <- id; ans
  }
  scalar <- function(x, default = NA_real_) {
    if (is.null(x) || !length(x)) default else as.numeric(x[1])
  }
  num_text <- function(x) {
    if (length(x) != 1L || is.na(x)) return(NA_character_)
    if (is.numeric(x)) return(format(x, digits = 17, scientific = NA, trim = TRUE))
    as.character(x)
  }
  finalize <- function() {
    if(audit$run_status=="incomplete" && !nzchar(audit$fatal))
      audit$fatal <- "Execution interrupted; inspect the R console and condition log. Partial outputs are not a completed analysis."
    ev <- bind_rows_base(audit$events)
    if (!nrow(ev)) ev <- data.frame(section = "run", type = "information", detail = "No captured conditions")
    emit("98_warnings_errors_and_messages", ev, "Conditions are retained, not silently discarded")
    pkg <- data.frame(package = deps, version = vapply(deps, function(p)
      as.character(utils::packageVersion(p)), character(1)))
    emit("99_package_versions", pkg)
    utils::capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
    meta <- data.frame(record_type = "metadata", analysis_table = "run_metadata",
      source_row = NA_integer_, measure = c("script_version", "run_status", "run_mode", "input_sha256",
      "input_file", "started_utc", "ended_utc", "R_version", "fatal_error", "confidentiality"),
      value = c(HCP_SCRIPT_VERSION, audit$run_status, audit$run_mode,
        digest::digest(file = input_file, algo = "sha256"), basename(input_file),
        audit$start, format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        R.version.string, audit$fatal,
        "Contains pseudonymised patient-profile data. Controlled research use; not public-ready."),
      is_missing = FALSE, stringsAsFactors = FALSE)
    pieces <- list(meta)
    for (nm in names(audit$tables)) {
      x <- audit$tables[[nm]]
      st <- data.frame(record_type = "table_status", analysis_table = nm,
        source_row = NA_integer_, measure = c("status", "n_rows", "n_columns", "description", "private"),
        value = c(if (nrow(x)) "exported" else "empty", nrow(x), ncol(x),
                  audit$descriptions[[nm]], audit$privacy[[nm]]), is_missing = FALSE)
      pieces[[length(pieces) + 1L]] <- st
      keys <- intersect(c("analysis_set","model","code","item","term","outcome",
                          "diagnosis","component","participant_id","scheme","replicate"),names(x))
      labels <- if(length(keys) && nrow(x)) vapply(seq_len(nrow(x)),function(i)
        paste(vapply(keys,function(k) paste0(k,"=",num_text(x[[k]][i])),character(1)),
              collapse="; "),character(1)) else paste0("row=",seq_len(nrow(x)))
      if (nrow(x) && ncol(x)) for (j in seq_along(x)) {
        vv <- vapply(seq_len(nrow(x)), function(i) num_text(x[[j]][i]), character(1))
        pieces[[length(pieces) + 1L]] <- data.frame(record_type = "result",
          analysis_table = nm, source_file = paste0(nm,".csv"),
          source_row = seq_len(nrow(x)), source_column=j,row_label=labels,measure = names(x)[j],
          value = vv, is_missing = is.na(vv), stringsAsFactors = FALSE)
      }
    }
    all_results <- bind_rows_base(pieces)
    all_results$value_numeric <- suppressWarnings(as.numeric(all_results$value))
    utils::write.csv(all_results, file.path(out, "PHAG_EVO_HCP_ALL_RESULTS_RAW_ONLY.csv"),
                     row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
    saveRDS(list(results = audit$tables, configuration = cfg, version = HCP_SCRIPT_VERSION,
                 run_status = audit$run_status), file.path(out, "all_tabulated_results.rds"))
    ff <- list.files(out, recursive = TRUE, full.names = TRUE)
    ff <- ff[!dir.exists(ff) & !grepl("SHA256SUMS", ff, fixed = TRUE)]
    hashes <- vapply(ff, function(f) digest::digest(file = f, algo = "sha256"), character(1))
    writeLines(paste(hashes, substring(ff, nchar(out) + 2L)), file.path(out, "SHA256SUMS.txt"))
    message("Run status: ", audit$run_status, "\nCombined results: ",
      normalizePath(file.path(out, "PHAG_EVO_HCP_ALL_RESULTS_RAW_ONLY.csv"), winslash = "/"))
  }
  on.exit(tryCatch(finalize(), error = function(e)
    message("Final export error: ", conditionMessage(e))), add = TRUE)
  step <- function(label) {audit$step <- audit$step + 1L; message("[", audit$step, "] ", label)}
  save_model <- function(fit, name) {
    if (isTRUE(cfg$save_models)) saveRDS(fit, file.path(out, "models", paste0(name, ".rds")))
  }

  # ---- Input validation: explicit raw-variable allow-list ----------------------
  step("Read and validate original responses")
  source <- utils::read.csv(input_file, colClasses = "character", check.names = FALSE,
                            na.strings = c("", "NA", "NaN"), fileEncoding = "UTF-8-BOM")
  stems <- c("ease_application", "body_contours", "wound_conformability", "exudate_management",
             "moist_wound_environment", "periwound_protection", "ease_removal",
             "nonadherence_removal", "stay_in_place", "overall_impression")
  raw_cols <- paste0("hcp_", stems, "_raw"); domain_cols <- raw_cols[1:9]
  item_names <- c("Ease of application", "Ability to conform to body contours",
    "Conformability to the wound", "Effective exudate management",
    "Ability to maintain moist wound environment", "Ability to protect peri-wound skin",
    "Ease of removal", "Non-adherence to wound at removal", "Ability to stay in place",
    "Overall impression of dressing")
  item_codes <- c(paste0("A", 1:9), "G1")
  source_cols <- paste0("hcp_", stems, "_source_response")
  context_cols <- c("participant_id", "analysis_population", "site_prefix", "dressing", "diagnosis",
                    "ecrf_user_id", "ecrf_account_label", "ecrf_site_prefix")
  outcomes <- c("war_pre", "war_post", "success_pre_20pct", "success_post_20pct", "pain_final",
                "final_exudate_volume", "final_daily_exudate_retention",
                "final_dressing_changes", "final_wear_time")
  area_cols <- c("wound_area_pre_baseline", "wound_area_pre_final_available",
                 "wound_area_post_baseline", "wound_area_post_final_available")
  needed <- c("participant_id", "dressing", "diagnosis", raw_cols)
  absent <- setdiff(needed, names(source))
  if (length(absent)) stop("Required raw-data columns absent: ", paste(absent, collapse = ", "))
  d <- source[, intersect(c(context_cols, raw_cols, source_cols, area_cols, outcomes), names(source)), drop = FALSE]
  for (v in setdiff(c("site_prefix", "ecrf_user_id", "ecrf_account_label"), names(d))) d[[v]] <- NA_character_
  if (anyNA(d$participant_id) || anyDuplicated(d$participant_id))
    stop("Participant IDs must be nonmissing and unique. No automatic deduplication is allowed.")
  numeric_cols <- intersect(c(raw_cols, area_cols, outcomes), names(d))
  for (v in numeric_cols) {
    tmp <- suppressWarnings(as.numeric(d[[v]]))
    if (any(!is.na(d[[v]]) & is.na(tmp))) stop("Non-numeric value in ", v)
    if (any(is.infinite(tmp))) stop("Infinite value in ", v, "; resolve before analysis")
    d[[v]] <- tmp
  }
  for (v in raw_cols) if (any(!is.na(d[[v]]) & !(d[[v]] %in% 1:5)))
    stop("Invalid questionnaire response in ", v, "; allowed: 1,2,3,4,5 or missing.")
  if (any(!d$dressing %in% c("RSSB", "RSSil")) || anyNA(d$dressing)) stop("Unexpected/missing dressing label")
  if (any(!d$diagnosis %in% c("VLU", "DFU")) || anyNA(d$diagnosis)) stop("Unexpected/missing diagnosis label")
  d$dressing <- factor(d$dressing, levels = c("RSSB", "RSSil"))
  d$diagnosis <- factor(d$diagnosis, levels = c("VLU", "DFU"))
  if ("analysis_population" %in% names(d) && any(d$analysis_population != "ITT", na.rm = TRUE))
    stop("Input contains non-ITT records. Supply the approved analytical population.")
  if (nrow(d) != 80L) event("data_validation", "warning", "Input size differs from the 80-patient supplied analytical dataset")
  emit("00_input_configuration", data.frame(setting = names(cfg), value = vapply(cfg, as.character, character(1))))
  emit("01_analysis_population", data.frame(characteristic = c("ITT patients", "RSSB", "RSSil", "VLU", "DFU"),
    n = c(nrow(d), sum(d$dressing == "RSSB"), sum(d$dressing == "RSSil"), sum(d$diagnosis == "VLU"), sum(d$diagnosis == "DFU"))))
  emit("02_item_dictionary", data.frame(code = item_codes, variable = raw_cols, item = item_names,
    role = c(rep("domain item", 9), "global judgement outcome; excluded from score derivation")))
  # Source-response check is a consistency check, not source-document adjudication.
  source_check <- list()
  for (j in seq_along(raw_cols)) if (source_cols[j] %in% names(d)) {
    text <- trimws(d[[source_cols[j]]]); has_code <- !is.na(text) & grepl("^[1-5]([[:space:]]*[-:]|$)", text)
    parsed <- rep(NA_real_, nrow(d)); parsed[has_code] <- as.numeric(substr(text[has_code], 1, 1))
    mismatch <- which(has_code & (is.na(d[[raw_cols[j]]]) | parsed != d[[raw_cols[j]]]))
    source_check[[j]] <- data.frame(code = item_codes[j], available_source_text = sum(!is.na(text)),
      parsed_1_to_5 = sum(has_code), mismatches = length(mismatch),
      category_1_retained = sum(d[[raw_cols[j]]] == 1, na.rm = TRUE))
    if (length(mismatch)) stop("Raw/source-response mismatch in ", raw_cols[j], "; verify input before analysis")
  }
  emit("03_source_response_consistency", bind_rows_base(source_check))
  # Validate supplied wound-area-derived variables, without silently replacing them.
  area_check <- list()
  for (s in c("pre", "post")) {
    vv <- c(paste0("wound_area_", s, "_baseline"), paste0("wound_area_", s, "_final_available"), paste0("war_", s))
    if (all(vv %in% names(d))) {
      ok <- complete.cases(d[, vv]) & d[[vv[1]]] > 0
      predicted <- 100 * (d[[vv[1]]][ok] - d[[vv[2]]][ok]) / d[[vv[1]]][ok]
      delta <- abs(predicted - d[[vv[3]]][ok])
      suc <- paste0("success_", s, "_20pct")
      km <- if (suc %in% names(d)) which(!is.na(d[[suc]]) & !is.na(d[[vv[3]]]) &
                                           d[[suc]] != as.numeric(d[[vv[3]]] > 20)) else integer()
      area_check[[s]] <- data.frame(measure = s, checked_n = sum(ok),
        max_absolute_WAR_difference = if (length(delta)) max(delta) else NA_real_,
        threshold_indicator_mismatches = length(km))
      if ((length(delta) && max(delta) > 1e-5) || length(km))
        stop("Supplied wound-area or >20% indicator inconsistency (", s, "); resolve before use")
    }
  }
  emit("04_wound_area_derivation_checks", bind_rows_base(area_check))
  emit("05_raw_analysis_dataset", d, "Approved input columns only; unmodified raw responses", private = TRUE)
  rm(source)

  # ---- Descriptives, completion and context -----------------------------------
  step("Describe responses, missingness and case-level assessment context")
  prop_ci <- function(k, n) {if (!n) return(c(NA_real_, NA_real_));
    z <- stats::qnorm(.975); p <- k/n; mid <- (p+z*z/(2*n))/(1+z*z/n)
    half <- z*sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n); c(mid-half, mid+half)}
  desc <- lapply(seq_along(raw_cols), function(j) {
    x <- d[[raw_cols[j]]]; y <- x[!is.na(x)]; n <- length(y)
    counts <- tabulate(as.integer(y), nbins = 5); ci <- prop_ci(sum(counts[4:5]), n)
    data.frame(code = item_codes[j], item = item_names[j], variable = raw_cols[j],
      valid_n = n, missing_n = sum(is.na(x)), mean = if(n) mean(y) else NA_real_,
      sd = if(n>1) stats::sd(y) else NA_real_, median = if(n) stats::median(y) else NA_real_,
      q1 = if(n) stats::quantile(y,.25,names=FALSE) else NA_real_,
      q3 = if(n) stats::quantile(y,.75,names=FALSE) else NA_real_,
      count_1 = counts[1], count_2 = counts[2], count_3 = counts[3], count_4 = counts[4], count_5 = counts[5],
      top_two_n = sum(counts[4:5]), top_two_pct = if(n) 100*sum(counts[4:5])/n else NA_real_,
      top_two_Wilson_low_pct = 100*ci[1], top_two_Wilson_high_pct = 100*ci[2],
      floor_1_pct = if(n) 100*counts[1]/n else NA_real_, ceiling_5_pct = if(n) 100*counts[5]/n else NA_real_)
  })
  desc <- emit("06_item_descriptives", bind_rows_base(desc),
    "Top two categories = scores 4 and 5. Wilson intervals are descriptive, unclustered intervals.")
  cc9 <- complete.cases(d[, domain_cols]); cc10 <- complete.cases(d[, raw_cols])
  X <- as.matrix(d[cc9, domain_cols]); colnames(X) <- item_codes[1:9]
  X10 <- as.matrix(d[cc10, raw_cols]); colnames(X10) <- item_codes
  if (nrow(X) <= ncol(X)+2L || any(apply(X,2,stats::sd)==0))
    stop("Insufficient complete cases or a constant domain item; factor analysis cannot proceed")
  emit("07_completeness", data.frame(analysis_set = c("9 domains", "10 items"),
    complete_n = c(sum(cc9),sum(cc10)), total_n=nrow(d),
    complete_pct=100*c(mean(cc9),mean(cc10))))
  emit("08_case_inclusion", data.frame(participant_id=d$participant_id, nine_complete=cc9,
    ten_complete=cc10, overall_available=!is.na(d[[raw_cols[10]]])), private=TRUE)
  context <- list()
  for (v in c("site_prefix","ecrf_user_id")) {
    tt <- table(d[[v]], useNA="ifany")
    context[[v]] <- data.frame(context_type=rep(v,length(tt)), context_level=names(tt), n=as.integer(tt))
  }
  emit("09_site_and_account_context", bind_rows_base(context),
    "User accounts indicate entry/modification, not verified HCP raters; denominators use this analytical dataset")

  # ---- Factorability and EFA ---------------------------------------------------
  step("Pearson structural analyses and internal consistency")
  R <- stats::cor(X); R10 <- stats::cor(X10)
  emit("10_pearson_correlation_9_domains", matrix_df(R))
  emit("11_pearson_correlation_10_items", matrix_df(R10))
  alpha_raw <- function(A) {k<-ncol(A); total<-stats::var(rowSums(A));
    if(k<2 || total<=0) return(NA_real_); k/(k-1)*(1-sum(apply(A,2,stats::var))/total)}
  standardized_alpha <- function(M) {p<-ncol(M); p/(p-1)*(1-p/sum(M))}
  factorability <- function(A, label) {
    M <- stats::cor(A); ee <- eigen(M,symmetric=TRUE,only.values=TRUE)$values
    km <- capture(paste0(label," KMO"), psych::KMO(M))
    ch <- -(nrow(A)-1-(2*ncol(A)+5)/6)*as.numeric(determinant(M,logarithm=TRUE)$modulus)
    df <- ncol(A)*(ncol(A)-1)/2
    data.frame(analysis_set=label, n=nrow(A), items=ncol(A), KMO=if(is.null(km$value)) NA_real_ else km$value$MSA,
      Bartlett_chi2=ch,Bartlett_df=df,Bartlett_p=stats::pchisq(ch,df,lower.tail=FALSE),
      Bartlett_log_p=stats::pchisq(ch,df,lower.tail=FALSE,log.p=TRUE),
      alpha_raw=alpha_raw(A),alpha_standardized=standardized_alpha(M),
      min_eigen=min(ee),condition_number=max(ee)/min(ee))
  }
  factab <- emit("12_factorability_and_alpha",rbind(factorability(X,"9-domain primary"),
    factorability(X10,"10-item descriptive sensitivity; includes outcome")))
  kmo9 <- capture("item KMO",psych::KMO(R))
  if(!is.null(kmo9$value)) emit("13_item_KMO",data.frame(code=item_codes[1:9],item=item_names[1:9],KMO=as.numeric(kmo9$value$MSAi)))
  fa_fit <- function(A, nf=1L, method="pa", matrix_input=FALSE, n_obs=nrow(A)) {
    args<-list(r=A,nfactors=nf,fm=method,rotate=if(nf>1L) "oblimin" else "none",
               scores=if(matrix_input) "none" else "regression",max.iter=cfg$pa_max_iter)
    if(method=="pa") args$min.err<-cfg$pa_min_err
    if(matrix_input) args$n.obs<-n_obs
    do.call(psych::fa,args)
  }
  load_table <- function(f, model_name) {
    L<-unclass(f$loadings); if(ncol(L)==1L && sum(L)<0) L<- -L
    data.frame(model=model_name,code=rownames(L),L,
      communality=as.numeric(f$communality),uniqueness=as.numeric(f$uniquenesses),check.names=FALSE)
  }
  efa1c <- capture("Pearson PAF one factor",fa_fit(X))
  if(is.null(efa1c$value)) stop("Primary PAF failed: ",efa1c$error)
  efa1<-efa1c$value; sign1<-if(sum(unclass(efa1$loadings))<0) -1 else 1
  L1<-as.numeric(unclass(efa1$loadings)[,1])*sign1
  score_unscaled<-as.numeric(efa1$scores[,1])*sign1
  if(anyNA(score_unscaled)||stats::sd(score_unscaled)==0) stop("Invalid primary regression factor scores")
  global_score<-as.numeric(scale(score_unscaled))
  emit("14_PAF_one_factor_loadings",load_table(efa1,"one-factor candidate"))
  save_model(efa1,"efa_primary_PAF_one_factor")
  candidates<-list(one_factor=efa1)
  for (nf in c(2L,3L)) {
    rr<-capture(paste0("PAF ",nf," factor exploratory"),fa_fit(X,nf))
    nm<-paste0("15_PAF_",nf,"_factor_exploratory_loadings")
    if(!is.null(rr$value)) {candidates[[paste0(nf,"_factors")]]<-rr$value
      emit(nm,load_table(rr$value,paste0(nf,"-factor exploratory")))
      if(!is.null(rr$value$Phi)) emit(paste0(nm,"_factor_correlations"),matrix_df(rr$value$Phi,"factor"))
      save_model(rr$value,paste0("efa_PAF_",nf,"_factors"))
    } else emit_failure(nm,rr)
  }
  minres<-capture("Pearson MINRES sensitivity",fa_fit(X,method="minres"))
  if(!is.null(minres$value)) emit("16_MINRES_one_factor_sensitivity",load_table(minres$value,"MINRES sensitivity")) else emit_failure("16_MINRES_one_factor_sensitivity",minres)
  rhoS<-stats::cor(X,method="spearman")
  emit("17_spearman_correlation_sensitivity",matrix_df(rhoS))
  sfa<-capture("Spearman PAF sensitivity",fa_fit(rhoS,matrix_input=TRUE,n_obs=nrow(X)))
  if(!is.null(sfa$value)) emit("18_spearman_PAF_sensitivity",load_table(sfa$value,"Spearman sensitivity")) else emit_failure("18_spearman_PAF_sensitivity",sfa)
  fit_metrics<-lapply(names(candidates),function(nm) {
    f<-candidates[[nm]]; LL<-as.matrix(unclass(f$loadings))
    Phi<-if(is.null(f$Phi)) diag(ncol(LL)) else f$Phi
    residual<-R-LL%*%Phi%*%t(LL); diag(residual)<-0
    data.frame(model=nm,n=nrow(X),df=scalar(f$dof),RMSR=scalar(f$rms),
      max_abs_offdiagonal_residual=max(abs(residual[lower.tri(residual)])),
      TLI=scalar(f$TLI),RMSEA=scalar(f$RMSEA),RMSEA_lower=if(length(f$RMSEA)>=2) f$RMSEA[2] else NA_real_,
      RMSEA_upper=if(length(f$RMSEA)>=3) f$RMSEA[3] else NA_real_,
      Heywood=any(f$communality>1|f$uniquenesses<0))
  })
  emit("19_PAF_model_fit_descriptive",bind_rows_base(fit_metrics),
    "Fit statistics are descriptive for ordinal, clustered data; positive loadings alone do not establish unidimensionality")
  residual<-R-tcrossprod(L1); diag(residual)<-0
  emit("20_one_factor_offdiagonal_residuals",matrix_df(residual))
  omega1<-(sum(L1)^2)/(sum(L1)^2+sum(efa1$uniquenesses))
  emit("21_one_factor_model_based_omega",data.frame(omega_congeneric=omega1,
    note="Model-implied congeneric omega for the candidate one-factor model; depends on its adequacy, not test-retest/inter-rater reliability"))
  item_diag<-lapply(seq_len(ncol(X)),function(j) data.frame(code=colnames(X)[j],
    corrected_item_total_r=stats::cor(X[,j],rowSums(X[,-j,drop=FALSE])),
    alpha_if_deleted=alpha_raw(X[,-j,drop=FALSE])))
  emit("22_item_total_and_alpha_if_deleted",bind_rows_base(item_diag))
  W<-as.matrix(efa1$weights)*sign1
  emit("23_factor_scoring_coefficients",data.frame(code=colnames(X),input_mean=colMeans(X),
    input_sd=apply(X,2,stats::sd),regression_weight=W[,1],
    score_mean=mean(score_unscaled),score_sd=stats::sd(score_unscaled)),
    "Apply weights to item z scores, then standardize to the recorded score mean and SD")

  # ---- Parallel analysis: PCA and common-factor diagnostics are distinct -------
  step("Parallel analysis with separate component and common-factor specifications")
  eig_pc<-function(M) eigen(M,symmetric=TRUE,only.values=TRUE)$values
  # Reduced-matrix SMC PA is an explicit common-factor sensitivity, NOT PCA.
  # SMC estimates may favour minor factors; retain the full diagnostics rather than
  # selecting a favourable retention rule after observing results.
  eig_common<-function(M) {h<-1-1/diag(solve(M)); diag(M)<-h; eig_pc(M)}
  pa_summary<-function(obs, sims, label, null_label) {
    valid<-complete.cases(sims); nr<-sum(valid); success<-nr/nrow(sims)
    qv<-if(nr) apply(sims[valid,,drop=FALSE],2,stats::quantile,probs=cfg$quantile,names=FALSE) else rep(NA_real_,length(obs))
    exceed<-obs>qv; sequential<-if(success<cfg$min_sim_success) rep(NA,length(obs)) else as.logical(cumprod(as.integer(exceed)))
    data.frame(analysis=label,null_model=null_label,component=seq_along(obs),observed_eigen=obs,
      random_mean=if(nr) colMeans(sims[valid,,drop=FALSE]) else NA_real_,random_p95=qv,
      exceeds_p95=exceed,retained_sequential=sequential,successful_replicates=nr,
      requested_replicates=nrow(sims),success_fraction=success)
  }
  pa_run<-function(A,B,seed,null="gaussian",common=FALSE) {
    set.seed(seed); p<-ncol(A); n<-nrow(A); fun<-if(common) eig_common else eig_pc
    obs<-fun(stats::cor(A)); sims<-matrix(NA_real_,B,p)
    for(b in seq_len(B)) {
      Z<-if(null=="gaussian") matrix(stats::rnorm(n*p),n,p) else
        vapply(seq_len(p),function(j) sample(A[,j],n,replace=FALSE),numeric(n))
      rr<-tryCatch(fun(stats::cor(Z)),error=function(e) rep(NA_real_,p)); sims[b,]<-rr
    }
    pa_summary(obs,sims,if(common) "common-factor reduced matrix with SMC diagonal" else "PCA full correlation matrix",null)
  }
  pa9<-emit("24_PCA_parallel_gaussian_9",pa_run(X,cfg$n_pa,cfg$seed))
  pa9perm<-emit("25_PCA_parallel_permutation_9",pa_run(X,cfg$n_pa,cfg$seed+1L,"permutation"))
  pa10<-emit("26_PCA_parallel_gaussian_10",pa_run(X10,cfg$n_pa,cfg$seed))
  pacom<-capture("common-factor SMC PA",pa_run(X,cfg$n_pa_common,cfg$seed+2L,"permutation",TRUE))
  if(!is.null(pacom$value)) emit("27_common_factor_parallel_SMC",pacom$value,
    "Separate common-factor retention sensitivity; potential minor-factor sensitivity of SMC is acknowledged") else emit_failure("27_common_factor_parallel_SMC",pacom)
  pc<-eigen(R,symmetric=TRUE); pcaload<-sweep(pc$vectors,2,sqrt(pmax(pc$values,0)),"*")
  rownames(pcaload)<-colnames(X); colnames(pcaload)<-paste0("PC",seq_len(ncol(X)))
  emit("28_PCA_loadings",matrix_df(pcaload))
  emit("29_PCA_variance",data.frame(component=seq_len(ncol(X)),eigenvalue=pc$values,
    variance_pct=100*pc$values/ncol(X),cumulative_pct=100*cumsum(pc$values)/ncol(X),
    Kaiser_gt_1=pc$values>1))

  # ---- Polychoric estimation, non-PD audit and matched ordinal null ------------
  step("Polychoric sensitivity with unsmoothed matrices and explicit repair audit")
  poly_estimate<-function(A,label) {
    p<-ncol(A)
    pp<-capture(label,psych::polychoric(as.data.frame(A),smooth=FALSE,
      correct=cfg$poly_correct,global=cfg$poly_global,progress=FALSE,max.cat=5))
    if(is.null(pp$value)) return(list(ok=FALSE,error=pp$error))
    M<-pp$value$rho
    if(!is.matrix(M)||any(!is.finite(M))||nrow(M)!=p)
      return(list(ok=FALSE,error="Nonfinite/incomplete polychoric matrix"))
    original<-M; ee<-eig_pc(M); repaired<-min(ee)<=cfg$poly_eigen_tolerance
    used<-M
    if(repaired) {
      sr<-capture(paste0(label," PD repair"),psych::cor.smooth(M,eig.tol=cfg$poly_eigen_tolerance))
      if(is.null(sr$value)) return(list(ok=FALSE,error=sr$error,original=original))
      used<-sr$value
    }
    eu<-eig_pc(used)
    if(min(eu)<=0 || any(!is.finite(used)))
      return(list(ok=FALSE,error="Matrix still non-positive-definite after explicit repair",original=original))
    list(ok=TRUE,original=original,used=used,thresholds=pp$value$tau,repaired=repaired,
      min_original=min(ee),min_used=min(eu),condition=max(eu)/min(eu),
      max_change=max(abs(used-original)),frobenius_change=sqrt(sum((used-original)^2)),
      warnings=paste(pp$warnings,collapse=" | "))
  }
  # Audit all pairwise ordinal tables, including sparse/zero cells.
  pairs<-utils::combn(seq_len(ncol(X)),2); cell_rows<-list(); pair_rows<-list()
  for(i in seq_len(ncol(pairs))) {
    a<-pairs[1,i]; b<-pairs[2,i]
    ct<-table(factor(X[,a],levels=1:5),factor(X[,b],levels=1:5))
    tmp<-as.data.frame(ct); names(tmp)<-c("response_a","response_b","n")
    tmp$item_a<-colnames(X)[a]; tmp$item_b<-colnames(X)[b]; cell_rows[[i]]<-tmp
    ct_obs<-ct[rowSums(ct)>0,colSums(ct)>0,drop=FALSE]
    pair_rows[[i]]<-data.frame(item_a=colnames(X)[a],item_b=colnames(X)[b],
      observed_levels_a=nrow(ct_obs),observed_levels_b=ncol(ct_obs),
      zero_cells_observed_margins=sum(ct_obs==0),cells_under_5=sum(ct_obs<5),
      total_cells=nrow(ct_obs)*ncol(ct_obs))
  }
  emit("30_polychoric_pairwise_counts",bind_rows_base(cell_rows))
  emit("31_polychoric_sparsity",bind_rows_base(pair_rows))
  poly<-poly_estimate(X,"observed polychoric")
  polyfa<-NULL; poly_pa<-NULL
  if(isTRUE(poly$ok)) {
    emit("32_polychoric_matrix_UNSMOOTHED",matrix_df(poly$original),
      "Original pairwise estimates. Do not use as a covariance matrix if non-positive-definite.")
    emit("33_polychoric_matrix_USED",matrix_df(poly$used),
      "Explicitly repaired only when required; see matrix audit. Not independent confirmation of validity.")
    emit("34_polychoric_matrix_adjustment",matrix_df(poly$used-poly$original))
    emit("35_polychoric_matrix_audit",data.frame(n=nrow(X),items=ncol(X),
      estimator="psych::polychoric two-step",global_thresholds=cfg$poly_global,
      zero_cell_correction=cfg$poly_correct,repair_required=poly$repaired,
      repair_method=if(poly$repaired) "psych::cor.smooth" else "none",
      repair_eigen_tolerance=cfg$poly_eigen_tolerance,min_eigen_unsmoothed=poly$min_original,
      min_eigen_used=poly$min_used,condition_used=poly$condition,
      max_absolute_correlation_change=poly$max_change,frobenius_change=poly$frobenius_change,
      determinant_used=det(poly$used),warnings=poly$warnings))
    if(!is.null(poly$thresholds)) {
      tt<-capture("polychoric thresholds export",as.data.frame(poly$thresholds))
      if(!is.null(tt$value)) {tt$value$source_row<-rownames(tt$value);emit("36_polychoric_thresholds",tt$value)}
    }
    pk<-capture("polychoric descriptive KMO",psych::KMO(poly$used))
    emit("37_polychoric_descriptive_indices",data.frame(
      KMO=if(is.null(pk$value)) NA_real_ else pk$value$MSA,
      ordinal_standardized_alpha=standardized_alpha(poly$used),
      interpretation=if(poly$repaired) "Indices depend on matrix repair; not primary factorability evidence" else "Ordinal sensitivity indices",
      Bartlett_status="Not used: ordinary Pearson Bartlett inference is not a validation test for an estimated/repaired polychoric matrix"))
    pf<-capture("polychoric one-factor PAF",fa_fit(poly$used,matrix_input=TRUE,n_obs=nrow(X)))
    if(!is.null(pf$value)) {
      polyfa<-pf$value;emit("38_polychoric_PAF_one_factor",load_table(polyfa,
        if(poly$repaired) "repaired-matrix exploratory one-factor candidate" else "ordinal one-factor candidate"))
      save_model(polyfa,"polychoric_PAF_one_factor")
    } else emit_failure("38_polychoric_PAF_one_factor",pf)
    set.seed(cfg$seed+3L); sims<-matrix(NA_real_,cfg$n_pa_poly,ncol(X)); sim_audit<-vector("list",cfg$n_pa_poly)
    for(b in seq_len(cfg$n_pa_poly)) {
      Z<-vapply(seq_len(ncol(X)),function(j) sample(X[,j],nrow(X),replace=FALSE),numeric(nrow(X)))
      colnames(Z)<-colnames(X); pp<-poly_estimate(Z,paste0("polychoric PA replicate ",b))
      if(isTRUE(pp$ok)) sims[b,]<-eig_pc(pp$used)
      sim_audit[[b]]<-data.frame(replicate=b,success=isTRUE(pp$ok),
        repaired=if(isTRUE(pp$ok)) pp$repaired else NA,
        min_eigen_original=if(isTRUE(pp$ok)) pp$min_original else NA_real_,
        max_adjustment=if(isTRUE(pp$ok)) pp$max_change else NA_real_,
        error=if(!is.null(pp$error)) pp$error else "")
      if(b%%100L==0L) message("  polychoric permutations: ",b,"/",cfg$n_pa_poly)
    }
    poly_pa<-pa_summary(eig_pc(poly$used),sims,"PCA eigenvalues of estimated polychoric matrix",
      "independent within-item permutations; same estimator/repair rule; NO Pearson fallback")
    emit("39_polychoric_parallel_analysis",poly_pa,
      "Retained_sequential is withheld if too few valid replicates; matrix repair limits interpretation")
    emit("40_polychoric_parallel_replication_audit",bind_rows_base(sim_audit))
    emit("41_polychoric_parallel_null_eigenvalues",data.frame(replicate=seq_len(nrow(sims)),sims))
  } else {
    emit("35_polychoric_matrix_audit",data.frame(status="not_estimable",reason=poly$error))
    if(!is.null(poly$original)) emit("32_polychoric_matrix_UNSMOOTHED",matrix_df(poly$original))
    event("polychoric sensitivity","not_estimable",poly$error)
  }

  # ---- Q-mode: case-level profile PCA -----------------------------------------
  step("Q-mode profile analysis; distinguish relative pattern from absolute ratings")
  uniform<-apply(X,1,function(z) max(z)==min(z)); nonuniform<-X[!uniform,,drop=FALSE]
  q_ids<-d$participant_id[cc9][!uniform]
  emit("42_Qmode_profile_counts",data.frame(complete_profiles=nrow(X),uniform_profiles=sum(uniform),
    nonuniform_profiles=nrow(nonuniform),uniform_pct=100*mean(uniform),
    all_4=sum(rowSums(X==4)==ncol(X)),all_5=sum(rowSums(X==5)==ncol(X))),
    "Patient-level rating profiles, not unique HCPs. Constant profiles cannot be Pearson-correlated across items.")
  ut<-table(factor(X[uniform,1],levels=1:5))
  emit("43_Qmode_uniform_scores",data.frame(score=1:5,n=as.integer(ut)))
  qpa<-NULL
  if(nrow(nonuniform)>=3L) {
    Q<-stats::cor(t(nonuniform)); qe<-eigen(Q,symmetric=TRUE)
    rank_limit<-min(nrow(nonuniform),ncol(nonuniform)-1L)
    emit("44_Qmode_rank_audit",data.frame(profiles=nrow(nonuniform),items=ncol(nonuniform),
      rank_upper_bound=rank_limit,numerical_rank=sum(qe$values>1e-8),trace=sum(qe$values)))
    q_null<-function(null,seed) {
      set.seed(seed); sims<-matrix(NA_real_,cfg$n_pa_q,rank_limit)
      for(b in seq_len(cfg$n_pa_q)) {
        Z<-if(null=="gaussian") matrix(stats::rnorm(length(nonuniform)),nrow(nonuniform),ncol(nonuniform)) else
          t(vapply(seq_len(nrow(nonuniform)),function(i) sample(nonuniform[i,],replace=FALSE),numeric(ncol(nonuniform))))
        sims[b,]<-eig_pc(stats::cor(t(Z)))[seq_len(rank_limit)]
      }
      pa_summary(qe$values[seq_len(rank_limit)],sims,"Q-mode case-profile PCA",null)
    }
    qpa<-emit("45_Qmode_parallel_gaussian",q_null("gaussian",cfg$seed))
    qperm<-emit("46_Qmode_parallel_profile_permutation",q_null("within-profile permutation",cfg$seed+4L),
      "Preserves each profile's ordinal score distribution while disrupting item-specific preferences")
    keep<-sum(qpa$retained_sequential %in% TRUE)
    # No forced retention when zero components exceed the sequential PA threshold.
    if(keep>0L) {
      QL<-sweep(qe$vectors[,seq_len(keep),drop=FALSE],2,sqrt(qe$values[seq_len(keep)]),"*")
      qsign<-ifelse(colSums(QL)<0,-1,1); QL<-sweep(QL,2,qsign,"*")
      colnames(QL)<-paste0("QPC",seq_len(keep))
      emit("47_Qmode_profile_loadings",data.frame(participant_id=q_ids,QL,check.names=FALSE),
        "Orientation chosen so sum of loadings is positive; signs are arbitrary, not good/bad satisfaction",private=TRUE)
      emit("48_Qmode_loading_summary",bind_rows_base(lapply(seq_len(keep),function(k)
        data.frame(component=k,eigenvalue=qe$values[k],variance_pct=100*qe$values[k]/nrow(nonuniform),
          positive_loadings=sum(QL[,k]>0),negative_loadings=sum(QL[,k]<0),
          min_loading=min(QL[,k]),max_loading=max(QL[,k]),
          abs_loading_ge_0_4=sum(abs(QL[,k])>=.4),
          loading_cutoff_status="0.4 is descriptive only; not a Q-sort significance/classification rule"))))
      Zrow<-t(scale(t(nonuniform),center=TRUE,scale=TRUE))
      rel<-t(Zrow)%*%sweep(qe$vectors[,seq_len(keep),drop=FALSE],2,qsign/sqrt(qe$values[seq_len(keep)]),"*")
      colnames(rel)<-paste0("QPC",seq_len(keep))
      emit("49_Qmode_relative_component_item_scores",data.frame(code=colnames(X),rel,check.names=FALSE),
        "Standardized relative item-pattern scores; not raw satisfaction, not a classical Q-sort factor array")
    } else emit("47_Qmode_profile_loadings",data.frame(status="no_component_retained"))
    emit("50_nonuniform_profile_item_descriptives",data.frame(code=colnames(X),
      n=nrow(nonuniform),mean=colMeans(nonuniform),sd=apply(nonuniform,2,stats::sd)),
      "Unweighted means across ALL non-uniform profiles, not a factor-weighted consensus array")
  } else event("Q-mode","not_estimable","Fewer than three non-uniform profiles")

  # ---- Score construction and regression cohort -------------------------------
  step("Construct raw-item scores and ordinal regression cohort")
  scored<-d[cc9,,drop=FALSE]; scored$global_factor<-global_score
  comp<-list(fit_stability=c(2,3,9),handling_removal=c(1,7,8),wound_environment=c(4,5,6))
  scales<-list()
  for(nm in names(comp)) {
    val<-rowMeans(X[,comp[[nm]],drop=FALSE]); scored[[paste0(nm,"_mean")]]<-val
    scored[[nm]]<-as.numeric(scale(val))
    scales[[nm]]<-data.frame(domain=nm,items=paste(colnames(X)[comp[[nm]]],collapse=";"),
      mean=mean(val),sd=stats::sd(val),scoring="Arithmetic item mean, then standardization; NOT an extracted factor score")
  }
  scored$overall_impression<-scored[[raw_cols[10]]]
  emit("51_domain_scoring",bind_rows_base(scales))
  emit("52_case_level_factor_and_domain_scores",scored,
    "Global factor and conceptual-domain scores; raw responses unchanged",private=TRUE)
  reg<-scored[complete.cases(scored[,c("overall_impression","global_factor","dressing","diagnosis",names(comp))]),,drop=FALSE]
  reg$overall_ord<-ordered(reg$overall_impression,levels=sort(unique(reg$overall_impression)))
  reg$overall_top_two<-as.integer(reg$overall_impression>=4)
  # Outcome grouping is a labelled sensitivity only. All original 1-5 observations
  # remain in primary analyses and are retained in the exported analytical dataset.
  reg$overall_three<-ordered(ifelse(reg$overall_impression<=3,"1-3",as.character(reg$overall_impression)),levels=c("1-3","4","5"))
  rt<-table(factor(reg$overall_impression,levels=1:5))
  emit("53_regression_outcome_distribution",data.frame(score=1:5,n=as.integer(rt),pct=100*as.integer(rt)/nrow(reg)))
  for(v in c("dressing","diagnosis")) {
    xx<-as.data.frame(table(score=factor(reg$overall_impression,levels=1:5),group=reg[[v]]))
    names(xx)[3]<-"n"; emit(paste0("54_sparse_outcome_by_",v),xx)
  }
  emit("55_regression_cohort",data.frame(n=nrow(reg),nine_item_complete=nrow(scored),
    overall_missing_among_nine_complete=sum(is.na(scored$overall_impression)),
    scale_reference_n=nrow(X),unit="patient-level HCP assessment profile"))
  D<-as.matrix(reg[,names(comp)]); emit("56_domain_correlations",matrix_df(stats::cor(D),"domain"))
  vif<-lapply(names(comp),function(v) {
    fit<-stats::lm(stats::reformulate(setdiff(names(comp),v),response=v),data=reg)
    data.frame(domain=v,VIF=1/(1-summary(fit)$r.squared))
  })
  emit("57_domain_VIF",bind_rows_base(vif))

  # ---- Ordinal models and explicit diagnostics --------------------------------
  step("Fit ordinal models and inspect every core model, including domain models")
  or_polr<-function(f,model) {
    b<-stats::coef(f); V<-stats::vcov(f)[names(b),names(b),drop=FALSE]; se<-sqrt(diag(V))
    data.frame(model=model,n=stats::nobs(f),term=names(b),beta=unname(b),se=unname(se),
      OR=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),
      p_value=2*stats::pnorm(abs(b/se),lower.tail=FALSE),
      inference="model-based Wald; factor-score estimation uncertainty not included")
  }
  formulas<-list(global_unadjusted=overall_ord~global_factor,
    global_adjusted=overall_ord~global_factor+dressing+diagnosis,
    domains_joint=overall_ord~fit_stability+handling_removal+wound_environment,
    fit_univariable=overall_ord~fit_stability,
    handling_univariable=overall_ord~handling_removal,
    environment_univariable=overall_ord~wound_environment)
  polr_models<-list(); polr_rows<-list()
  for(nm in names(formulas)) {
    rr<-capture(paste0("polr ",nm),MASS::polr(formulas[[nm]],data=reg,Hess=TRUE,
      method="logistic",control=list(maxit=2000,reltol=1e-10)))
    if(!is.null(rr$value)) {
      polr_models[[nm]]<-rr$value
      et<-capture(paste0("polr ",nm," coefficient extraction"),or_polr(rr$value,nm))
      polr_rows[[nm]]<-if(!is.null(et$value)) et$value else data.frame(model=nm,status="coefficient extraction failed",reason=et$error)
      save_model(rr$value,paste0("polr_",nm))
    } else polr_rows[[nm]]<-data.frame(model=nm,status="failed",reason=rr$error)
  }
  emit("58_ordinal_primary_and_domain_models",bind_rows_base(polr_rows))

  clm_diagnostics<-list(); clm_models<-list(); clm_param_rows<-list()
  convergence_code<-function(f) {
    x<-f$convergence
    if(is.list(x)) {
      if(!is.null(x$code)) return(scalar(x$code))
      if(!is.null(x$convergence)) return(scalar(x$convergence))
      return(NA_real_)
    }
    scalar(x)
  }
  clm_check<-function(f,label,dat,condition_text="") {
    H<-f$Hessian; ee<-tryCatch(eigen(H,symmetric=TRUE,only.values=TRUE)$values,error=function(e) NA_real_)
    g<-if(length(f$gradient)) max(abs(f$gradient)) else scalar(f$maxGradient)
    mn<-if(all(is.finite(ee))) min(ee) else NA_real_
    cn<-if(is.finite(mn)&&mn>0) max(ee)/mn else Inf
    code<-convergence_code(f); aliases<-any(unlist(f$aliased),na.rm=TRUE)
    resp<-all.vars(stats::formula(f))[1]
    pred<-capture(paste0(label," probability check"),stats::predict(f,
      newdata=dat[,setdiff(names(dat),resp),drop=FALSE],type="prob")$fit)
    probs<-pred$value
    pvalid<-is.matrix(probs)&&nrow(probs)==nrow(dat)&&all(is.finite(probs))&&
      all(probs>=0)&&all(probs<=1)&&max(abs(rowSums(probs)-1))<1e-6
    good<-is.finite(code)&&code==0&&is.finite(g)&&g<cfg$diagnostic_gradient_tolerance&&
      is.finite(mn)&&mn>0&&is.finite(cn)&&cn<cfg$diagnostic_condition_limit&&!aliases&&pvalid
    data.frame(model=label,n=stats::nobs(f),convergence_code=code,
      max_abs_gradient=g,min_Hessian_eigen=mn,Hessian_condition=cn,
      aliased_parameters=aliases,valid_training_probabilities=pvalid,
      logLik=as.numeric(stats::logLik(f)),parameters=attr(stats::logLik(f),"df"),
      AIC=stats::AIC(f),usable_for_asymptotic_LRT=good,
      status=if(good) "numerically_acceptable_not_proof_of_assumptions" else "inconclusive_or_unstable",
      conditions=condition_text)
  }
  clm_parameters<-function(f,label,scale_present=FALSE) {
    # Location and scale effects can have IDENTICAL names (e.g. diagnosisDFU).
    # Work positionally in alpha/beta/zeta order; name-based matching would select
    # the wrong diagnosis coefficient in a location-scale model.
    b<-c(f$alpha,f$beta,f$zeta)
    type<-c(rep("threshold_or_nominal",length(f$alpha)),rep("location",length(f$beta)),
            rep("log_scale",length(f$zeta)))
    V<-tryCatch(stats::vcov(f),error=function(e) NULL)
    se<-rep(NA_real_,length(b))
    if(!is.null(V)&&nrow(V)==length(b)) se<-sqrt(diag(V))
    meaningful_or<-type=="location" & !grepl("ns\\(",names(b))
    data.frame(model=label,parameter_index=seq_along(b),term=names(b),parameter_type=type,
      estimate=unname(b),SE=unname(se),lower95=unname(b-1.96*se),upper95=unname(b+1.96*se),
      OR_if_location=ifelse(meaningful_or,exp(b),NA_real_),
      OR_low95=ifelse(meaningful_or,exp(b-1.96*se),NA_real_),
      OR_high95=ifelse(meaningful_or,exp(b+1.96*se),NA_real_),
      p_value=2*stats::pnorm(abs(b/se),lower.tail=FALSE),
      OR_interpretation=if(scale_present) "Reference latent-scale OR only; not common across scale strata" else "Cumulative OR only for location parameters")
  }
  clm_fit<-function(form,dat,label,nominal=NULL,scale=NULL,record=TRUE) {
    # do.call uses evaluated formulas and a data frame: update-based diagnostics
    # do not depend on vanished helper-local names for formula/data arguments.
    args<-list(formula=form,data=dat,link="logit",threshold="flexible",model=TRUE,
      control=ordinal::clm.control(maxIter=200,maxLineIter=100,gradTol=1e-7))
    if(!is.null(nominal)) args$nominal<-nominal
    if(!is.null(scale)) args$scale<-scale
    rr<-capture(label,do.call(ordinal::clm,args))
    if(is.null(rr$value)) {
      dd<-data.frame(model=label,n=nrow(dat),status="fit_failed",
        usable_for_asymptotic_LRT=FALSE,conditions=rr$error)
      if(record) clm_diagnostics[[label]]<<-dd
      return(list(fit=NULL,diagnostics=dd))
    }
    f<-rr$value; dd<-clm_check(f,label,dat,paste(c(rr$warnings,rr$messages),collapse=" | "))
    if(record) {
      clm_diagnostics[[label]]<<-dd;clm_models[[label]]<<-f
      pp<-capture(paste0(label," parameter extraction"),clm_parameters(f,label,!is.null(scale)))
      clm_param_rows[[label]]<<-if(!is.null(pp$value)) pp$value else
        data.frame(model=label,status="parameter extraction failed",reason=pp$error)
      save_model(f,label)
      cv<-capture(paste0(label," convergence detail"),ordinal::convergence(f))
      utils::capture.output(print(cv$value),file=file.path(out,"models",paste0(label,"_convergence.txt")))
    }
    list(fit=f,diagnostics=dd)
  }
  lrt<-function(base,alternative,kind,tested_term) {
    ok<-!is.null(base$fit)&&!is.null(alternative$fit)&&
      isTRUE(base$diagnostics$usable_for_asymptotic_LRT)&&isTRUE(alternative$diagnostics$usable_for_asymptotic_LRT)
    lr<-df<-pv<-NA_real_; status<-"inconclusive: alternative/base did not meet numerical checks"
    if(ok) {
      df<-attr(stats::logLik(alternative$fit),"df")-attr(stats::logLik(base$fit),"df")
      lr<-2*(as.numeric(stats::logLik(alternative$fit))-as.numeric(stats::logLik(base$fit)))
      if(df>0&&lr>=-1e-6&&stats::nobs(base$fit)==stats::nobs(alternative$fit)) {
        lr<-max(lr,0);pv<-stats::pchisq(lr,df,lower.tail=FALSE)
        status<-"computed; exploratory asymptotic LRT with sparse ordinal categories"
      } else {pv<-NA_real_;status<-"inconclusive: nonnested/different-N/numerically inferior alternative"}
    }
    data.frame(base_model=base$diagnostics$model,alternative_model=alternative$diagnostics$model,
      test=kind,term=tested_term,LR=lr,df=df,p_value=pv,status=status)
  }
  model_formula<-function(terms,response="overall_ord") {
    if(length(terms)) stats::reformulate(terms,response=response) else stats::as.formula(paste(response,"~1"))
  }
  clm_base<-list(); tests<-list(); official_tests<-list()
  for(nm in names(formulas)) {
    base<-clm_fit(formulas[[nm]],reg,paste0("clm_",nm));clm_base[[nm]]<-base
    if(is.null(base$fit)) next
    for(fun in c("nominal_test","scale_test")) {
      rr<-capture(paste0(nm," official ",fun),getExportedValue("ordinal",fun)(base$fit))
      if(!is.null(rr$value)) {tab<-as.data.frame(rr$value);tab$term<-rownames(tab);tab$model<-nm;tab$test<-fun
        official_tests[[paste(nm,fun)]]<-tab}
      else official_tests[[paste(nm,fun)]]<-data.frame(model=nm,test=fun,status="failed",reason=rr$error)
    }
    terms<-attr(stats::terms(formulas[[nm]]),"term.labels")
    for(tt in terms) {
      # Exclude tt from location when giving it threshold-specific nominal effects;
      # this avoids a redundant location/nominal parameterization.
      nf<-stats::reformulate(tt)
      alt<-clm_fit(model_formula(setdiff(terms,tt)),reg,paste0("nominal_",nm,"_",tt),nominal=nf)
      tests[[paste0("nom_",nm,"_",tt)]]<-lrt(base,alt,"partial proportional-odds nominal extension",tt)
      sc<-clm_fit(formulas[[nm]],reg,paste0("scale_",nm,"_",tt),scale=stats::reformulate(tt))
      tests[[paste0("scale_",nm,"_",tt)]]<-lrt(base,sc,"latent-scale extension (not identical to nominal test)",tt)
    }
    if(length(terms)>1L) {
      full<-clm_fit(model_formula(character()),reg,paste0("nominal_omnibus_",nm),nominal=stats::reformulate(terms))
      tests[[paste0("nominal_omnibus_",nm)]]<-lrt(base,full,"omnibus non-proportional-odds extension",paste(terms,collapse="+"))
    }
  }
  emit("59_official_nominal_scale_tests",bind_rows_base(official_tests),
    "NA is inconclusive, not evidence of proportional odds; see explicit model refits and convergence logs")
  # Nonlinearity and sparse-category sensitivities supplement, not replace, 1-5 analysis.
  nl<-clm_fit(overall_ord~splines::ns(global_factor,df=3)+dressing+diagnosis,reg,"clm_global_spline")
  tests[["nonlinearity"]]<-lrt(clm_base$global_adjusted,nl,"natural-spline nonlinearity sensitivity","global_factor")
  three<-clm_fit(overall_three~global_factor+dressing+diagnosis,reg,"clm_three_category_sensitivity")
  if(!is.null(three$fit)) {
    nominal3<-clm_fit(overall_three~dressing+diagnosis,reg,"clm_three_category_nominal_global",nominal=~global_factor)
    tests[["grouped_PO"]]<-lrt(three,nominal3,"nominal extension after grouping 1-3; sensitivity only","global_factor")
  }
  emit("60_explicit_ordinal_diagnostic_LRTs",bind_rows_base(tests))
  emit("61_all_ordinal_model_numerical_diagnostics",bind_rows_base(clm_diagnostics))
  emit("62_all_clm_parameters",bind_rows_base(clm_param_rows),
    "Nominal/threshold estimates are not exponentiated as common cumulative ORs; scale-model ORs require stratum-specific interpretation")

  # Expose predicted probabilities and in-sample calibration for acceptable models.
  prediction_rows<-list(); calibration_rows<-list()
  pred_models<-intersect(c("clm_global_adjusted","scale_global_adjusted_diagnosis",
                           "clm_domains_joint","clm_three_category_sensitivity"),names(clm_models))
  for(nm in pred_models) {
    f<-clm_models[[nm]];dg<-clm_diagnostics[[nm]]
    if(!isTRUE(dg$usable_for_asymptotic_LRT)) next
    response<-all.vars(stats::formula(f))[1]
    pp<-capture(paste0(nm," predictions"),stats::predict(f,newdata=reg[,setdiff(names(reg),response),drop=FALSE],type="prob")$fit)
    if(is.null(pp$value)) next
    pr<-pp$value; yy<-as.integer(reg[[response]]); K<-ncol(pr)
    for(k in seq_len(K)) prediction_rows[[paste(nm,k)]]<-data.frame(model=nm,
      participant_id=reg$participant_id,category=colnames(pr)[k],probability=pr[,k])
    for(k in seq_len(K-1L)) {
      cp<-rowSums(pr[,seq_len(k),drop=FALSE]); obs<-as.integer(yy<=k)
      br<-unique(stats::quantile(cp,seq(0,1,length.out=6),names=FALSE))
      bin<-if(length(br)>1L) cut(cp,breaks=br,include.lowest=TRUE,labels=FALSE) else rep(1L,length(cp))
      for(b in sort(unique(bin))) {ix<-which(bin==b)
        calibration_rows[[paste(nm,k,b)]]<-data.frame(model=nm,cumulative_cut=colnames(pr)[k],
          bin=b,n=length(ix),mean_predicted=mean(cp[ix]),observed_fraction=mean(obs[ix]))}
    }
  }
  emit("63_ordinal_predicted_probabilities",bind_rows_base(prediction_rows),private=TRUE)
  emit("64_ordinal_apparent_calibration",bind_rows_base(calibration_rows),
    "In-sample calibration only, not out-of-sample predictive validation")

  # Diagnosis-scale model: location coefficient is on the reference VLU scale.
  # Group-specific global-score log OR is beta/exp(zeta*I(DFU)); not an interaction test.
  sname<-"scale_global_adjusted_diagnosis"
  if(sname%in%names(clm_models)&&isTRUE(clm_diagnostics[[sname]]$usable_for_asymptotic_LRT)) {
    f<-clm_models[[sname]];vc<-stats::vcov(f)
    bpos<-match("global_factor",names(f$beta));zpos<-match("diagnosisDFU",names(f$zeta))
    if(!is.na(bpos)&&!is.na(zpos)&&nrow(vc)==length(c(f$alpha,f$beta,f$zeta))) {
      z<-unname(f$zeta[zpos]);b<-unname(f$beta[bpos]);ans<-list()
      indices<-c(length(f$alpha)+bpos,length(f$alpha)+length(f$beta)+zpos)
      for(dx in c("VLU","DFU")) {
        iz<-as.numeric(dx=="DFU"); ell<-b*exp(-iz*z)
        gradient<-c(exp(-iz*z),-iz*ell)
        se<-sqrt(as.numeric(t(gradient)%*%vc[indices,indices,drop=FALSE]%*%gradient))
        ans[[dx]]<-data.frame(diagnosis=dx,global_score_logOR=ell,SE=se,OR=exp(ell),
          CI_low=exp(ell-1.96*se),CI_high=exp(ell+1.96*se),
          interpretation="Diagnosis-scale sensitivity only; strongly model-dependent, not causal effect modification")
      }
      emit("65_diagnosis_scale_model_global_association",bind_rows_base(ans))
    } else event(sname,"information","Could not unambiguously locate diagnosis scale coefficient; inspect saved parameter table")
  }

  # Threshold-specific binary regressions diagnose instability; they are not a
  # joint ordinal model and overlapping cut-point comparisons are not independent.
  glm_table<-function(f,model,vc=NULL,df=Inf) {
    b<-stats::coef(f);if(is.null(vc)) vc<-stats::vcov(f)
    se<-sqrt(diag(vc));crit<-if(is.finite(df)) stats::qt(.975,df) else stats::qnorm(.975)
    p<-if(is.finite(df)) 2*stats::pt(abs(b/se),df,lower.tail=FALSE) else 2*stats::pnorm(abs(b/se),lower.tail=FALSE)
    data.frame(model=model,n=stats::nobs(f),term=names(b),beta=unname(b),SE=unname(se),
      OR=exp(b),CI_low=exp(b-crit*se),CI_high=exp(b+crit*se),p_value=p,
      reference_df=df,converged=isTRUE(f$converged))
  }
  threshold_rows<-list()
  for(k in 1:4) {
    dd<-reg;dd$above_cut<-as.integer(dd$overall_impression>k)
    if(length(unique(dd$above_cut))<2L) next
    for(type in c("ML","mean_bias_reduced")) {
      ff<-if(type=="ML") capture(paste0("threshold ",k," ML"),stats::glm(above_cut~global_factor+dressing+diagnosis,data=dd,family=stats::binomial())) else
        capture(paste0("threshold ",k," mean BR"),stats::glm(above_cut~global_factor+dressing+diagnosis,
          data=dd,family=stats::binomial(),method=brglm2::brglmFit,type="AS_mean"))
      if(!is.null(ff$value)) {
        tab<-glm_table(ff$value,paste0("score > ",k,"; ",type));tab$lower_n<-sum(dd$above_cut==0);tab$higher_n<-sum(dd$above_cut==1)
        threshold_rows[[paste(k,type)]]<-tab
      }
    }
  }
  emit("66_threshold_specific_binary_diagnostics",bind_rows_base(threshold_rows),
    "Sensitivity diagnostics, not four independent tests of proportional odds; bias reduction helps sparse binary fits but does not validate a cumulative model")

  # ---- OLS restored, top-two binary sensitivity and clustering -----------------
  step("OLS, binary sensitivities, clustering and leave-one-site-out influence")
  ols_rows<-list()
  for(adj in c(FALSE,TRUE)) {
    ff<-if(adj) overall_impression~global_factor+dressing+diagnosis else overall_impression~global_factor
    lmfit<-stats::lm(ff,data=reg);save_model(lmfit,paste0("OLS_",adj))
    for(vtype in c("model_based","HC3")) {
      V<-if(vtype=="HC3") sandwich::vcovHC(lmfit,type="HC3") else stats::vcov(lmfit)
      b<-stats::coef(lmfit);se<-sqrt(diag(V));df<-stats::df.residual(lmfit);crit<-stats::qt(.975,df)
      ols_rows[[paste(adj,vtype)]]<-data.frame(model=if(adj) "OLS adjusted" else "OLS unadjusted",
        variance=vtype,n=stats::nobs(lmfit),term=names(b),beta=unname(b),SE=unname(se),
        CI_low=b-crit*se,CI_high=b+crit*se,p_value=2*stats::pt(abs(b/se),df,lower.tail=FALSE),
        R_squared=summary(lmfit)$r.squared,df=df,
        interpretation="Descriptive numeric-score sensitivity; ordinal scale not assumed interval in primary analysis")
    }
  }
  emit("67_OLS_sensitivity",bind_rows_base(ols_rows))
  tb<-capture("top-two binary ML",stats::glm(overall_top_two~global_factor+dressing+diagnosis,data=reg,family=stats::binomial()))
  tbbr<-capture("top-two binary BR",stats::glm(overall_top_two~global_factor+dressing+diagnosis,data=reg,
    family=stats::binomial(),method=brglm2::brglmFit,type="AS_mean"))
  if(!is.null(tb$value)) emit("68_top_two_binary_ML",glm_table(tb$value,"top-two adjusted")) else emit_failure("68_top_two_binary_ML",tb)
  if(!is.null(tbbr$value)) emit("69_top_two_binary_bias_reduced",glm_table(tbbr$value,"top-two adjusted mean BR")) else emit_failure("69_top_two_binary_bias_reduced",tbbr)
  cluster_rows<-list(); cluster_counts<-list()
  for(cv in c("site_prefix","ecrf_user_id")) {
    dd<-reg[!is.na(reg[[cv]])&nzchar(reg[[cv]]),,drop=FALSE];G<-length(unique(dd[[cv]]))
    tt<-table(dd[[cv]]);cluster_counts[[cv]]<-data.frame(cluster_type=rep(cv,length(tt)),cluster=names(tt),n=as.integer(tt))
    if(G<2L) {event(cv,"not_estimable","Fewer than two clusters");next}
    gf<-capture(paste0("cluster binary ",cv),stats::glm(overall_top_two~global_factor+dressing+diagnosis,data=dd,family=stats::binomial()))
    if(!is.null(gf$value)) {
      vv<-capture(paste0("cluster binary covariance ",cv),sandwich::vcovCL(gf$value,cluster=dd[[cv]],type="HC1",cadjust=TRUE))
      if(!is.null(vv$value)) for(ref in c("normal","t_G_minus_1")) {
        z<-glm_table(gf$value,paste0("binary cluster ",cv),vv$value,if(ref=="normal") Inf else G-1L)
        z$clusters<-G;z$cluster_type<-cv;z$reference_distribution<-ref
        cluster_rows[[paste(cv,"binary",ref)]]<-z
      }
    }
    dd$.analysis_weight<-1
    sf<-capture(paste0("survey ordinal ",cv),{
      des<-survey::svydesign(ids=stats::reformulate(cv),weights=~.analysis_weight,data=dd)
      survey::svyolr(overall_ord~global_factor+dressing+diagnosis,design=des,method="logistic")
    })
    if(!is.null(sf$value)) {
      f<-sf$value;b_all<-stats::coef(f)
      wanted<-intersect(c("global_factor","dressingRSSil","diagnosisDFU"),names(b_all))
      b<-b_all[wanted];V<-stats::vcov(f)[wanted,wanted,drop=FALSE];se<-sqrt(diag(V))
      for(ref in c("normal","t_G_minus_1")) {
        df<-if(ref=="normal") Inf else G-1L;crit<-if(is.finite(df)) stats::qt(.975,df) else stats::qnorm(.975)
        z<-data.frame(model=paste0("ordinal cluster ",cv),n=nrow(dd),term=names(b),
          beta=unname(b),SE=unname(se),OR=exp(b),CI_low=exp(b-crit*se),CI_high=exp(b+crit*se),
          p_value=if(is.finite(df)) 2*stats::pt(abs(b/se),df,lower.tail=FALSE) else 2*stats::pnorm(abs(b/se),lower.tail=FALSE),
          reference_df=df,clusters=G,cluster_type=cv,reference_distribution=ref)
        cluster_rows[[paste(cv,"ordinal",ref)]]<-z
      }
      save_model(f,paste0("survey_ordinal_",cv))
    }
  }
  emit("70_cluster_sensitivity_models",bind_rows_base(cluster_rows),
    "Eight sites/ten entry accounts in supplied data. t(G-1) is a sensitivity, not a guaranteed small-cluster correction. Account is not verified HCP.")
  emit("71_cluster_sizes",bind_rows_base(cluster_counts))
  loso<-list()
  for(site in sort(unique(reg$site_prefix[!is.na(reg$site_prefix)]))) {
    dd<-reg[reg$site_prefix!=site&!is.na(reg$site_prefix),,drop=FALSE]
    # Fixed original factor score: isolates site influence on regression, not EFA.
    ff<-clm_fit(overall_ord~global_factor+dressing+diagnosis,dd,
                 paste0("leave_site_",gsub("[^A-Za-z0-9]","",site)),record=FALSE)
    if(!is.null(ff$fit)&&isTRUE(ff$diagnostics$usable_for_asymptotic_LRT)) {
      b<-ff$fit$beta["global_factor"];se<-sqrt(stats::vcov(ff$fit)["global_factor","global_factor"])
      loso[[site]]<-data.frame(omitted_site=site,n=nrow(dd),OR=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),status="computed")
    } else loso[[site]]<-data.frame(omitted_site=site,n=nrow(dd),OR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,status=ff$diagnostics$status)
  }
  emit("72_leave_one_site_out",bind_rows_base(loso),
    "Regression influence sensitivity with full-sample factor score held fixed; not independent factor validation")

  # ---- Full-pipeline stability: resample -> PAF -> score -> ordinal model --------
  step("Bootstrap stability, including factor-score re-estimation")
  boot_one<-function(dd,label) {
    ans<-data.frame(replicate=NA_integer_,scheme=label,n=nrow(dd),success=FALSE,
      alpha=NA_real_,global_beta=NA_real_,OR=NA_real_,reason="")
    A<-as.matrix(dd[,domain_cols,drop=FALSE]);colnames(A)<-item_codes[1:9]
    ans$alpha<-alpha_raw(A)
    if(any(apply(A,2,stats::sd)==0)) {ans$reason<-"constant resampled item";return(ans)}
    ef<-capture(paste0("bootstrap ",label," PAF"),fa_fit(A))
    if(is.null(ef$value)) {ans$reason<-ef$error;return(ans)}
    ll<-as.numeric(unclass(ef$value$loadings)[,1]);sgn<-if(sum(ll)<0) -1 else 1;ll<-ll*sgn
    for(j in seq_along(ll)) ans[[paste0("loading_",item_codes[j])]]<-ll[j]
    dd$global_factor<-as.numeric(scale(as.numeric(ef$value$scores[,1])*sgn))
    dd<-dd[is.finite(dd$global_factor)&!is.na(dd$overall_impression),,drop=FALSE]
    dd$overall_ord<-ordered(dd$overall_impression,levels=sort(unique(dd$overall_impression)))
    if(length(unique(dd$dressing))<2||length(unique(dd$diagnosis))<2||nlevels(dd$overall_ord)<3) {
      ans$reason<-"insufficient response/covariate variation";return(ans)}
    ff<-capture(paste0("bootstrap ",label," ordinal"),MASS::polr(overall_ord~global_factor+dressing+diagnosis,
      data=dd,Hess=TRUE,method="logistic",control=list(maxit=2000,reltol=1e-10)))
    if(is.null(ff$value)) {ans$reason<-ff$error;return(ans)}
    f<-ff$value;eh<-tryCatch(eigen(f$Hessian,symmetric=TRUE,only.values=TRUE)$values,error=function(e) NA_real_)
    okay<-isTRUE(f$convergence==0)&&all(is.finite(eh))&&min(eh)>0&&max(eh)/min(eh)<cfg$diagnostic_condition_limit&&
      is.finite(stats::coef(f)["global_factor"])
    ans$n<-nrow(dd);ans$success<-okay
    if(okay) {ans$global_beta<-stats::coef(f)["global_factor"];ans$OR<-exp(ans$global_beta)} else ans$reason<-"nonconvergence or ill-conditioned ordinal Hessian"
    ans
  }
  boot_rows<-list(); count<-0L
  for(scheme in c("patient_iid","site_cluster")) {
    B<-if(scheme=="patient_iid") cfg$n_boot else cfg$n_boot_site
    base<-scored
    if(scheme=="site_cluster") base<-base[!is.na(base$site_prefix)&nzchar(base$site_prefix),,drop=FALSE]
    sites<-unique(base$site_prefix)
    if(scheme=="site_cluster"&&length(sites)<2L) {event("site bootstrap","not_estimable","Insufficient sites");next}
    set.seed(cfg$seed+if(scheme=="patient_iid") 10L else 11L)
    for(b in seq_len(B)) {
      ii<-if(scheme=="patient_iid") sample.int(nrow(base),nrow(base),replace=TRUE) else
        unlist(lapply(sample(sites,length(sites),replace=TRUE),function(g) which(base$site_prefix==g)),use.names=FALSE)
      rr<-boot_one(base[ii,,drop=FALSE],paste0(scheme,"_",b));rr$scheme<-scheme;rr$replicate<-b
      count<-count+1L;boot_rows[[count]]<-rr
      if(b%%100L==0L) message("  ",scheme," bootstrap: ",b,"/",B)
    }
  }
  boot<-emit("73_full_pipeline_bootstrap_replicates",bind_rows_base(boot_rows),
    "Factor model and scale re-estimated each replicate. Patient IID ignores clustering; site bootstrap has few clusters. Failed replicates retained in audit.")
  bsumm<-list()
  if(nrow(boot)) for(scheme in unique(boot$scheme)) {
    z<-boot[boot$scheme==scheme,,drop=FALSE]
    for(v in intersect(c("alpha","global_beta","OR",paste0("loading_",item_codes[1:9])),names(z))) {
      valid<-is.finite(z[[v]])
      if(v%in%c("global_beta","OR")) valid<-valid & z$success
      fraction<-mean(valid);quant<-if(sum(valid)>1&&fraction>=cfg$min_sim_success)
        stats::quantile(z[[v]][valid],c(.025,.5,.975),names=FALSE) else rep(NA_real_,3)
      bsumm[[paste(scheme,v)]]<-data.frame(scheme=scheme,statistic=v,requested=nrow(z),
        successful=sum(valid),success_fraction=fraction,lower025=quant[1],median=quant[2],upper975=quant[3],
        interval_status=if(fraction>=cfg$min_sim_success) "percentile stability interval" else "withheld: too many non-estimable replicates")
    }
  }
  emit("74_full_pipeline_bootstrap_summary",bind_rows_base(bsumm))

  # ---- Clinical plausibility: exploratory associations, not causal validation --
  step("Clinical associations and indication-specific structural sensitivity")
  hypotheses<-data.frame(outcome=outcomes,
    expected_direction=c("positive","positive","positive","positive","negative","negative","negative","unspecified","unspecified"),
    note=c(rep("Clinically motivated expectation; not claimed prospectively preregistered",7),
           "Exploratory care-burden association","Exploratory care-burden association"))
  corr_one<-function(dd,outcome,label) {
    if(!outcome%in%names(dd)) return(data.frame(analysis=label,outcome=outcome,n=0,status="variable unavailable"))
    z<-dd[complete.cases(dd[,c("global_factor",outcome)]),c("global_factor",outcome),drop=FALSE]
    if(nrow(z)<5||stats::sd(z[[outcome]])==0||stats::sd(z$global_factor)==0)
      return(data.frame(analysis=label,outcome=outcome,n=nrow(z),status="insufficient variation"))
    rr<-capture(paste0(label," Spearman ",outcome),stats::cor.test(z$global_factor,z[[outcome]],method="spearman",exact=FALSE))
    if(is.null(rr$value)) return(data.frame(analysis=label,outcome=outcome,n=nrow(z),status=rr$error))
    set.seed(cfg$seed+200L+match(outcome,c(outcomes,"overall_impression")))
    vals<-rep(NA_real_,cfg$n_boot_cor)
    for(b in seq_len(cfg$n_boot_cor)) {
      ii<-sample.int(nrow(z),nrow(z),replace=TRUE)
      if(stats::sd(z$global_factor[ii])>0&&stats::sd(z[[outcome]][ii])>0)
        vals[b]<-stats::cor(z$global_factor[ii],z[[outcome]][ii],method="spearman")
    }
    valid<-is.finite(vals);qc<-if(mean(valid)>=cfg$min_sim_success) stats::quantile(vals[valid],c(.025,.975),names=FALSE) else c(NA_real_,NA_real_)
    data.frame(analysis=label,outcome=outcome,n=nrow(z),rho=unname(rr$value$estimate),p_value=rr$value$p.value,
      CI_low=qc[1],CI_high=qc[2],bootstrap_valid=sum(valid),status="computed",
      CI_interpretation="IID bootstrap with original factor score fixed; ignores rater/site clustering")
  }
  ext<-bind_rows_base(lapply(outcomes,function(v) corr_one(scored,v,"pooled raw")))
  ext$expected_direction<-hypotheses$expected_direction[match(ext$outcome,hypotheses$outcome)]
  ext$hypothesis_status<-hypotheses$note[match(ext$outcome,hypotheses$outcome)]
  ext$p_BH_sensitivity<-stats::p.adjust(ext$p_value,method="BH")
  emit("75_clinical_plausibility_associations",ext,
    "Main p-values unadjusted/exploratory; BH column supplementary. Direction consistency alone does not prove construct validity.")
  emit("76_global_judgement_Spearman",corr_one(scored,"overall_impression","within-instrument association"),
    "Overall impression is not an external gold standard or criterion-validity standard")
  substats<-list(); subloads<-list(); subcors<-list(); subpa<-list()
  for(dx in levels(d$diagnosis)) {
    ii<-which(d$diagnosis==dx & cc9);A<-as.matrix(d[ii,domain_cols]);colnames(A)<-item_codes[1:9]
    if(nrow(A)<=ncol(A)+2L||any(apply(A,2,stats::sd)==0)) {
      substats[[dx]]<-data.frame(analysis_set=dx,n=nrow(A),status="insufficient variation/sample");next}
    substats[[dx]]<-factorability(A,dx)
    ff<-capture(paste0("subgroup ",dx," PAF"),fa_fit(A))
    if(!is.null(ff$value)) {z<-load_table(ff$value,paste0(dx," one-factor candidate"));z$diagnosis<-dx;subloads[[dx]]<-z}
    subpa[[paste0(dx,"_PCA")]]<-transform(pa_run(A,cfg$n_pa,cfg$seed,"gaussian"),diagnosis=dx)
    pp<-capture(paste0("subgroup ",dx," common PA"),pa_run(A,cfg$n_pa_common,cfg$seed+2L,"permutation",TRUE))
    if(!is.null(pp$value)) subpa[[paste0(dx,"_common")]]<-transform(pp$value,diagnosis=dx)
    ss<-scored[scored$diagnosis==dx,,drop=FALSE]
    subcors[[dx]]<-bind_rows_base(lapply(c("war_pre","war_post","final_exudate_volume"),
                                      function(v) corr_one(ss,v,paste0(dx," using pooled global score"))))
  }
  emit("77_subgroup_factorability",bind_rows_base(substats),"Indication-specific exploratory results; not measurement invariance")
  emit("78_subgroup_one_factor_loadings",bind_rows_base(subloads))
  emit("79_subgroup_parallel_analysis",bind_rows_base(subpa),"PCA and reduced-SMC common-factor diagnostics are separate analyses")
  emit("80_subgroup_clinical_associations",bind_rows_base(subcors),"Associations use the pooled global score; subgroup EFA scores are not substituted")

  # ---- Programmatic plots; no hand-entered statistical labels ------------------
  if(isTRUE(cfg$make_figures)) {
    step("Create numerical figures directly from this run")
    plot_file<-function(name,draw,width=2400,height=1600) {
      ff<-file.path(out,"figures",paste0(name,".png"))
      pp<-capture(paste0("figure ",name),{
        grDevices::png(ff,width=width,height=height,res=300)
        tryCatch(force(draw),finally=grDevices::dev.off())
      })
      invisible(pp$value)
    }
    plot_file("Figure1_raw_satisfaction",{
      graphics::par(mar=c(5,14,3,7));yy<-rev(seq_len(nrow(desc)))
      graphics::plot(desc$mean,yy,xlim=c(1,5.55),ylim=c(.5,10.7),yaxt="n",xlab="Mean satisfaction score (1-5)",ylab="",pch=16,bty="l")
      graphics::axis(2,at=yy,labels=item_names,las=2,cex.axis=.66)
      graphics::text(5.08,yy,sprintf("%.2f; %.1f%%; n=%d",desc$mean,desc$top_two_pct,desc$valid_n),pos=4,cex=.62,xpd=TRUE)
      graphics::mtext("Mean; scores 4-5; valid n",side=3,at=5.1,cex=.6)
    })
    plot_file("Figure2_primary_factor_loadings",{
      graphics::par(mar=c(5,14,3,3));yy<-rev(seq_len(9))
      graphics::plot(L1,yy,xlim=c(0,1),ylim=c(.5,9.5),yaxt="n",xlab="One-factor principal-axis loading",ylab="",pch=16,bty="l")
      graphics::axis(2,at=yy,labels=item_names[1:9],las=2,cex.axis=.65)
      graphics::text(L1,yy,sprintf("%.3f",L1),pos=4,cex=.75)
    })
    plot_file("FigureS1_PCA_parallel_analysis",{
      graphics::plot(pa9$component,pa9$observed_eigen,type="b",pch=16,xlab="Component",ylab="Eigenvalue",bty="l")
      graphics::lines(pa9$component,pa9$random_p95,type="b",pch=1,lty=2)
      graphics::abline(h=1,lty=3)
      graphics::legend("topright",legend=c("Observed PCA","Gaussian PA 95th percentile","Kaiser = 1"),lty=c(1,2,3),pch=c(16,1,NA),bty="n")
    })
    if(!is.null(qpa)) plot_file("FigureS2_Qmode_parallel_analysis",{
      graphics::plot(qpa$component,qpa$observed_eigen,type="b",pch=16,xlab="Q-mode component",ylab="Eigenvalue",bty="l")
      graphics::lines(qpa$component,qpa$random_p95,type="b",pch=1,lty=2)
      graphics::legend("topright",legend=c("Observed case-profile PCA","Gaussian PA 95th percentile"),lty=c(1,2),pch=c(16,1),bty="n")
    })
    if(length(polr_rows)) {
      fr<-bind_rows_base(polr_rows);fr<-fr[fr$term%in%c("global_factor",names(comp)) &
        fr$model%in%c("global_unadjusted","global_adjusted","domains_joint") & is.finite(fr$OR),,drop=FALSE]
      if(nrow(fr)) plot_file("Figure3_ordinal_associations",{
        graphics::par(mar=c(5,12,3,3));yy<-rev(seq_len(nrow(fr)))
        graphics::plot(fr$OR,yy,log="x",xlim=range(c(fr$CI_low,fr$CI_high)),
          ylim=c(.5,nrow(fr)+.5),pch=16,yaxt="n",xlab="Cumulative odds ratio (95% CI), per SD",ylab="",bty="l")
        graphics::segments(fr$CI_low,yy,fr$CI_high,yy);graphics::abline(v=1,lty=2)
        graphics::axis(2,at=yy,labels=paste(fr$model,fr$term,sep=": "),las=2,cex.axis=.65)
      })
    }
  }
  # Include author-facing cautions in the same consolidated file.
  cautions<-data.frame(topic=c("Scoring","Measurement","Primary dimensionality","Polychoric","Ordinal diagnostics",
    "Multiplicity","Clustering","Q-mode","Exploratory expectations","Source data"),
    interpretation=c(
      "All legitimate 1-5 responses retained; no report-based low-score recoding implemented.",
      "This is an initial quantitative evaluation; content validity, reproducibility and responsiveness are not established.",
      "PCA parallel analysis and common-factor sensitivity are distinct; do not call PCA variance EFA explained variance or force agreement.",
      "Inspect unsmoothed matrix, repairs and replicate failures before interpreting ordinal sensitivity; no Pearson fallback.",
      "Check diagnostic status for each model; NA/failure is inconclusive. Diagnosis scale model does not have a common score OR across indications.",
      "All added diagnostics and sensitivities exploratory; no automatic model selection by p value.",
      "Few site/entry-account clusters; account IDs do not identify verified HCP raters. Bootstrap/t-reference sensitivities do not remove this limitation.",
      "Patient-profile PCA examines relative item-pattern similarity, not absolute satisfaction or confirmed HCP types.",
      "Clinical directions are interpretive expectations, not evidence of prospective prespecification.",
      "Analytical CSV consistency checks do not substitute for source-document/clinical database audit."))
  emit("97_interpretation_and_reporting_cautions",cautions)
  audit$run_status<-"analysis_completed_with_diagnostic_statuses_to_review"
  message("Analysis sections completed. Review all diagnostic/status tables before manuscript use.")
  invisible(out)
}

# ---- Execution entry point ----------------------------------------------------
# Options override defaults for RStudio source() without changing working directory:
# options(hcp.input = "C:/project/PHAG_EVO_HCP_raw_only_analytical_dataset.csv",
#         hcp.output = "C:/project/raw_only_results")
# options(hcp.config = list(install_missing = TRUE))  # explicit opt-in only
# Smoke testing only (NOT publication outputs):
# options(hcp.config = list(n_pa=20L,n_pa_common=20L,n_pa_poly=10L,n_pa_q=20L,
#                          n_boot=10L,n_boot_site=10L,n_boot_cor=20L))
# Restore production defaults with options(hcp.config = NULL).
if(!identical(getOption("hcp.autorun"),FALSE)) {
  cli <- if(sys.nframe()==0L) commandArgs(trailingOnly=TRUE) else character()
  here <- if(!is.na(HCP_SCRIPT_PATH)&&nzchar(HCP_SCRIPT_PATH)) dirname(HCP_SCRIPT_PATH) else getwd()
  defaults <- c(file.path(here,"PHAG_EVO_HCP_raw_only_analytical_dataset.csv"),
                file.path(here,"PHAG_EVO_HCP_final_analytical_dataset.csv"),
                "PHAG_EVO_HCP_raw_only_analytical_dataset.csv",
                "PHAG_EVO_HCP_final_analytical_dataset.csv")
  found <- defaults[file.exists(defaults)]
  input <- if(length(cli)>=1L) cli[1] else getOption("hcp.input",if(length(found)) found[1] else defaults[1])
  output <- if(length(cli)>=2L) cli[2] else getOption("hcp.output",file.path(here,"PHAG_EVO_HCP_raw_only_outputs"))
  hcp_analysis(input,output,config=getOption("hcp.config",list()))
}
