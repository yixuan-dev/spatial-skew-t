normalize_yes_no <- function(x) {
  txt <- tolower(trimws(as.character(x)))
  txt[!nzchar(txt)] <- NA_character_
  txt
}

parse_mrts_int <- function(x) {
  raw <- trimws(as.character(x))
  raw[!nzchar(raw)] <- NA_character_
  num <- suppressWarnings(as.numeric(raw))
  out <- suppressWarnings(as.integer(round(num)))
  out[!is.finite(num)] <- NA_integer_
  out
}

build_mrts_target_pairs <- function(settings) {
  settings_local <- settings
  if (!("setting_num" %in% names(settings_local)) && "setting" %in% names(settings_local)) {
    settings_local$setting_num <- suppressWarnings(as.integer(settings_local$setting))
  }

  needed_cols <- c("setting_num", "method", "knots", "thresh", "CMAQ", "TS", "mrts")
  missing_cols <- needed_cols[!needed_cols %in% names(settings_local)]
  if (length(missing_cols) > 0) {
    stop("settings.csv is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  settings_local <- settings_local[!is.na(settings_local$setting_num), , drop = FALSE]
  settings_local$method_key <- tolower(trimws(as.character(settings_local$method)))
  settings_local$knots_num <- suppressWarnings(as.integer(settings_local$knots))
  settings_local$thresh_num <- suppressWarnings(as.numeric(settings_local$thresh))
  settings_local$CMAQ_key <- normalize_yes_no(settings_local$CMAQ)
  settings_local$TS_key <- normalize_yes_no(settings_local$TS)
  settings_local$mrts_k <- parse_mrts_int(settings_local$mrts)

  baseline_pool <- settings_local[
    settings_local$CMAQ_key == "yes" &
      settings_local$TS_key == "yes" &
      settings_local$knots_num == 1 &
      is.na(settings_local$mrts_k),
    ,
    drop = FALSE
  ]

  mrts_pool <- settings_local[
    settings_local$CMAQ_key == "yes" &
      settings_local$TS_key == "yes" &
      settings_local$knots_num == 1 &
      !is.na(settings_local$mrts_k) &
      settings_local$mrts_k > 0,
    ,
    drop = FALSE
  ]

  if (nrow(mrts_pool) == 0) {
    stop("No MRTS settings found in settings.csv (expected rows like 201:209).")
  }

  mrts_pool <- mrts_pool[order(mrts_pool$mrts_k, mrts_pool$thresh_num, mrts_pool$setting_num), , drop = FALSE]

  rows <- list()
  row_id <- 1

  for (i in seq_len(nrow(mrts_pool))) {
    proposed_row <- mrts_pool[i, , drop = FALSE]
    same_basis <- baseline_pool[
      baseline_pool$method_key == proposed_row$method_key &
        baseline_pool$knots_num == proposed_row$knots_num &
        baseline_pool$thresh_num == proposed_row$thresh_num &
        baseline_pool$CMAQ_key == proposed_row$CMAQ_key &
        baseline_pool$TS_key == proposed_row$TS_key,
      ,
      drop = FALSE
    ]

    if (nrow(same_basis) == 0) {
      next
    }

    same_basis <- same_basis[order(same_basis$setting_num), , drop = FALSE]
    baseline_row <- same_basis[1, , drop = FALSE]

    rows[[row_id]] <- data.frame(
      pair_id = row_id,
      method = proposed_row$method[1],
      knots = proposed_row$knots_num[1],
      thresh = proposed_row$thresh_num[1],
      CMAQ = proposed_row$CMAQ[1],
      TS = proposed_row$TS[1],
      mrts_k = proposed_row$mrts_k[1],
      baseline_setting = baseline_row$setting_num[1],
      mrts_setting = proposed_row$setting_num[1],
      stringsAsFactors = FALSE
    )
    row_id <- row_id + 1
  }

  if (length(rows) == 0) {
    stop("Unable to construct MRTS target pairs from settings.csv")
  }

  out <- do.call(rbind, rows)
  out <- out[order(out$mrts_k, out$thresh, out$mrts_setting), , drop = FALSE]
  rownames(out) <- NULL
  out
}

