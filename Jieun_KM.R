#!/usr/bin/env Rscript
# Kaplan–Meier 생존분석 워크플로우 (R) - 첨부 NSCLC 코호트용
#
# 이 스크립트가 하는 일(재현 가능한 end-to-end):
# 1) 첨부 TSV 로드 (844명; 임상 + 치료 + WES 파생 특성)
# 2) time-to-event 분석용 데이터셋 구성:
#    - index_date, event_date, censor_date, last_followup 및 옵션 data_cutoff 정의
#    - CNSR 플래그에서 event 지표 생성 (event = 1 - CNSR)
# 3) QC 체크:
#    - 환자 ID 유일성
#    - 결측률 프로파일링
#    - 검열 플래그 유효성(0/1)
#    - 추적 일관성(OS >= PFS), 음수 시간 여부
#    - 날짜 일관성(end date >= index date, event/censor 날짜가 event 지표와 일치)
# 4) 임상/치료/유전체 하위군으로 층화한 KM 분석 수행
# 5) log-rank test로 생존곡선 비교
# 6) 해석/가설생성에 적합한 표/그림 산출물 자동 생성
# 7) 옵션: ~400 유전자 변이 log-rank 스캔 + BH-FDR + 상위 유전자 KM 플롯 생성
#
# 중요: 날짜 관련 주의
# - 첨부 TSV는 OS/PFS '기간(duration)'만 제공하고 실제 달력 날짜(index/event/last follow-up)가 없습니다.
# - 날짜 기반 워크플로우 데모 및 날짜 QC를 위해, 제공된 duration과 검열 플래그에 맞춰
#   'pseudo 날짜'를 결정적으로 생성합니다.
# - KM 추정 자체는 제공된 duration(OS/PFS)을 그대로 사용합니다.
#
# 실행:
#   Rscript km_workflow_R_ko.R ML-ready-oak-poplar.tsv km_outputs
#
# 옵션(환경변수):
#   DATA_CUTOFF=2016-01-01   (행정 검열 cutoff; 데모)
#   INDEX_DATE=2013-01-01    (pseudo index date; 데모)
#   MIN_MUT=10               (유전자 스캔 최소 변이 환자 수)
#   TOP_GENES=10             (BH q 기준 상위 유전자 플롯 개수)
#   BTMB_CUTPOINT=16         (bTMB 고정 컷포인트; 미지정 시 median split)
#

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("사용법: Rscript Jieun_KM.R <input_tsv> <outdir>\n")
  quit(status = 1)
}

input_tsv <- args[1]
outdir <- args[2]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 패키지 설치/로드
# -------------------------
pkgs <- c("readr", "dplyr", "tidyr", "stringr", "survival", "survminer", "ggplot2", "purrr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(survival)
library(survminer)
library(ggplot2)
library(purrr)

# -------------------------
# 헬퍼
# -------------------------

bh_fdr <- function(p) {
  n <- length(p)
  o <- order(p)
  p_sorted <- p[o]
  q <- rep(NA_real_, n)
  prev <- 1.0
  for (i in n:1) {
    rank <- i
    val <- p_sorted[i] * n / rank
    prev <- min(prev, val)
    q[i] <- prev
  }
  out <- rep(NA_real_, n)
  out[o] <- pmin(pmax(q, 0), 1)
  out
}

safe_name <- function(x) {
  x %>% str_replace_all("[/ ]+", "_")
}

# -------------------------
# 데이터 로드
# -------------------------
df <- read_tsv(input_tsv, show_col_types = FALSE)

mut_cols <- names(df)[str_detect(names(df), "^mutations_")]
vaf_cols <- names(df)[str_detect(names(df), "^max_vaf_")]
base_cols <- setdiff(names(df), c(mut_cols, vaf_cols))

# -------------------------
# time-to-event + 변이 0/1
# -------------------------
df <- df %>%
  mutate(
    OS_event  = ifelse(`OS.CNSR` == 0, 1, 0),
    PFS_event = ifelse(`PFS.CNSR` == 0, 1, 0)
  )

if (length(mut_cols) > 0) {
  df <- df %>% mutate(across(all_of(mut_cols), ~ ifelse(is.na(.x), 0L, 1L)))
}

# -------------------------
# 파생 그룹 변수 생성(임상/WES/QC)
# -------------------------
btmb_cutpoint <- Sys.getenv("BTMB_CUTPOINT", unset = NA_character_)
btmb_cutpoint <- suppressWarnings(as.numeric(btmb_cutpoint))

if ("METSITES" %in% names(df)) {
  df <- df %>%
    mutate(extent_group = cut(
      METSITES,
      breaks = c(-Inf, 0, 2, Inf),
      labels = c("0 sites", "1-2 sites", "3+ sites")
    ))
}

if ("btmb" %in% names(df)) {
  if (is.na(btmb_cutpoint)) {
    med <- median(df$btmb, na.rm = TRUE)
    df <- df %>% mutate(btmb_group = ifelse(btmb >= med, paste0("High (>=median=", round(med, 1), ")"), paste0("Low (<median=", round(med, 1), ")")))
  } else {
    df <- df %>% mutate(btmb_group = ifelse(btmb >= btmb_cutpoint, paste0("High (>= ", btmb_cutpoint, ")"), paste0("Low (< ", btmb_cutpoint, ")")))
  }
}

if ("Median_exon_coverage" %in% names(df)) {
  med <- median(df$Median_exon_coverage, na.rm = TRUE)
  df <- df %>% mutate(coverage_group = ifelse(Median_exon_coverage >= med, paste0("High (>=median=", round(med, 0), ")"), paste0("Low (<median=", round(med, 0), ")")))
}

if ("MSAF" %in% names(df)) {
  med <- median(df$MSAF, na.rm = TRUE)
  df <- df %>% mutate(msaf_group = ifelse(MSAF >= med, paste0("High (>=median=", sprintf("%.3f", med), ")"), paste0("Low (<median=", sprintf("%.3f", med), ")")))
}

if ("cfDNA_Input_ng" %in% names(df)) {
  med <- median(df$cfDNA_Input_ng, na.rm = TRUE)
  df <- df %>% mutate(cfdna_input_group = ifelse(cfDNA_Input_ng >= med, paste0("High (>=median=", sprintf("%.2f", med), ")"), paste0("Low (<median=", sprintf("%.2f", med), ")")))
}

if ("TC3IC3" %in% names(df)) {
  df <- df %>%
    mutate(
      TC3IC3_clean = case_when(
        is.na(TC3IC3) ~ "Missing",
        TC3IC3 == "UNKNOWN" ~ "Missing",
        TRUE ~ TC3IC3
      ),
      pdl1_tc3ic3_binary = case_when(
        TC3IC3_clean == "TC3 or IC3" ~ "High (TC3/IC3)",
        TC3IC3_clean == "TC0/1/2 and IC0/1/2" ~ "Low (TC0-2/IC0-2)",
        TRUE ~ "Missing/Other"
      )
    ) %>%
    select(-TC3IC3_clean)
}

# 예시 유전자 하위군 플래그 (컬럼이 있으면 생성)
for (g in c("TP53", "KEAP1", "STK11", "EGFR")) {
  col <- paste0("mutations_", g)
  if (col %in% names(df)) {
    df[[paste0(g, "_mut")]] <- df[[col]]
  }
}

# -------------------------
# pseudo 날짜 생성(데모용; OS/PFS duration과 일치)
# -------------------------
index_date <- Sys.getenv("INDEX_DATE", unset = "2013-01-01")
index_date <- as.Date(index_date)

data_cutoff <- Sys.getenv("DATA_CUTOFF", unset = NA_character_)
if (!is.na(data_cutoff)) {
  data_cutoff <- as.Date(data_cutoff)
}

days_per_month <- 365.25 / 12.0

df <- df %>%
  mutate(index_date = index_date)

df <- df %>%
  mutate(
    OS_end_date  = index_date + as.difftime(OS * days_per_month, units = "days"),
    PFS_end_date = index_date + as.difftime(PFS * days_per_month, units = "days")
  )

df <- df %>%
  mutate(
    OS_event_date   = ifelse(OS_event == 1, OS_end_date, as.Date(NA)),
    OS_censor_date  = ifelse(OS_event == 0, OS_end_date, as.Date(NA)),
    PFS_event_date  = ifelse(PFS_event == 1, PFS_end_date, as.Date(NA)),
    PFS_censor_date = ifelse(PFS_event == 0, PFS_end_date, as.Date(NA)),
    OS_last_followup_date  = OS_end_date,
    PFS_last_followup_date = PFS_end_date
  )

if (!is.na(data_cutoff)) {
  over_os <- df$OS_end_date > data_cutoff
  df$OS_end_date[over_os] <- data_cutoff
  df$OS_event[over_os] <- 0
  df$OS_event_date[over_os] <- as.Date(NA)
  df$OS_censor_date[over_os] <- data_cutoff

  over_pfs <- df$PFS_end_date > data_cutoff
  df$PFS_end_date[over_pfs] <- data_cutoff
  df$PFS_event[over_pfs] <- 0
  df$PFS_event_date[over_pfs] <- as.Date(NA)
  df$PFS_censor_date[over_pfs] <- data_cutoff

  df$OS  <- as.numeric(difftime(df$OS_end_date, df$index_date, units = "days")) / days_per_month
  df$PFS <- as.numeric(difftime(df$PFS_end_date, df$index_date, units = "days")) / days_per_month
}

# -------------------------
# QC
# -------------------------
qc <- list()

if ("PtID" %in% names(df)) {
  qc$n_patients <- n_distinct(df$PtID)
  qc$n_rows <- nrow(df)
  qc$ptid_duplicates <- nrow(df) - n_distinct(df$PtID)
} else {
  qc$n_patients <- nrow(df)
  qc$n_rows <- nrow(df)
  qc$ptid_duplicates <- NA
}

base_missing <- df %>%
  summarize(across(all_of(base_cols), ~ mean(is.na(.x)))) %>%
  pivot_longer(cols = everything(), names_to = "col", values_to = "missing_rate") %>%
  arrange(desc(missing_rate))

qc$top_missing_base_cols <- head(base_missing, 15)

qc$OS_CNSR_unique_values  <- sort(unique(df$`OS.CNSR`))
qc$PFS_CNSR_unique_values <- sort(unique(df$`PFS.CNSR`))
qc$OS_CNSR_invalid_values  <- sum(!(df$`OS.CNSR` %in% c(0, 1)))
qc$PFS_CNSR_invalid_values <- sum(!(df$`PFS.CNSR` %in% c(0, 1)))

qc$OS_missing <- sum(is.na(df$OS))
qc$PFS_missing <- sum(is.na(df$PFS))
qc$OS_negative <- sum(df$OS < 0, na.rm = TRUE)
qc$PFS_negative <- sum(df$PFS < 0, na.rm = TRUE)
qc$OS_lt_PFS_count <- sum(df$OS < df$PFS, na.rm = TRUE)

qc$OS_end_before_index <- sum(df$OS_end_date < df$index_date, na.rm = TRUE)
qc$PFS_end_before_index <- sum(df$PFS_end_date < df$index_date, na.rm = TRUE)

qc$OS_event_censor_date_mismatch <- sum((df$OS_event == 1 & is.na(df$OS_event_date)) | (df$OS_event == 0 & is.na(df$OS_censor_date)))
qc$PFS_event_censor_date_mismatch <- sum((df$PFS_event == 1 & is.na(df$PFS_event_date)) | (df$PFS_event == 0 & is.na(df$PFS_censor_date)))

qc$n_mut_cols <- length(mut_cols)
qc$n_vaf_cols <- length(vaf_cols)
qc$data_cutoff <- ifelse(is.na(data_cutoff), as.character(max(df$OS_end_date, df$PFS_end_date, na.rm = TRUE)), as.character(data_cutoff))

qc_path <- file.path(outdir, "QC_report.txt")
sink(qc_path)
cat("QC 리포트 (R)\n")
cat("====================\n\n")
cat("환자/행:\n")
cat("  n_patients:", qc$n_patients, "\n")
cat("  n_rows:", qc$n_rows, "\n")
cat("  ptid_duplicates:", qc$ptid_duplicates, "\n\n")

cat("검열 플래그 값:\n")
cat("  OS.CNSR unique:", paste(qc$OS_CNSR_unique_values, collapse = ", "), "\n")
cat("  PFS.CNSR unique:", paste(qc$PFS_CNSR_unique_values, collapse = ", "), "\n")
cat("  OS.CNSR invalid:", qc$OS_CNSR_invalid_values, "\n")
cat("  PFS.CNSR invalid:", qc$PFS_CNSR_invalid_values, "\n\n")

cat("기간(duration) 체크:\n")
cat("  OS missing:", qc$OS_missing, "\n")
cat("  PFS missing:", qc$PFS_missing, "\n")
cat("  OS negative:", qc$OS_negative, "\n")
cat("  PFS negative:", qc$PFS_negative, "\n")
cat("  OS < PFS count:", qc$OS_lt_PFS_count, "\n\n")

cat("날짜(pseudo) 체크:\n")
cat("  OS_end_before_index:", qc$OS_end_before_index, "\n")
cat("  PFS_end_before_index:", qc$PFS_end_before_index, "\n")
cat("  OS event/censor mismatch:", qc$OS_event_censor_date_mismatch, "\n")
cat("  PFS event/censor mismatch:", qc$PFS_event_censor_date_mismatch, "\n\n")

cat("특성 컬럼 개수:\n")
cat("  n_mut_cols:", qc$n_mut_cols, "\n")
cat("  n_vaf_cols:", qc$n_vaf_cols, "\n")
cat("  data_cutoff:", qc$data_cutoff, "\n\n")

cat("결측률 상위(base cols):\n")
print(qc$top_missing_base_cols)
sink()

write_tsv(df, file.path(outdir, "analysis_ready_with_dates_and_mut0_1.tsv"))

# -------------------------
# KM suite (표 + 플롯)
# -------------------------

km_by_group <- function(data, time_col, event_col, group_col, timepoints) { # new version

  # 1) filter rows with complete essentials (OK to use .data here)
  d <- data %>%
    dplyr::filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[event_col]]),
      !is.na(.data[[group_col]])
    )

  # 2) keep groups with >=2 rows
  counts <- table(d[[group_col]])
  keep <- names(counts[counts >= 2])
  d <- d %>% dplyr::filter(.data[[group_col]] %in% keep)

  # If nothing left, fail early
  if (nrow(d) == 0L) {
    return(list(
      table = tibble::tibble(),
      logrank_stat = NA_real_,
      logrank_p = NA_real_
    ))
  }

  # 3) standardize columns for survival functions (NO .data in formulas)
  d <- d %>%
    dplyr::mutate(
      .time  = .data[[time_col]],
      .event = as.integer(.data[[event_col]]),   # ensure 0/1 integer
      .group = as.factor(.data[[group_col]])
    )

  # If only one level after filtering, log-rank is undefined
  do_logrank <- nlevels(d$.group) >= 2

  # 4) log-rank
  stat <- NA_real_
  p <- NA_real_
  if (do_logrank) {
    lr <- survival::survdiff(survival::Surv(.time, .event) ~ .group, data = d)
    stat <- unname(lr$chisq)
    df_lr <- length(lr$n) - 1
    p <- 1 - stats::pchisq(stat, df = df_lr)
  }

  # 5) KM fit
  fit <- survival::survfit(survival::Surv(.time, .event) ~ .group, data = d)

  # 6) median survival per group
  med <- survminer::surv_median(fit) %>%
    dplyr::rename(group = strata)

  # 7) survival at requested timepoints
  surv_tp <- purrr::map_dfr(timepoints, function(t) {
    s <- summary(fit, times = t)
    tibble::tibble(
      group = s$strata,      # will look like ".group=SOC"
      time = t,
      surv = s$surv
    )
  }) %>%
    dplyr::mutate(group = sub("^\\.group=", "", group))  # strip prefix

  surv_wide <- surv_tp %>%
    dplyr::mutate(time = paste0("S(", time, ")")) %>%
    tidyr::pivot_wider(names_from = time, values_from = surv)

  # 8) n and events per group (stable; no NSE tricks)
  counts2 <- d %>%
    dplyr::count(.group, name = "n") %>%
    dplyr::mutate(group = sub("^\\.group=", "", paste0(".group=", .group))) %>%
    dplyr::select(group, n)

  events2 <- d %>%
    dplyr::group_by(.group) %>%
    dplyr::summarise(events = sum(.event, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(group = sub("^\\.group=", "", paste0(".group=", .group))) %>%
    dplyr::select(group, events)

  # But med$group is already like ".group=SOC" (from surv_median), so match that:
  counts2 <- d %>%
    dplyr::count(.group, name = "n") %>%
    dplyr::mutate(group = paste0(".group=", .group)) %>%
    dplyr::select(group, n)

  events2 <- d %>%
    dplyr::group_by(.group) %>%
    dplyr::summarise(events = sum(.event, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(group = paste0(".group=", .group)) %>%
    dplyr::select(group, events)

  out <- med %>%
    dplyr::left_join(counts2, by = "group") %>%
    dplyr::left_join(events2, by = "group") %>%
    dplyr::left_join(
      surv_wide %>% dplyr::mutate(group = paste0(".group=", group)),
      by = "group"
    ) %>%
    dplyr::arrange(dplyr::desc(n))

  list(table = out, logrank_stat = stat, logrank_p = p)
} # end of new km_by group
plot_km <- function(data, time_col, event_col, group_col, title, out_png,
                    conf.int = FALSE, show_p = TRUE, show_risktable = TRUE) {

  # Minimal dependencies
  stopifnot(is.data.frame(data))
  for (nm in c(time_col, event_col, group_col)) {
    if (!nm %in% names(data)) stop("Missing column: ", nm)
  }

  # Build a small working frame with standard column names
  d <- data.frame(
    .time  = data[[time_col]],
    .event = data[[event_col]],
    .group = data[[group_col]]
  )

  # Drop missing essentials
  d <- d[!is.na(d$.time) & !is.na(d$.event) & !is.na(d$.group), , drop = FALSE]

  # Ensure event is 0/1 integer (Surv tolerates logical too, but make it clean)
  if (is.logical(d$.event)) d$.event <- as.integer(d$.event)
  d$.event <- as.integer(d$.event)

  # Ensure group is a factor (pretty legend ordering, stable behavior)
  d$.group <- as.factor(d$.group)

  # Guardrails: need at least 2 groups to plot by group meaningfully
  if (nlevels(d$.group) < 1) stop("No non-missing groups to plot.")
  if (nlevels(d$.group) == 1) {
    # still can plot a single KM curve
    fit <- survival::survfit(survival::Surv(.time, .event) ~ 1, data = d)
  } else {
    fit <- survival::survfit(survival::Surv(.time, .event) ~ .group, data = d)
  }

  # p-value only makes sense with >=2 groups
  pval_flag <- isTRUE(show_p) && nlevels(d$.group) >= 2

  # Plot with survminer
  g <- survminer::ggsurvplot(
    fit,
    data = d,
    risk.table = isTRUE(show_risktable),
    pval = pval_flag,
    conf.int = isTRUE(conf.int),
    xlab = "Time (dataset-provided duration)",
    ylab = "Survival probability",
    title = title
  )

  # Save: include risk table if present
  # ggsurvplot returns a list; use arrange_ggsurvplots to combine safely
  if (isTRUE(show_risktable)) {
    ggpubr::ggexport(g, filename = out_png, width = 7, height = 5, res = 200)
  } else {
    ggplot2::ggsave(out_png, g$plot, width = 7, height = 5, dpi = 200)
  }

  invisible(list(fit = fit, plot = g))
}
## end of plot_km

run_km_suite <- function(data, endpoint, time_col, event_col, timepoints, strata, out_subdir) {
  out_path <- file.path(outdir, out_subdir)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  master <- tibble()
  for (s in strata) {
    label <- s[[1]]
    col <- s[[2]]
    if (!(col %in% names(data))) next

    safe <- safe_name(col)
    tsv <- file.path(out_path, paste0("KM_", endpoint, "_by_", safe, ".tsv"))
    png <- file.path(out_path, paste0("KM_", endpoint, "_by_", safe, ".png"))

    res <- tryCatch({
      km_by_group(data, time_col, event_col, col, timepoints)
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(res)) {
      master <- bind_rows(master, tibble(endpoint = endpoint, stratifier = label, group_col = col, error = "KM failed"))
      next
    }

    write_tsv(res$table, tsv)
    plot_km(data, time_col, event_col, col, paste0(endpoint, " by ", label), png)

    master <- bind_rows(master, tibble(
      endpoint = endpoint,
      stratifier = label,
      group_col = col,
      logrank_stat = res$logrank_stat,
      logrank_p = res$logrank_p,
      table_tsv = basename(tsv),
      plot_png = basename(png)
    ))
  }

  write_tsv(master, file.path(out_path, paste0("KM_", endpoint, "_MASTER.tsv")))
}

strata <- list(
  list("Treatment", "Treatment"),
  list("Study cohort", "trial"),
  list("Histology", "HIST"),
  list("ECOG", "ECOGGR"),
  list("Race", "race2"),
  list("Sex", "SEX"),
  list("Smoking history", "TOBHX"),
  list("Initial extent (proxy: METSITES)", "extent_group"),
  list("bTMB group", "btmb_group"),
  list("Exon coverage group", "coverage_group"),
  list("MSAF group", "msaf_group"),
  list("cfDNA input group", "cfdna_input_group"),
  list("PD-L1 (TC3/IC3) binary", "pdl1_tc3ic3_binary"),
  list("TP53 mutation", "TP53_mut"),
  list("KEAP1 mutation", "KEAP1_mut")
)

run_km_suite(df, "OS", "OS", "OS_event", c(6, 12, 18, 24), strata, "KM_OS")
run_km_suite(df, "PFS", "PFS", "PFS_event", c(3, 6, 12), strata, "KM_PFS")

# -------------------------
# 유전자 log-rank 스캔 + BH-FDR
# -------------------------

min_mut <- Sys.getenv("MIN_MUT", unset = "10")
min_mut <- as.integer(min_mut)
top_genes <- Sys.getenv("TOP_GENES", unset = "10")
top_genes <- as.integer(top_genes)

gene_outdir <- file.path(outdir, "GeneScan_OS")
dir.create(gene_outdir, recursive = TRUE, showWarnings = FALSE)

scan_rows <- list()

if (length(mut_cols) > 0) {
  for (c in mut_cols) {
    gene <- str_replace(c, "^mutations_", "")
    x <- df %>% select(OS, OS_event, all_of(c)) %>% drop_na()

    n_mut <- sum(x[[c]] == 1)
    n_wt  <- sum(x[[c]] == 0)

    if (n_mut < min_mut || n_wt < 2) next
    x$.temp <- x[[c]]

    lr <- survdiff(Surv(OS, OS_event) ~ .temp, data = x)
    stat <- unname(lr$chisq)
    p <- 1 - pchisq(stat, df = 1)

    scan_rows[[length(scan_rows) + 1]] <- tibble(
      gene = gene,
      mut_col = c,
      n_mut = n_mut,
      n_wt = n_wt,
      logrank_stat = stat,
      p = p
    )
  }
}

scan <- bind_rows(scan_rows)
scan_tsv <- file.path(gene_outdir, "gene_logrank_scan_OS.tsv")
cat("scan")
cat(nrow(scan))

if (nrow(scan) == 0) {
  write_tsv(scan, scan_tsv)
} else {
  scan <- scan %>%
    mutate(q_BH = bh_fdr(p)) %>%
    arrange(q_BH, p, gene)

  write_tsv(scan, scan_tsv)

  top <- head(scan, top_genes)
  for (i in seq_len(nrow(top))) {
    gene <- top$gene[i]
    c <- top$mut_col[i]

    d <- df %>% select(OS, OS_event, all_of(c)) %>% drop_na() %>%
      mutate(mut_status = ifelse(.data[[c]] == 1, "Mutated", "Wild-type"))

    fit <- survfit(Surv(OS, OS_event) ~ mut_status, data = d)

    p_plot <- ggsurvplot(
      fit, data = d,
      risk.table = TRUE,
      pval = TRUE,
      conf.int = FALSE,
      xlab = "Time (months; dataset-provided duration)",
      ylab = "Survival probability",
      title = paste0(gene, ": OS (Mutated vs Wild-type)")
    )

    out_png <- file.path(gene_outdir, paste0("KM_OS_", gene, "_mut_vs_wt.png"))
    ggsave(out_png, p_plot$plot, width = 7, height = 5, dpi = 200)

    res <- km_by_group(d %>% mutate(mut_status = factor(mut_status)), "OS", "OS_event", "mut_status", c(6, 12, 18, 24))
    out_tsv <- file.path(gene_outdir, paste0("KM_OS_", gene, "_mut_vs_wt.tsv"))
    write_tsv(res$table, out_tsv)
  }
}

cat("\n완료.\n")
cat("Outputs written to:", normalizePath(outdir), "\n")
cat("Key folders:\n")
cat("  -", file.path(outdir, "KM_OS"), "\n")
cat("  -", file.path(outdir, "KM_PFS"), "\n")
cat("  -", file.path(outdir, "GeneScan_OS"), "\n")
