#############################################
# 0. setup ##################################
#############################################
library(dplyr)
library(mpra)
library(stringr)
library(lme4)
library(lmerTest)
library(parallel)
library(argparse)
library(ggplot2)

parser <- ArgumentParser()
parser$add_argument("--input_file", required=TRUE, help="Directory containing pooled RNA/DNA/meta CSV files")
parser$add_argument("--prefix", required=TRUE, help="Prefix used for pooled input and output files")
parser$add_argument("--config_json", required=TRUE, help="Path to JSON configuration file")
parser$add_argument("--output_dir", required=TRUE, help="Output directory for MPRA results")
parser$add_argument("--num_cores", type="integer", default=1, help="Number of parallel cores")
args <- parser$parse_args()

input_file <- args$input_file
prefix <- args$prefix
output_dir <- args$output_dir
nc <- args$num_cores
config_json <- jsonlite::fromJSON(args$config_json)
thre_active <- config_json$thresholds$mpralm$logfc_thre_active
thre_allelic <- config_json$thresholds$mpralm$logfc_thre_allelic

dir.create(file.path(output_dir, "mpralm"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "MPRA_QC"), recursive = TRUE, showWarnings = FALSE)

rna_mat <- read.csv(file.path(input_file, paste0(prefix, "_RNA.csv")), row.names = 1)
dna_mat <- read.csv(file.path(input_file, paste0(prefix, "_DNA.csv")), row.names = 1)
meta_mat <- read.csv(file.path(input_file, paste0(prefix, "_meta.csv")), row.names = 1)
rna_mat[rna_mat == 0] <- NA
dna_mat[dna_mat == 0] <- NA

var_meta <- data.frame(oligo = rownames(rna_mat)) %>%
  mutate(
    Class = case_when(
      str_detect(oligo, "chr") ~ "variant",
      str_detect(oligo, "Scramble") ~ "negative_control",
      TRUE ~ "positive_control"
    ),
    Element = case_when(
      Class != "variant" ~ NA_character_,
      TRUE ~ {
        ref <- str_match(oligo, ".*:([ACGTN]+):([ACGTN]+)_[ACGTN]+$")[, 2]
        alt <- str_match(oligo, ".*:([ACGTN]+):([ACGTN]+)_[ACGTN]+$")[, 3]
        obs <- str_match(oligo, ".*:([ACGTN]+):([ACGTN]+)_([ACGTN]+)$")[, 4]

        ifelse(obs == ref, "Ref", ifelse(obs == alt, "Alt", NA_character_))
      }
    ),
    vid = if_else(
      Class == "variant",
      sub("_[^_]+$", "", oligo),
      NA_character_
    )
  )

#############################################
# 1. Allelic analysis #######################
#############################################

allelic_matrix <- function(mat, meta) {
  stopifnot(all(rownames(mat) %in% meta$oligo))
  meta2 <- meta %>%
    filter(Class == "variant", !is.na(vid), Element %in% c("Ref", "Alt")) %>%
    select(oligo, vid, Element)

  keep_vid <- meta2 %>%
    distinct(vid, Element, oligo) %>%
    dplyr::count(vid, Element) %>%
    tidyr::pivot_wider(names_from = Element, values_from = n, values_fill = 0) %>%
    filter(Ref == 1, Alt == 1) %>%
    pull(vid)

  meta2 <- meta2 %>% filter(vid %in% keep_vid)

  ref_map <- meta2 %>% filter(Element == "Ref") %>% distinct(vid, oligo)
  alt_map <- meta2 %>% filter(Element == "Alt") %>% distinct(vid, oligo)
  pair_map <- inner_join(ref_map, alt_map, by = "vid", suffix = c("_Ref", "_Alt"))

  ref_mat <- mat[pair_map$oligo_Ref, , drop = FALSE]
  alt_mat <- mat[pair_map$oligo_Alt, , drop = FALSE]

  rownames(ref_mat) <- pair_map$vid
  rownames(alt_mat) <- pair_map$vid
  colnames(ref_mat) <- paste0(colnames(ref_mat), "_Ref")
  colnames(alt_mat) <- paste0(colnames(alt_mat), "_Alt")

  cbind(ref_mat, alt_mat)
}

rna_allelic <- allelic_matrix(rna_mat, var_meta)
dna_allelic <- allelic_matrix(dna_mat, var_meta)

variant_id <- rownames(rna_allelic)
n_samples <- ncol(rna_mat)
allelic_predictor <- rep(c("Ref", "Alt"), each = n_samples)
sex_predictor <- rep(meta_mat$Sex, 2)
pnd_predictor <- rep(meta_mat$PND, 2)
rin_predictor <- rep(meta_mat$RIN, 2)
batch_predictor <- rep(meta_mat$batch, 2)
replicates <- rep(colnames(rna_mat), 2)

design_allelic <- model.matrix(~ allelic_predictor + sex_predictor + rin_predictor + batch_predictor)

colnames(rna_allelic) <- NULL
colnames(dna_allelic) <- NULL

mpra_set <- MPRASet(DNA = as.matrix(dna_allelic), RNA = as.matrix(rna_allelic), eid = variant_id)

tr <- mpralm(
  object = mpra_set,
  design = design_allelic,
  aggregate = "none",
  normalize = TRUE,
  block = replicates,
  model_type = "corr_groups"
)

mpra_variant <- topTable(tr, coef = 2, number = Inf)
se_allelic <- tr$stdev.unscaled[, 2] * sqrt(tr$s2.post)
mpra_variant$SE_allelic <- se_allelic[rownames(mpra_variant)]
mpra_variant$fdr <- p.adjust(mpra_variant$P.Value, "BH")
mpra_variant$postProbEffect <- 1 / (1 + exp(-mpra_variant$B))
mpra_variant$CI_lower_95 <- mpra_variant$logFC - 1.96 * mpra_variant$SE_allelic
mpra_variant$CI_upper_95 <- mpra_variant$logFC + 1.96 * mpra_variant$SE_allelic
mpra_variant$Sig_allelic <- ifelse(
  (mpra_variant$fdr < 0.05) & (abs(mpra_variant$logFC) > log2(thre_allelic)),
  "MPRA allelic",
  ifelse(mpra_variant$P.Value > 0.05, "MPRA non-allelic", "Allelic no decision")
)

mpra_variant <- mpra_variant %>%
  dplyr::rename(
    logFC_allelic = logFC,
    AveExpr_allelic = AveExpr,
    t_allelic = t,
    pVal_allelic = P.Value,
    FDR_allelic = adj.P.Val,
    B_allelic = B,
    fdr_mpralm = fdr
  )

write.csv(mpra_variant, file = file.path(output_dir, "mpralm", paste0(prefix, "_mpra_variant.csv")), row.names = TRUE)

get_all_coef_tables <- function(tr, design, number = Inf, adjust = "BH") {
  coef_names <- colnames(design)
  stopifnot(!is.null(coef_names))

  out_list <- lapply(seq_along(coef_names), function(j) {
    tt <- limma::topTable(tr, coef = j, number = number, sort.by = "none")
    tt <- tibble::rownames_to_column(as.data.frame(tt), var = "id")
    tt$term <- coef_names[j]
    tt
  })

  out <- bind_rows(out_list)
  out %>%
    group_by(term) %>%
    mutate(FDR = p.adjust(P.Value, method = adjust)) %>%
    ungroup()
}

all_coef_df <- get_all_coef_tables(tr, design_allelic, number = Inf)

all_coef_df %>%
  select(id, term, logFC, AveExpr, t, P.Value, adj.P.Val, FDR, B) %>%
  write.csv(file = file.path(output_dir, "mpralm", paste0(prefix, "_mpra_variant_all_coef.csv")), row.names = FALSE)

#############################################
# 2. Active analysis ########################
#############################################

negative_id <- var_meta %>% filter(Class == "negative_control") %>% pull(oligo)
logRD_mat <- log2(rna_mat / dna_mat)
neg_ctrl_median <- apply(logRD_mat[negative_id, ], 2, function(x) median(x, na.rm = TRUE))
target_list <- var_meta %>% filter(Class == "variant") %>% pull(oligo)

formula_active <- as.formula("expr ~ group + Sex + RIN + batch + (1|replicate)")

one_target <- function(target) {
  tryCatch({
    rd_vector <- t(logRD_mat[target, ]) %>% as.data.frame()

    df <- data.frame(
      expr_target = rd_vector,
      expr_control = neg_ctrl_median,
      check.names = FALSE
    )
    colnames(df) <- c("expr_target", "expr_control")
    df <- merge(df, meta_mat, by.x = "row.names", by.y = "row.names")
    df <- df[!is.na(df$expr_target) & !is.na(df$expr_control), ]

    df_long <- rbind(
      transform(df, expr = expr_target, group = "target"),
      transform(df, expr = expr_control, group = "control")
    )
    df_long$group <- factor(df_long$group, levels = c("control", "target"))
    colnames(df_long)[1] <- "replicate"

    fit <- lmer(formula_active, data = df_long)
    s <- summary(fit)$coefficients

    data.frame(
      target = target,
      beta = s["grouptarget", "Estimate"],
      p_value = s["grouptarget", "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(target = target, beta = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
  })
}

res_list <- mclapply(target_list, one_target, mc.cores = nc)
mpra_element <- bind_rows(res_list)

mpra_element$FDR <- p.adjust(mpra_element$p_value, method = "BH")
mpra_element$Sig_active <- ifelse(
  (mpra_element$FDR < 0.05) & (mpra_element$beta > log2(thre_active)),
  "Enhancing CRE",
  ifelse(
    (mpra_element$FDR < 0.05) & (mpra_element$beta < log2(1 / thre_active)),
    "Silencing CRE",
    ifelse((mpra_element$p_value > 0.05), "Non Active CRE", "nodecision")
  )
)

colnames(mpra_element) <- c("target", "logFC_active", "pVal_active", "FDR_active", "Sig_active")

mpra_element %>%
  write.csv(file = file.path(output_dir, "mpralm", paste0(prefix, "_mpra_element.csv")), row.names = FALSE)

mpra_element <- mpra_element %>%
  mutate(
    vid = sub("_[^_]+$", "", target),
    allele = sub(".*_([^_]+)$", "\\1", target),
    base = str_remove(target, "_[^_]+$"),
    A1 = str_match(base, ":([^:]+):([^:]+)$")[, 2],
    A2 = str_match(base, ":([^:]+):([^:]+)$")[, 3],
    Element = case_when(
      allele == A1 ~ "Ref",
      allele == A2 ~ "Alt",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-base, -A1, -A2, -allele)

mpra_element_arranged <- merge(
  mpra_element %>% filter(Element == "Ref"),
  mpra_element %>% filter(Element == "Alt"),
  by = "vid",
  suffixes = c("_Ref", "_Alt")
) %>%
  mutate(
    num_active = (Sig_active_Ref == "Enhancing CRE") + (Sig_active_Alt == "Enhancing CRE"),
    num_repressive = (Sig_active_Ref == "Silencing CRE") + (Sig_active_Alt == "Silencing CRE"),
    Sig_active_concensus = case_when(
      (num_active >= 1) & (num_repressive == 0) ~ "Enhancing CRE",
      (num_active == 0) & (num_repressive >= 1) ~ "Silencing CRE",
      (num_active >= 1) & (num_repressive >= 1) ~ "Other CRE",
      TRUE ~ "Non Active CRE"
    )
  ) %>%
  dplyr::select(
    vid, logFC_active_Ref, pVal_active_Ref, FDR_active_Ref, Sig_active_Ref,
    logFC_active_Alt, pVal_active_Alt, FDR_active_Alt, Sig_active_Alt,
    Sig_active_concensus
  )

p <- mpra_element %>%
  ggplot(aes(x = logFC_active, y = -log(FDR_active))) +
  geom_point() +
  theme_bw() +
  labs(x = "logFC_active", y = "-log(FDR_active)")

ggsave(
  filename = file.path(output_dir, "MPRA_QC", paste0(prefix, "_mpra_element_volcano_plot.png")),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

mpra_element %>%
  write.csv(file = file.path(output_dir, "mpralm", paste0(prefix, "_mpra_element.csv")), row.names = FALSE)

mpra_element_arranged %>%
  write.csv(file = file.path(output_dir, "mpralm", paste0(prefix, "_mpra_element_concensus.csv")), row.names = FALSE)

#############################################
# 4. EmVar analysis #########################
#############################################

mpra_element_arranged <- read.csv(file = file.path(output_dir, "mpralm", paste0(prefix, "_mpra_element_concensus.csv")))
mpra_variant <- read.csv(file = file.path(output_dir, "mpralm", paste0(prefix, "_mpra_variant.csv")))

if ("X" %in% colnames(mpra_variant)) {
  mpra_variant$vid <- mpra_variant$X
} else if ("eid" %in% colnames(mpra_variant)) {
  mpra_variant$vid <- mpra_variant$eid
} else {
  stop("mpra_variant must contain either X or eid column")
}

EmVar_summary <- mpra_variant %>%
  inner_join(mpra_element_arranged, by = "vid") %>%
  mutate(EmVar = case_when(
    Sig_allelic == "MPRA allelic" & Sig_active_concensus == "Enhancing CRE" ~ "Enhancing emVar",
    Sig_allelic == "MPRA allelic" & Sig_active_concensus == "Silencing CRE" ~ "Silencing emVar",
    Sig_allelic == "MPRA allelic" & Sig_active_concensus %in% c("Other CRE", "Other") ~ "Other emVar",
    Sig_allelic == "MPRA non-allelic" & Sig_active_concensus == "Non Active CRE" ~ "Non-emVar",
    TRUE ~ "No decision"
  ))

EmVar_summary <- EmVar_summary %>%
  dplyr::rename(eid = vid) %>%
  dplyr::select(
    eid, logFC_allelic, pVal_allelic, FDR_allelic, Sig_allelic, logFC_active_Ref,
    pVal_active_Ref, FDR_active_Ref, Sig_active_Ref, logFC_active_Alt, pVal_active_Alt,
    FDR_active_Alt, Sig_active_Alt, Sig_active_concensus, EmVar
  )

EmVar_summary %>%
  write.csv(file = file.path(output_dir, "mpralm", paste0(prefix, "_EmVar_summary.csv")), row.names = FALSE)

table_summary <- function() {
  print("Active Variants")
  print(EmVar_summary$Sig_active_concensus %>% table())
  print("Allelic Variants")
  print(EmVar_summary$Sig_allelic %>% table())
  print("EmVars")
  print(EmVar_summary$EmVar %>% table())
}

table_summary()
