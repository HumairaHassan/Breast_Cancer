library(readxl)
library(ggplot2)

# --------- INPUTS ---------
infile <- "miRNA_DE (1).xlsx"   # change path if needed
sheet <- 1                      # first sheet in R is 1
out_png <- "miRNA_volcano.png"  # output image name
p_cutoff <- 0.05                # used only if Regulation column missing
logfc_cutoff <- 1.0             # used only if Regulation column missing
label_top_n <- 10               # reserved for future labeling
# --------------------------

df <- read_excel(infile, sheet = sheet)
colnames(df) <- trimws(colnames(df))

pick_col <- function(df, cands) {
  hit <- cands[cands %in% colnames(df)]
  if (length(hit) > 0) hit[1] else NA_character_
}

logfc_col <- pick_col(df, c("logFC", "log2FC", "log2FoldChange", "log_fold_change", "logFC."))
p_col <- pick_col(df, c("adj.P.Val", "adj_p_val", "FDR", "padj", "adj.P.Val.", "P.Value", "pvalue", "PValue"))
name_col <- pick_col(df, c("miRNA_ID", "miRNA", "mirna", "Gene", "Gene.symbol", "Symbol", "ID"))
reg_col <- pick_col(df, c("Regulation", "regulation", "Status", "Group"))

if (is.na(logfc_col) || is.na(p_col)) {
  stop(sprintf("Could not find required columns. Found columns: %s",
               paste(colnames(df), collapse = ", ")))
}

# Clean numeric columns
df[[logfc_col]] <- suppressWarnings(as.numeric(df[[logfc_col]]))
df[[p_col]] <- suppressWarnings(as.numeric(df[[p_col]]))

# Drop invalid rows
df <- df[!is.na(df[[logfc_col]]) & !is.na(df[[p_col]]) & df[[p_col]] > 0, ]
df$neglog10 <- -log10(df[[p_col]])

# Build regulation groups
if (!is.na(reg_col)) {
  rg <- trimws(as.character(df[[reg_col]]))
  rg[is.na(rg) | rg %in% c("nan", "NaN", "None", "")] <- "NotSig"
  rg[!rg %in% c("Up", "Down", "NotSig")] <- "NotSig"
  df$RegGroup <- rg
} else {
  df$RegGroup <- "NotSig"
  df$RegGroup[df[[p_col]] < p_cutoff & df[[logfc_col]] >= logfc_cutoff] <- "Up"
  df$RegGroup[df[[p_col]] < p_cutoff & df[[logfc_col]] <= -logfc_cutoff] <- "Down"
}

df$RegGroup <- factor(df$RegGroup, levels = c("NotSig", "Down", "Up"))

p <- ggplot(df, aes_string(x = logfc_col, y = "neglog10", shape = "RegGroup", color = "RegGroup")) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = -log10(p_cutoff)) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff)) +
  labs(
    title = "miRNA Volcano Plot",
    x = "log2 Fold Change (logFC)",
    y = "-log10(adjusted p-value)"
  ) +
  theme_minimal(base_size = 12)

ggsave(out_png, p, width = 8.5, height = 6.5, dpi = 300)
print(p)

cat("Saved:", out_png, "| rows plotted:", nrow(df), "\n")
print(table(df$RegGroup))
