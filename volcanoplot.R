library(readxl)
library(ggplot2)

# =========================
# INPUTS (edit if needed)
# =========================
DEG_FILE <- "GSE45827_regulation.xlsx"
OUT_PNG <- "volcano.png"
P_COL <- "adj.P.Val"
FC_COL <- "logFC"
P_CUTOFF <- 0.05
FC_CUTOFF <- 1

# =========================
# LOAD + CLEAN
# =========================
df <- read_excel(DEG_FILE)

df[[P_COL]] <- suppressWarnings(as.numeric(df[[P_COL]]))
df[[FC_COL]] <- suppressWarnings(as.numeric(df[[FC_COL]]))

df <- df[!is.na(df[[P_COL]]) & !is.na(df[[FC_COL]]) & df[[P_COL]] > 0, ]
df$minus_log10_p <- -log10(df[[P_COL]])

# =========================
# CLASSIFY POINTS
# =========================
df$class <- "Notsig"
df$class[df[[P_COL]] < P_CUTOFF & df[[FC_COL]] >= FC_CUTOFF] <- "Up"
df$class[df[[P_COL]] < P_CUTOFF & df[[FC_COL]] <= -FC_CUTOFF] <- "Down"
df$class <- factor(df$class, levels = c("Notsig", "Up", "Down"))

# =========================
# PLOT
# =========================
p <- ggplot(df, aes_string(x = FC_COL, y = "minus_log10_p", color = "class", shape = "class")) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed") +
  geom_hline(yintercept = -log10(P_CUTOFF), linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    x = FC_COL,
    y = paste0("-log10(", P_COL, ")")
  ) +
  theme_minimal(base_size = 12)

ggsave(OUT_PNG, p, width = 10, height = 7, dpi = 300)
print(p)

cat("Saved:", OUT_PNG, "\n")
print(table(df$class))
