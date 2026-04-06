library(readxl)
library(ggplot2)

# Load KEGG data
df <- read_excel("DAVIDChartReport_Up_Genes_KEGG.xlsx")

# Clean columns
df$Count <- as.numeric(df$Count)
df$PValue <- as.numeric(df$PValue)
df$`List Total` <- as.numeric(df$`List Total`)

# Derive metrics
df$GeneRatio <- df$Count / df$`List Total`
df$logP <- -log10(df$PValue)

# Remove invalid rows and sort
df <- df[is.finite(df$GeneRatio) & is.finite(df$logP), ]
df <- df[order(df$GeneRatio), ]
df$Term <- factor(df$Term, levels = df$Term)

# Plot
p <- ggplot(df, aes(x = GeneRatio, y = Term, size = Count, color = logP)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "KEGG",
    x = "Gene Ratio",
    y = NULL,
    color = "-log10(P Value)",
    size = "Gene"
  ) +
  theme_minimal(base_size = 12)

print(p)

# Optional save
# ggsave("KEGG_plot.png", p, width = 7, height = 7, dpi = 300)
