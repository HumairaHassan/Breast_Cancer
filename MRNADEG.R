FILE <- "GSE65194_series_matrix.txt"
OUT  <- "gse65194_samples.csv"

lines <- readLines(FILE, warn = FALSE, encoding = "UTF-8")

sample_ids <- NULL
char_lines <- list()

for (line in lines) {
  if (startsWith(line, "!Sample_geo_accession")) {
    sample_ids <- strsplit(line, "\t", fixed = TRUE)[[1]][-1]
  } else if (startsWith(line, "!Sample_characteristics_ch1")) {
    char_lines[[length(char_lines) + 1]] <- strsplit(line, "\t", fixed = TRUE)[[1]][-1]
  }
}

if (is.null(sample_ids)) {
  stop("Could not find !Sample_geo_accession in the file.")
}

# Combine all characteristics into one text per sample
if (length(char_lines) == 0) {
  combined_text <- rep("", length(sample_ids))
} else {
  char_mat <- do.call(rbind, char_lines)
  combined_text <- apply(char_mat, 2, function(x) tolower(paste(x, collapse = " | ")))
}

label_sample <- function(text) {
  if (grepl("normal|control|adjacent|non-tumor|nontumor", text, ignore.case = TRUE)) {
    return("Normal")
  }
  if (grepl("tumor|tumour|cancer|carcinoma|tnbc|triple negative|triple-negative", text, ignore.case = TRUE)) {
    return("Tumor")
  }
  return("Other")
}

groups <- vapply(combined_text, label_sample, character(1))

samples <- data.frame(
  sample_id = sample_ids,
  group = groups,
  stringsAsFactors = FALSE
)

# Keep only tumor/normal
samples <- samples[samples$group %in% c("Tumor", "Normal"), ]

write.csv(samples, OUT, row.names = FALSE)

cat("Saved:", OUT, "\n")
print(table(samples$group))
print(utils::head(samples, 10))
