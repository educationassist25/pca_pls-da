# YouTube Tutorial Script:
# 3D PCA and PLS-DA for Metabolomics Data in R
# Author: Your Name / BioAnalytics Lab
##############################################################

# =========================
# 1️⃣ Install Required Packages
# =========================
install.packages("plotly")        # for interactive 3D plots
install.packages("BiocManager")   # for Bioconductor packages
BiocManager::install("mixOmics")  # for PLS-DA

# =========================
# 2️⃣ Load Libraries
# =========================
library(plotly)
library(mixOmics)
library(readr)
library(dplyr)

# =========================
# 3️⃣ Load CSV Data
# Replace filename with your CSV file
# Format: Sample | Group | Metabolites...
# =========================
data <- read_csv("Example_MetaPeak Area.csv", show_col_types = FALSE)

# View the first rows
View(data)

# =========================
# 4️⃣ Separate Metadata and Features
# =========================
group <- as.factor(data$Group)        # class/group info
sample_id <- data$Sample              # sample IDs

feature_data <- data[, -(1:2)]        # remove Sample and Group columns
feature_matrix <- as.matrix(feature_data)

# =========================
# 5️⃣ Handle Missing Values
# Replace NA with half of minimum value (common in metabolomics)
# =========================
min_val <- min(feature_matrix, na.rm = TRUE)
feature_matrix[is.na(feature_matrix)] <- min_val / 2

# =========================
# 6️⃣ Log2 Transformation
# Recommended for peak area data
# =========================
log_data <- log2(feature_matrix + 1)

# =========================
# 7️⃣ Autoscale Data
# Mean-center + unit variance scaling
# =========================
scaled_data <- scale(log_data)

# =========================
# 8️⃣ PCA Analysis (3D)
# =========================
pca_res <- prcomp(scaled_data, center = TRUE, scale. = FALSE)

pca_3d <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  PC3 = pca_res$x[,3],
  Group = group,
  Sample = sample_id
)

# =========================
# 9️⃣ 3D PCA Plot
# Interactive plot using plotly
# =========================
pca_plot_3d <- plot_ly(
  pca_3d,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Group,
  text = ~Sample,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(title = "3D PCA Plot")

pca_plot_3d

# =========================
# 10️⃣ PLS-DA Analysis (3 Components)
# =========================
plsda_res <- plsda(scaled_data, group, ncomp = 3)

pls_scores <- plsda_res$variates$X

pls_3d <- data.frame(
  Comp1 = pls_scores[,1],
  Comp2 = pls_scores[,2],
  Comp3 = pls_scores[,3],
  Group = group,
  Sample = sample_id
)

# =========================
# 11️⃣ 3D PLS-DA Plot
# Interactive plot using plotly
# =========================
plsda_plot_3d <- plot_ly(
  pls_3d,
  x = ~Comp1,
  y = ~Comp2,
  z = ~Comp3,
  color = ~Group,
  text = ~Sample,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(title = "3D PLS-DA Plot")

plsda_plot_3d

# =========================
# ✅ Optional: Save Plots as HTML
# =========================
library(htmlwidgets)
saveWidget(pca_plot_3d, "PCA_3D.html")
saveWidget(plsda_plot_3d, "PLSDA_3D.html")

# =========================
# Tutorial Complete
# Now you can rotate plots interactively and screen-record for YouTube
# =========================

