# Set the directory where you want to save plots
setwd("/Users/deniztasgin/Downloads/R_clustring/R/plots")

# ========================================
# Save UMAP Plots
# ========================================

# 1. UMAP colored by condition
ggsave(
  filename = "UMAP_by_condition.png",
  plot = ggplot(normdata_umap, aes(x = UMAP_1, y = UMAP_2, colour = condition)) + 
    geom_point(size = 2) +
    theme_minimal() +
    ggtitle("UMAP by Condition"),
  width = 10, height = 8, dpi = 300
)

# 2. UMAP colored by clusters
ggsave(
  filename = "UMAP_by_clusters.png",
  plot = ggplot(normdata_umap, aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster_results$membership))) + 
    geom_point(size = 2) + 
    guides(colour = guide_legend(title = "Clusters")) +
    theme_minimal() +
    ggtitle("UMAP by Phenograph Clusters"),
  width = 10, height = 8, dpi = 300
)

# 3. UMAP with annotated clusters (with labels)
ggsave(
  filename = "UMAP_annotated_clusters.png",
  plot = DimPlot(seu, reduction = "umap", group.by = "Annotated_Cluster", 
                 pt.size = 1.5, label = TRUE, label.size = 5) + 
    theme_gray() + 
    ggtitle("UMAP with Annotated Clusters"),
  width = 10, height = 8, dpi = 300
)

# 4. UMAP split by condition (4 panels)
ggsave(
  filename = "UMAP_split_by_condition.png",
  plot = DimPlot(seu, reduction = "umap", group.by = "Annotated_Cluster", 
                 split.by = "condition", pt.size = 1.5, label = FALSE) + 
    theme_gray() + 
    ggtitle("UMAP with Annotated Clusters split by condition"),
  width = 20, height = 5, dpi = 300
)

# ========================================
# Save Dot Plots
# ========================================

# 5. Dotplot by cluster
ggsave(
  filename = "DotPlot_by_cluster.png",
  plot = Seurat::DotPlot(seu, 
                         features = c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act'), 
                         group.by = "pheno_cluster",
                         cols = c("blue","orange"), 
                         assay = "Pixel") +
    RotatedAxis() + 
    coord_flip() + 
    theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=8)),
  width = 8, height = 6, dpi = 300
)

# 6. Dotplot by annotated cluster
ggsave(
  filename = "DotPlot_by_annotated_cluster.png",
  plot = Seurat::DotPlot(seu, 
                         features = c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act'), 
                         group.by = "Annotated_Cluster",
                         cols = c("blue","orange"), 
                         assay = "Pixel") +
    RotatedAxis() + 
    coord_flip() + 
    theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=8)) +
    ggtitle("Metabolic Features by Annotated Cluster"),
  width = 8, height = 6, dpi = 300
)

# 7. Dotplot by condition
ggsave(
  filename = "DotPlot_by_condition.png",
  plot = Seurat::DotPlot(seu, 
                         features = c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act'), 
                         group.by = "condition",
                         assay = "Pixel") +
    RotatedAxis() + 
    coord_flip() + 
    theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=8)) +
    ggtitle("Metabolic Features by Condition"),
  width = 8, height = 6, dpi = 300
)

# ========================================
# Save Coxcomb Charts
# ========================================

# 8. Coxcomb dotplot (faceted by feature)
ggsave(
  filename = "Coxcomb_by_feature.png",
  plot = p_coxcomb_dotplot,
  width = 12, height = 10, dpi = 300
)

# 9. Nightingale rose chart (faceted by cluster)
ggsave(
  filename = "Nightingale_rose_by_cluster.png",
  plot = p_radial_coxcomb,
  width = 12, height = 10, dpi = 300
)

# ========================================
# Save Feature Plots
# ========================================

# 10. Feature plots (all 7 metabolic features)
ggsave(
  filename = "FeaturePlot_all_metabolic_features.png",
  plot = FeaturePlot(seu, 
                     features = c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act'), 
                     reduction = "umap", 
                     min.cutoff = "q1", 
                     max.cutoff = "q99", 
                     ncol = 3),
  width = 12, height = 10, dpi = 300
)

# ========================================
# Save Stacked Bar Charts
# ========================================

# 11. Cluster composition by condition
stacked_data_condition <- seu@meta.data %>%
  select(condition, Annotated_Cluster) %>%
  group_by(condition, Annotated_Cluster) %>%
  summarise(Cell_Count = n(), .groups = 'drop_last') %>%
  mutate(Fraction = Cell_Count / sum(Cell_Count)) %>%
  ungroup()

ggsave(
  filename = "Stacked_bar_by_condition.png",
  plot = ggplot(stacked_data_condition, aes(x = condition, y = Fraction, fill = Annotated_Cluster)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Cluster Composition by Condition",
         x = "Condition", y = "Fraction of Cells", fill = "Annotated Cluster") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  width = 10, height = 7, dpi = 300
)

# Print confirmation
print("✅ All plots saved successfully!")
print(paste("Location:", getwd()))
print("Files created:")
print("  - UMAP_by_condition.png")
print("  - UMAP_by_clusters.png")
print("  - UMAP_annotated_clusters.png")
print("  - UMAP_split_by_condition.png")
print("  - DotPlot_by_cluster.png")
print("  - DotPlot_by_annotated_cluster.png")
print("  - DotPlot_by_condition.png")
print("  - Coxcomb_by_feature.png")
print("  - Nightingale_rose_by_cluster.png")
print("  - FeaturePlot_all_metabolic_features.png")
print("  - Stacked_bar_by_condition.png")
