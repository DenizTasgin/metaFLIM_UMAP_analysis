library(Seurat)
library(dplyr)
library(reshape2)
library(Rphenograph) #https://github.com/JinmiaoChenLab/Rphenograph
library(igraph) 
library(GGally)
library(uwot)
library(scales) 

####
# Import Spatial Data ####
####

metadata <- read.csv("/Users/deniztasgin/Downloads/R_clustring/R/10_time_point_clustering.csv")

# ===============================================
# FIX: Clean condition names - remove FOV suffixes
# ===============================================
metadata$condition <- gsub("_FOV[0-9]+$", "", metadata$condition)
metadata$condition <- as.factor(metadata$condition)
# ===============================================

metadata$pixels <- as.factor(metadata$pixels)

# ===============================================
# FIX: Include ALL 7 features in Seurat object creation
# Columns 2-8 = free, FAO, anaerob, TCA, aerob, oxburst, act
# ===============================================
rawdata <- t(metadata[,c(2,3,4,5,6,7,8)])
seu <- CreateSeuratObject(rawdata, assay = "Pixel")


####
# phenograph ####
####

# ===============================================
# FIX: Include ALL 7 metabolic features in normalization
# Columns 2-8 = free, FAO, anaerob, TCA, aerob, oxburst, act
# ===============================================
rawdata <- apply(as.matrix(metadata[,c(2,3,4,5,6,7,8)]),2,as.numeric)
rawdata[is.na(rawdata)] <- 0
normdata <- apply(rawdata, 2, function(x) {
    log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x))))
})



# ----------------- FIX: ADD CUSTOM NORM DATA TO SEURAT OBJECT -----------------
# 1. Transpose normdata back to the gene x cell format Seurat expects
normdata_t <- t(normdata) 

# 2. Assign the custom normalized data to the 'data' slot of the 'Pixel' assay
seu[["Pixel"]] <- SetAssayData(object = seu[["Pixel"]], layer = "data", new.data = normdata_t)

# 3. Ensure the cell and feature names match (crucial)
# This usually happens automatically if 'seu' was created from the same data
# and normdata has correctly named columns/rows.

# ----------------- END FIX -----------------



####
## cluster ####
####
set.seed(1)
Rphenograph_out <- Rphenograph(normdata, k = 20)
cluster_results <- igraph::cluster_louvain(Rphenograph_out[[1]], weights = NULL, resolution = 0.2)


####
## train UMAP ####
####

set.seed(1)
normdata_umap <- uwot::umap(normdata)
normdata_umap <- data.frame(normdata_umap)
colnames(normdata_umap) <- c("UMAP_1","UMAP_2")


# ----------------- ADDED CODE: Integrate Results into Seurat Object -----------------

# 1. Add UMAP Coordinates to Seurat Object
# The cell barcodes must match the row names of the UMAP matrix
rownames(normdata_umap) <- colnames(seu) # Assuming the order is already correct (based on your script flow)
umap_matrix <- as.matrix(normdata_umap[, c("UMAP_1", "UMAP_2")])

# Create the DimReducObject and add it
umap_dr <- CreateDimReducObject(
    embeddings = umap_matrix,
    key = "UMAP_",            
    assay = DefaultAssay(seu) 
)
seu[["umap"]] <- umap_dr # Now FeaturePlot can find the 'umap' reduction!

# 2. Add Cluster and Mouse Metadata
seu <- AddMetaData(seu, metadata = metadata$condition, col.name = "condition")
seu <- AddMetaData(seu, metadata = metadata$timepoint, col.name = "timepoint")
seu <- AddMetaData(seu, metadata = metadata$dataset, col.name = "dataset")
seu <- AddMetaData(seu, metadata = metadata$pixels, col.name = "pixels")
seu <- AddMetaData(seu, metadata = factor(cluster_results$membership), col.name = "pheno_cluster")

# ----------------- END ADDED CODE -----------------



# visualize umap

normdata_umap$condition <- metadata$condition # $condition 
normdata_umap$timepoint <- metadata$timepoint
normdata_umap$dataset <- metadata$dataset
normdata_umap$pixels <- metadata$pixels

# condition 
ggplot(normdata_umap, mapping = aes(x = UMAP_1, y = UMAP_2, colour = condition)) + 
    geom_point() 

# timepoint
ggplot(normdata_umap, mapping = aes(x = UMAP_1, y = UMAP_2, colour = timepoint)) + 
    geom_point()

# dataset
ggplot(normdata_umap, mapping = aes(x = UMAP_1, y = UMAP_2, colour = dataset)) + 
    geom_point()

# pixles
ggplot(normdata_umap, mapping = aes(x = UMAP_1, y = UMAP_2, colour = pixels)) + 
    geom_point()

# leiden clusters
ggplot(normdata_umap, mapping = aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster_results$membership))) + 
    geom_point() + guides(colour = guide_legend(title = "clusters"))




####
### dot plots ####
####

# Dotplot cluster
Seurat::DotPlot(seu, 
                features = c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act'), 
                group.by = "pheno_cluster", # Use the new metadata column name
                cols = c("blue","orange"), 
                assay = "Pixel"
) +
    RotatedAxis() + 
    coord_flip() + 
    theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=8)) + 
    NoLegend()
########################################################################################

#Coxcomb analogon

# 1. Define the features and the clustering variable
features_to_plot <- c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act')
clustering_column <- "pheno_cluster" # Replace with the actual metadata column name

# 2. Use Seurat's functions to calculate the DotPlot metrics
# The built-in DotPlot function uses these metrics, but you have to extract them.
# The following code is an efficient way to get the data structure you need:

coxcomb_data <- Seurat::DotPlot(
    object = seu,
    features = features_to_plot,
    group.by = clustering_column,
    # This makes the data frame tidy (long format)
    cols = c("lightgrey", "blue"), # Dummy colors, just for calculation
    dot.scale = 6
)$data %>%
    dplyr::rename(
        Cluster = id,
        Feature = features.plot,
        avg_exp = avg.exp,
        pct_exp = pct.exp
    )

coxcomb_data <- coxcomb_data %>%
    # Ensure the data is ordered correctly for the polar plot
    arrange(Feature, Cluster) %>%
    # Map a categorical variable (Cluster) to the X-axis for the angle slices
    mutate(x_axis_var = Cluster)

p_coxcomb_dotplot <- ggplot(
    data = coxcomb_data,
    aes(
        x = x_axis_var,          # Mapped to Cluster (angle)
        y = avg_exp,             # Mapped to avg_exp (radius)
        fill = pct_exp           # Mapped to pct_exp (color)
    )
) +
    # Create the bars (wedges)
    geom_bar(width = 1, stat = "identity", color = "black", linewidth = 0.3) +

    # --- CRITICAL: Apply polar coordinates ---
    coord_polar(theta = "x", start = 0) +

    # --- CRITICAL: Facet by Feature ---
    # Each facet is a single Coxcomb plot for one Feature
    facet_wrap(~Feature, ncol = 3) +

    # --- Customization ---
    # Use a continuous color scale (analogous to the color in DotPlot)
    scale_fill_viridis_c(name = "% Cell proportion") +


    #Adjust radial axis (y-axis in Cartesian, radius in polar)
    scale_y_continuous(
        name = "Average Expression",
        limits = c(0, max(coxcomb_data$avg_exp)), # Set max limit
        breaks = seq(0, max(coxcomb_data$avg_exp), by = 0.1) # Example breaks
    ) +

    theme_minimal() +
    theme(
        # Remove titles and labels that are misleading in polar view
        axis.title = element_blank(),
        # Customize the X-axis text (cluster labels)
        axis.text.x = element_text(size = 8, face = "bold"),
        # Remove radial Y-axis ticks/text
        axis.text.y = element_blank(),
        panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
        panel.spacing = unit(1, "lines"), # Add space between facets
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +

    labs(
        title = "Coxcomb Chart: Activity Level (Radius) by Cell Proportion (Color)"
    )

print(p_coxcomb_dotplot)

# Feature Coxcombs for each individual clusters


feature_colors <- c(
    "free" = "#293276",   
    "FAO" = "#2CA9E0",
    "anaerob" = "#4CB648", 
    "TCA" = "#F9E300",    
    "aerob" = "#EF6724",  
    "oxburst" = "#E7138A",
    "act" = "#750079"     
    # ... include all unique levels of your Feature variable
)



coxcomb_data_radial <- coxcomb_data %>%
    #mutate(avg_exp_log = avg_exp + 0.001) %>%#######################  replace zeros
    group_by(Cluster) %>%
    arrange(Cluster, Feature) %>% # Ensure consistent order for stacking
    mutate(
        # Calculate the fraction of the total pct_exp for each cluster
        pct_exp_norm = pct_exp / sum(pct_exp),

        # Calculate start and end angles (in radians, 0 to 2*pi)
        end_angle = 2 * pi * cumsum(pct_exp_norm),
        start_angle = lag(end_angle, default = 0)
    ) %>%
    ungroup()

# --- Plotting the Coxcomb Chart ---
p_radial_coxcomb <- ggplot(
    data = coxcomb_data_radial,
    aes(
        # X and Y are dummy for geom_rect, values for start/end angles and radius are set in geom_rect
        fill = Feature,
        color = Feature # Use same color for line border for consistency
    )
) +
    # Use geom_rect to draw each wedge with custom start/end angles and radius
    geom_rect(
        aes(
            xmin = start_angle,  # Start angle of the wedge
            xmax = end_angle,    # End angle of the wedge
            ymin = 0,            # Min radius
            ymax = avg_exp       # Max radius (avg_exp)
        ),
        linewidth = 0.5
    ) +

    # --- CRITICAL: Apply polar coordinates with theta="x" ---
    # The x-axis is now representing the angle range (0 to 2*pi)
    coord_polar(theta = "x", start = 0) +

    # --- Facet by Cluster ---
    facet_wrap(~Cluster, ncol = 3) +

    # Use the manual scale for the fill aesthetic (Feature)
    scale_fill_manual(
        name = "Feature",
        values = feature_colors
    ) +

    # Use the manual scale for the color aesthetic (Feature border)
    scale_color_manual(
        name = "Feature",
        values = feature_colors
    ) +

    #Adjust radial axis (y-axis in Cartesian, radius in polar)
    scale_y_continuous(
        name = "Average Expression",
        limits = c(0, max(coxcomb_data_radial$avg_exp)), # Set max limit
        breaks = seq(0, max(coxcomb_data_radial$avg_exp), by = 0.1) # Example breaks
    ) +


    # scale_y_continuous(
    #     name = "Average Expression (Log10 Scale)",
    #     # Use log10 for breaks and labels, but NO trans="log10" on the scale itself
    #     breaks = trans_breaks("log10", function(x) 10^x),
    #     labels = function(x) paste0("10^", round(log10(x - 0.001), 1))
    # ) +

    theme_minimal() +
    theme(
        axis.title.x = element_blank(), # Remove X-axis label (angle)
        axis.text.x = element_blank(),  # Remove X-axis text (angle values)
        axis.ticks.x = element_blank(), # Remove X-axis ticks
        axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = -10)), # Custom Y-axis title
        axis.text.y = element_text(size = 8), # Y-axis text (radius values)
        panel.grid.major.x = element_blank(), # Remove major grid lines for X (angle)
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5), # Keep radial grid lines
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
    ) +
    labs(
        title = "Nightingale Rose Chart: Pct Expressing (Angle) & Avg Expression (Radius)"
    )

print(p_radial_coxcomb)

########################################################################################

# Dotplot timepoint
Seurat::DotPlot(seu, 
                features = c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act'), 
                group.by = "timepoint", # Use the new metadata column name
                cols = c("blue","orange"), 
                assay = "Pixel"
) +
    RotatedAxis() + 
    coord_flip() + 
    theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=8)) + 
    NoLegend()

# ===============================================
# FIX: Add 4th color for 4 conditions
# ===============================================
# Dotplot timepoint for different conditions
Seurat::DotPlot(seu, 
                features = c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act'), 
                group.by = "timepoint", # Use the new metadata column name
                split.by = "condition", 
                cols = c("orange","cyan", "blue", "red"), # Added 4th color
                assay = "Pixel"
) +
    RotatedAxis() + 
    coord_flip() + 
    theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=8)) + 
    NoLegend()


# Dotplot condition
Seurat::DotPlot(seu, 
                features = c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act'), 
                group.by = "condition", # Use the new metadata column name
                #cols = c("blue","red"), 
                assay = "Pixel"
) +
    RotatedAxis() + 
    coord_flip() + 
    theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=8)) + 
    NoLegend()


####
### feature plots ####
####


# Feature plots
FeaturePlot(seu, 
            features = c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act'), 
            reduction = "umap", 
            min.cutoff = "q1", 
            max.cutoff = "q99", 
            ncol = 3
)




####
### violin points ####
####

# rawdata violinplot (variable vs timepoint)
# rawdata_vis <- data.frame(rawdata, clusters = factor(cluster_results$membership), 
#                           condition = metadata$condition, timepoint = metadata$timepoint)
# rawdata_melt_vis <- reshape2::melt(rawdata_vis, id.vars = c("clusters", "condition", "timepoint"))
# ggplot(rawdata_melt_vis, mapping = aes(x = variable, y = value, color = variable, fill = variable)) + 
#     geom_point(position = position_jitter()) + 
#     facet_grid(~timepoint) + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 


# rawdata violinplot (variable vs condition)
# rawdata_vis <- data.frame(rawdata, clusters = factor(cluster_results$membership), 
#                           condition = metadata$condition, dataset = metadata$dataset)
# rawdata_melt_vis <- reshape2::melt(rawdata_vis, id.vars = c("clusters", "condition", "dataset"))
# ggplot(rawdata_melt_vis, mapping = aes(x = variable, y = value, color = variable, fill = variable)) + 
#     geom_point(position = position_jitter()) + 
#     facet_grid(~condition) + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# 


####



### annotation of clusters ####
####

# ===============================================
# FIX: Handle all 21 clusters from Rphenograph
# First, check how many clusters we actually have
# ===============================================
print("Detected cluster IDs:")
print(unique(seu$pheno_cluster))
print(paste("Total number of clusters:", length(unique(seu$pheno_cluster))))

# 1. CRITICAL: Set the pheno_cluster column as the active identity.
# This ensures RenameIdents knows WHICH cluster labels to look for.
Idents(seu) <- "pheno_cluster"

# 2. Create annotation for ALL clusters
# Map the 21 clusters to meaningful categories
# You'll need to adjust these based on your biological interpretation
new_cluster_names <- c(
    "1" = "Canonical", 
    "2" = "Canonical", 
    "3" = "Canonical",
    "4" = "AAM_dominant",
    "5" = "AAM_dominant",
    "6" = "Cytosolic_glycolysis",
    "7" = "Canonical",
    "8" = "Unassigned_8",
    "9" = "Unassigned_9",
    "10" = "Unassigned_10",
    "11" = "Unassigned_11",
    "12" = "Unassigned_12",
    "13" = "Unassigned_13",
    "14" = "Unassigned_14",
    "15" = "Unassigned_15",
    "16" = "Unassigned_16",
    "17" = "Unassigned_17",
    "18" = "Unassigned_18",
    "19" = "Unassigned_19",
    "20" = "Unassigned_20",
    "21" = "Unassigned_21"
)


# 3. Apply RenameIdents and store the RESULTING SEURAT OBJECT.
seu <- RenameIdents(seu, new_cluster_names)

# 4. Store the new identity labels in a NEW metadata column for stability:
seu$Annotated_Cluster <- Idents(seu) 

# View the result with Seurat's DimPlot (using your UMAP reduction)
DimPlot(seu, reduction = "umap", group.by = "Annotated_Cluster", , pt.size = 1.5, label = TRUE) + theme_gray() +
    ggtitle("UMAP with Annotated Clusters")


# View the result with Seurat's DimPlot split by condition (using your UMAP reduction)
DimPlot(seu, reduction = "umap", group.by = "Annotated_Cluster", split.by = "condition", pt.size = 1.5, label = FALSE) + theme_gray() +
    ggtitle("UMAP with Annotated Clusters split by condition")


# p <- DimPlot(seu, reduction = "umap", group.by = "Annotated_Cluster", split.by = "condition", pt.size = 1, label = FALSE) + theme_gray()
# ggsave(
#     filename = "my_dimplot.png",
#     plot = p,
#     width = 30, # Width in units (e.g., inches)
#     height = 12, # Height in units (e.g., inches)
#     units = "cm", # Units can be "in", "cm", or "mm"
#     dpi = 300 # Resolution (dots per inch)
# )

features_list <- c('free','FAO','anaerob','TCA', 'aerob', 'oxburst', 'act')

# Dotplot Annotated cluster (NEW p6)
Seurat::DotPlot(seu, features = features_list, group.by = "Annotated_Cluster", cols = c("blue","orange"), assay = "Pixel") +
    RotatedAxis() + coord_flip() + theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=8)) + NoLegend() + ggtitle("Dot Plot: Grouped by Annotated Cluster")





####
### fractional bar charts ####
####


# Compare the fractional composition across all 'timepoint' groups

stacked_data <- seu@meta.data %>%
    # Select the two columns of interest: the grouping variable and the clusters
    select(timepoint, Annotated_Cluster) %>%

    # Count the cells in each combination (timepoint x cluster)
    group_by(timepoint, Annotated_Cluster) %>%
    summarise(Cell_Count = n(), .groups = 'drop_last') %>%

    # Calculate the fraction of total cells *within each timepoint*
    mutate(Fraction = Cell_Count / sum(Cell_Count)) %>%
    ungroup()


stacked_composition_plot <- ggplot(stacked_data, aes(x = timepoint, y = Fraction, fill = Annotated_Cluster)) +
    geom_bar(stat = "identity", position = "stack") +

    # Add labels to the plot
    labs(
        title = "Cluster Composition by Timepoint",
        x = "Timepoint",
        y = "Fraction of Cells",
        fill = "Annotated Cluster"
    ) +
    scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(stacked_composition_plot)




# Compare the fractional composition across all 'condition' groups

stacked_data <- seu@meta.data %>%
    # Select the two columns of interest: the grouping variable and the clusters
    select(condition, Annotated_Cluster) %>%

    # Count the cells in each combination (condition x cluster)
    group_by(condition, Annotated_Cluster) %>%
    summarise(Cell_Count = n(), .groups = 'drop_last') %>%

    # Calculate the fraction of total cells *within each condition*
    mutate(Fraction = Cell_Count / sum(Cell_Count)) %>%
    ungroup()


stacked_composition_plot <- ggplot(stacked_data, aes(x = condition, y = Fraction, fill = Annotated_Cluster)) +
    geom_bar(stat = "identity", position = "stack") +

    # Add labels to the plot
    labs(
        title = "Cluster Composition by condition",
        x = "condition",
        y = "Fraction of Cells",
        fill = "Annotated Cluster"
    ) +
    scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(stacked_composition_plot)



####
### Coxcomb charts ####
####
# 
# M0_counts <- stacked_data %>%
#     filter(condition == "M0") %>%
#     # Ensure the factor levels are correct after filtering
#     droplevels()
# 
# M1_counts <- stacked_data %>%
#     filter(condition == "M1") %>%
#     # Ensure the factor levels are correct after filtering
#     droplevels()
# 
# M2_counts <- stacked_data %>%
#     filter(condition == "M2") %>%
#     # Ensure the factor levels are correct after filtering
#     droplevels()
# 
# p_M0_coxcomb <- ggplot(M0_counts, 
#     aes(x = Annotated_Cluster, y = Fraction, fill = Annotated_Cluster)
# ) +
#     # Use geom_bar to create the wedges. stat="identity" tells it to use the 'y' value directly.
#     geom_bar(width = 1, stat = "identity") +                #, color = "white", linewidth = 0.5
#     
#     # --- CRITICAL: Apply polar coordinates ---
#     coord_polar(theta = "x", start = 0) +
#     
#     # --- Customization ---
#     #scale_fill_brewer(palette = "Set1") + # Use a nice color palette
#     theme_minimal() +
#     
#     # Remove unnecessary labels and background elements
#     theme(
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"),
#         panel.grid.minor = element_blank(),
#         axis.text.y = element_blank(), # Remove Y-axis labels (radial tick marks)
#         axis.text.x = element_text(size = 10, face = "bold"), # Keep cluster labels
#         legend.position = "right",
#         plot.title = element_text(hjust = 0.5, face = "bold")
#     ) 
#     
#     # Add titles and labels
#     # labs(
#     #     title = "Coxcomb Chart of metabolic Cluster Cell Counts",
#     #     fill = "Cluster ID"
#     # ) 
#     
#     # # Add labels for the count on the outside of the chart
#     # geom_text(
#     #     aes(label = Fraction), 
#     #     y = stacked_data$Fraction + (max(stacked_data$Fraction) * 0.05), # Position slightly outside the bar
#     #     # Use polar coordinates to adjust text position:
#     #     position = position_stack(vjust = 0.5), # Center text if bars were stacked, but here used for polar
#     #     inherit.aes = TRUE, 
#     #     size = 3
#     # )
# 
# # Print the plot
#     print(p_M0_coxcomb)
# 
# p_M1_coxcomb <- ggplot(M1_counts, 
#                        aes(x = Annotated_Cluster, y = Fraction, fill = Annotated_Cluster)
# ) +
#     # Use geom_bar to create the wedges. stat="identity" tells it to use the 'y' value directly.
#     geom_bar(width = 1, stat = "identity") +                #, color = "white", linewidth = 0.5
#     
#     # --- CRITICAL: Apply polar coordinates ---
#     coord_polar(theta = "x", start = 0) +
#     
#     # --- Customization ---
#     #scale_fill_brewer(palette = "Set1") + # Use a nice color palette
#     theme_minimal() +
#     
#     # Remove unnecessary labels and background elements
#     theme(
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"),
#         panel.grid.minor = element_blank(),
#         axis.text.y = element_blank(), # Remove Y-axis labels (radial tick marks)
#         axis.text.x = element_text(size = 10, face = "bold"), # Keep cluster labels
#         legend.position = "right",
#         plot.title = element_text(hjust = 0.5, face = "bold")
#     ) 
#     
#     # Add titles and labels
#     # labs(
#     #     title = "Coxcomb Chart of metabolic Cluster Cell Counts",
#     #     fill = "Cluster ID"
#     # ) 
#     
#     # # Add labels for the count on the outside of the chart
#     # geom_text(
#     #     aes(label = Fraction), 
#     #     y = stacked_data$Fraction + (max(stacked_data$Fraction) * 0.05), # Position slightly outside the bar
#     #     # Use polar coordinates to adjust text position:
#     #     position = position_stack(vjust = 0.5), # Center text if bars were stacked, but here used for polar
#     #     inherit.aes = TRUE, 
#     #     size = 3
#     # )
#     
#     # Print the plot
#     print(p_M1_coxcomb)
# 
# p_M2_coxcomb <- ggplot(M2_counts, 
#                        aes(x = Annotated_Cluster, y = Fraction, fill = Annotated_Cluster)
# ) +
#     # Use geom_bar to create the wedges. stat="identity" tells it to use the 'y' value directly.
#     geom_bar(width = 1, stat = "identity") +                #, color = "white", linewidth = 0.5
#     
#     # --- CRITICAL: Apply polar coordinates ---
#     coord_polar(theta = "x", start = 0) +
#     
#     # --- Customization ---
#     #scale_fill_brewer(palette = "Set1") + # Use a nice color palette
#     theme_minimal() +
#     
#     # Remove unnecessary labels and background elements
#     theme(
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major.y = element_line(color = "gray80", linetype = "dashed"),
#         panel.grid.minor = element_blank(),
#         axis.text.y = element_blank(), # Remove Y-axis labels (radial tick marks)
#         axis.text.x = element_text(size = 10, face = "bold"), # Keep cluster labels
#         legend.position = "right",
#         plot.title = element_text(hjust = 0.5, face = "bold")
#     ) 
#     
#     # Add titles and labels
#     # labs(
#     #     title = "Coxcomb Chart of metabolic Cluster Cell Counts",
#     #     fill = "Cluster ID"
#     # ) 
#     
#     # # Add labels for the count on the outside of the chart
#     # geom_text(
#     #     aes(label = Fraction), 
#     #     y = stacked_data$Fraction + (max(stacked_data$Fraction) * 0.05), # Position slightly outside the bar
#     #     # Use polar coordinates to adjust text position:
#     #     position = position_stack(vjust = 0.5), # Center text if bars were stacked, but here used for polar
#     #     inherit.aes = TRUE, 
#     #     size = 3
#     # )
#     
#     # Print the plot
#     print(p_M2_coxcomb)






####
### Save #####
####

# ===============================================
# FIX: Update save paths to Mac-style
# ===============================================
seurat_save_path <- "/Users/deniztasgin/Downloads/R_clustring/R/seu.rds" 
csv_save_path <- "/Users/deniztasgin/Downloads/R_clustring/R/Annotated_Cluster.csv" 
# ===============================================

# ===============================================
## Saving the Seurat Object
# ===============================================

# It's best practice to save the entire Seurat object for future use
saveRDS(seu, file = seurat_save_path)
print(paste("Seurat object saved to:", seurat_save_path))


cell_ids <- colnames(seu)

cluster_annotations <- seu@meta.data$Annotated_Cluster

annotations_df <- data.frame(
    CellID = cell_ids,
    Cluster = cluster_annotations
)

# 4. Save the Data Frame to a CSV File
write.csv(annotations_df, file = csv_save_path, row.names = FALSE)

print(paste("Cluster annotations saved to:", csv_save_path))

# Display the first few rows of the saved data frame for verification
print(head(annotations_df))
