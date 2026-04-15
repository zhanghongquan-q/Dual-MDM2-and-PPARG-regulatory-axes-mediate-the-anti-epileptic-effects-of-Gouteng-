library(dplyr)        
library(CIBERSORT)    
library(ggplot2)      
library(pheatmap)     
library(ggpubr)       
library(reshape2)
library(tidyr)
library(limma)        
library(ggthemes)
library(RColorBrewer) 
library(corrplot) 
expr <- read.csv("exp2.csv", row.names = 1) 
group <- read.csv("group2.csv", row.names = 1)    
expr <- normalizeBetweenArrays(expr) 
LM22 <- read.table("LM22.txt", header = TRUE, row.names = 1, sep = "\t")
LM22 <- as.matrix(LM22)
mode(LM22) <- "numeric"
expr <- as.matrix(expr)
mode(expr) <- "numeric"
expr <- expr[complete.cases(expr), ]
result <- cibersort(sig_matrix = LM22, mixture_file = expr, perm = 1000, QN = FALSE)
result1 <- result[, 1:22]
non_zero_cols <- colSums(result1 == 0) != nrow(result1)  
result1 <- result1[, non_zero_cols, drop = FALSE]  


cor <- cor(result1) 
corrplot(cor, method = "number") 
corrplot(cor, order = "AOE", type = "upper", tl.pos = "d")


data1 <- cbind(result1, group)  
data1$sample <- rownames(data1)  
data1 <- data1[, c("sample", setdiff(colnames(data1), "sample"))]  
data <- gather(data1, key = CIBERSORT, value = Proportion, -c(group, sample))



data_summary <- data %>%
  group_by(group, CIBERSORT) %>%
  summarise(Mean_Proportion = mean(Proportion, na.rm = TRUE)) %>%
  ungroup()
group_colors <- c("#E69F00", "#56B4E9")  
ggplot(data_summary, aes(x = CIBERSORT, y = Mean_Proportion, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +  
  scale_fill_manual(values = group_colors, name = "Group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "top",  
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = ggplot2::margin(0,0,15,0))
  ) +
  labs(
    title = "Immune Cell Proportion Comparison Between Groups",
    subtitle = "Grouped Bar Plot (Direct Comparison per Cell Type)",
    x = "Immune Cell Type",
    y = "Mean Proportion"
  )

ggplot(data, aes(x = group, y = Proportion, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
stat_compare_means(
  method = "wilcox.test",  
  label = "p.signif",      
  hide.ns = TRUE,          
  size = 5,                
  vjust = 0.8              
) +
facet_wrap(~CIBERSORT, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = ggplot2::margin(0,0,15,0))
  ) +
  labs(
    title = "Immune Cell Proportion Differences Between Groups",
    x = "Group",
    y = "Relative Proportion"
  )
ggsave("Immune_Cell_Boxplot_Comparison.tiff", width = 18, height = 12, dpi = 600, device = "tiff")


key_genes <- c("MDM2", "PPARG") 
key_expr <- expr[rownames(expr) %in% key_genes, , drop = FALSE]  
key_expr_t <- as.data.frame(t(key_expr))
combined_data <- cbind(key_expr_t, result1) 

corr_analysis <- function(gene_name, data, method = "pearson") {
  if (!gene_name %in% colnames(data)) {
    warning(paste("gene", gene_name, "no"))
    return(NULL)
  }
  gene_expr <- data[, gene_name, drop = TRUE]
  if (length(gene_expr) == 0 || all(is.na(gene_expr)) || length(unique(gene_expr)) <= 1) {
    warning(paste("gene", gene_name, "no"))
    return(NULL)
  }
  immune_cols <- !colnames(data) %in% key_genes
  immune_cells <- data[, immune_cols, drop = FALSE]
  if (ncol(immune_cells) == 0) {
    warning("no")
    return(NULL)
  }
  corr_results <- lapply(colnames(immune_cells), function(cell_name) {
    cell <- immune_cells[, cell_name, drop = TRUE]

    if (length(cell) == 0 || all(is.na(cell)) || length(unique(cell)) <= 1) {
      warning(paste("type", cell_name, "no"))
      return(NULL)
    }
    if (length(gene_expr) != length(cell)) {
      warning(paste("gene", gene_name, "cell", cell_name, "no"))
      return(NULL)
    }
    tryCatch({
      test <- cor.test(gene_expr, cell, method = method)
      if (is.null(test$estimate) || is.null(test$p.value)) {
        warning(paste("cell", cell_name, "no"))
        return(NULL)
      }
      data.frame(
        CellType = cell_name,
        Correlation = test$estimate,
        Pvalue = test$p.value,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning(paste("cell", cell_name, "F:", e$message))
      return(NULL)
    })
  })
  valid_results <- Filter(Negate(is.null), corr_results)
  if (length(valid_results) == 0) {
    warning(paste("gene", gene_name, "F"))
    return(NULL)
  }
  
  do.call(rbind, valid_results) %>% 
    arrange(desc(abs(Correlation)))
}

gene_corr_list <- lapply(key_genes, function(gene) {
  result <- corr_analysis(gene, combined_data, method = "pearson")
  if (!is.null(result)) {
    result %>% mutate(Gene = gene)
  } else {
    NULL
  }
})
gene_corr_list <- Filter(Negate(is.null), gene_corr_list)

if (length(gene_corr_list) > 0) {
  all_corr <- do.call(rbind, gene_corr_list)
  all_corr$FDR <- p.adjust(all_corr$Pvalue, method = "fdr")  # BH校正
  print("done")
} else {
  warning("nothing")
}

all_corr_unique <- all_corr %>%
  dplyr::distinct(Gene, CellType, .keep_all = TRUE) %>%
  mutate(Sig = case_when(
    Pvalue < 0.001 ~ "***",
    Pvalue < 0.01  ~ "**",
    Pvalue < 0.05  ~ "*",
    TRUE ~ ""
  ))

sig_marker_df <- all_corr_unique %>%
  select(Gene, CellType, Sig) %>%
  pivot_wider(
    id_cols = Gene,
    names_from = CellType,
    values_from = Sig,
    values_fill = ""
  )
sig_marker <- as.matrix(sig_marker_df[, -1])  
rownames(sig_marker) <- sig_marker_df$Gene     

heatmap_df <- all_corr_unique %>%
  select(Gene, CellType, Correlation) %>%
  pivot_wider(
    id_cols = Gene,
    names_from = CellType,
    values_from = Correlation,
    values_fill = 0
  )

heatmap_data <- as.matrix(heatmap_df[, -1])
rownames(heatmap_data) <- heatmap_df$Gene
sig_marker <- sig_marker[rownames(heatmap_data), colnames(heatmap_data)]

pheatmap(
  mat = heatmap_data,
  display_numbers = sig_marker,
  number_color = "black",
  number_fontsize = 10,
  color = colorRampPalette(brewer.pal(10, "RdBu"))(100),
  scale = "none",
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Correlation Heatmap (P-value Significance Labeled)",
  treeheight_row = 0,
  treeheight_col = 12,
  fontsize_col = 8,
  angle_col = 90
)


sig_strong_pairs <- all_corr_unique %>%
  filter(Pvalue < 0.05 & abs(Correlation) > 0.5) %>%
  arrange(desc(abs(Correlation)))

if (nrow(sig_strong_pairs) == 0) {
  message("no")
} else {
  message(paste("ok", nrow(sig_strong_pairs), "drawing"))
  
  for (i in 1:nrow(sig_strong_pairs)) {
    current_gene <- sig_strong_pairs$Gene[i]
    current_cell <- sig_strong_pairs$CellType[i]
    current_r <- round(sig_strong_pairs$Correlation[i], 2)
    current_pval <- round(sig_strong_pairs$Pvalue[i], 4)
    plot_df <- data.frame(
      Gene_Expression = combined_data[[current_gene]],
      Immune_Abundance = combined_data[[current_cell]]
    )
    p <- ggplot(plot_df, aes(x = Gene_Expression, y = Immune_Abundance)) +
      geom_point(alpha = 0.6, size = 2.5, color = "#2E86AB") +
      geom_smooth(method = "lm", se = TRUE, color = "#A23B72", 
                  linewidth = 1.2, fill = "#F18F01", alpha = 0.3) +
      labs(
        x = paste0(current_gene, " Expression"),
        y = paste0(current_cell, " Abundance"),
        title = paste0(current_gene, " vs ", current_cell),
        subtitle = paste0("Pearson r = ", current_r, " | P = ", current_pval)
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)
      )
    
    print(p)
    ggsave(
      filename = paste0("1", current_gene, "_", gsub(" ", "_", current_cell), ".pdf"),
      plot = p,
      width = 8,
      height = 6,
      dpi = 600
    )
  }
}


combined_data <- combined_data %>%
  mutate(
    MDM2_Group = ifelse(MDM2 > median(MDM2, na.rm = TRUE), "High", "Low"),
    PPARG_Group = ifelse(PPARG > median(PPARG, na.rm = TRUE), "High", "Low")
  )
for (i in 1:nrow(sig_strong_pairs)) {
  current_gene <- sig_strong_pairs$Gene[i]
  current_cell <- sig_strong_pairs$CellType[i]
  group_col <- paste0(current_gene, "_Group") 
  
  p_box <- ggplot(combined_data, aes(x = .data[[group_col]], y = .data[[current_cell]], fill = .data[[group_col]])) +
    geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2, color = "#333333") +
    scale_fill_manual(values = c("#3498DB", "#E74C3C")) + 
    labs(
      x = paste0(current_gene, " Level"),
      y = paste0(current_cell, " Abundance"),
      title = paste0(current_cell, " in ", current_gene, " High vs Low Groups")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none",
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10)
    ) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",    
      label.x = 1.5,
      size = 6,
      vjust = 0.5,
      hide.ns = TRUE         
    )
  
  print(p_box)
  
  ggsave(
    filename = paste0("2", current_gene, "_", gsub(" ", "_", current_cell), ".pdf"),
    plot = p_box,
    width = 7,
    height = 6,
    dpi = 600
  )
}
