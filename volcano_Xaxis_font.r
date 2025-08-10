library(ggplot2)
library(dplyr)
library(argparse)
library(ggrepel)

make_volcano <- function(res, height, width, outfile, degsfile, title, upcol, downcol, pcut, foldcut, size, user_genes, subset) {

  data <- read.delim(degsfile, header = TRUE, stringsAsFactors = FALSE)

  if (subset != "None" && file.exists(subset)) {
    subset_genes <- readLines(subset)
    data <- data[data$GeneSymbol %in% subset_genes, ]
  }

  pvalue_cutoff <- pcut
  fold_change_cutoff <- foldcut

  data$category <- ifelse(data$log2FoldChange > fold_change_cutoff & data$padj < pvalue_cutoff, "Upper right",
                 ifelse(data$log2FoldChange < -fold_change_cutoff & data$padj < pvalue_cutoff, "Upper left",
                 ifelse(data$log2FoldChange < fold_change_cutoff & data$log2FoldChange >= 0 & data$padj < pvalue_cutoff, "Middle right",
                 ifelse(data$log2FoldChange > -fold_change_cutoff & data$log2FoldChange < 0 & data$padj < pvalue_cutoff, "Middle left", "Other"))))

  # Handle cases where no genes are categorized as "Upper right" or "Upper left"
  upregulated_count <- sum(data$category == "Upper right", na.rm = TRUE)
  downregulated_count <- sum(data$category == "Upper left", na.rm = TRUE)

  upreg_label <- paste("Upregulated (", upregulated_count, ")", sep = "")
  downreg_label <- paste("Downregulated (", downregulated_count, ")", sep = "")

  max_log2fc <- max(abs(data$log2FoldChange), na.rm = TRUE)

  top_upregulated <- data[data$category == "Upper right", ]
  top_upregulated <- top_upregulated[order(-abs(top_upregulated$log2FoldChange)), ][1:5, ]
  
  top_downregulated <- data[data$category == "Upper left", ]
  top_downregulated <- top_downregulated[order(-abs(top_downregulated$log2FoldChange)), ][1:5, ]

  label_genes <- rbind(top_upregulated, top_downregulated)

  if (!is.null(user_genes)) {
    user_genes_df <- data[data$Gene %in% user_genes, ]
    label_genes <- rbind(label_genes, user_genes_df)
  }

  # Filter for significant genes (padj <= 0.05)
  significant_data <- data %>% filter(padj <= 0.05)

  # Calculate the maximum absolute log2FoldChange for significant genes
  max_log2fc <- max(abs(significant_data$log2FoldChange), na.rm = TRUE)

  # Set symmetrical x-axis limits based on max absolute log2FoldChange
  x_limit_left <- -max_log2fc  
  x_limit_right <- max_log2fc  

  print(paste("x_limit_left is:", x_limit_left, "x_limit_right is:", x_limit_right))

  # Create the volcano plot
  volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = category, alpha = category), size = 2) +  # Modify point size if needed
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = "black") +  
    geom_vline(xintercept = c(-fold_change_cutoff, fold_change_cutoff), linetype = "dashed", color = "black") + 
    scale_x_continuous(limits = c(x_limit_left, x_limit_right)) +  
    #labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)", title = title) +
    labs(x = expression(log[2]*"(Fold Change)"),y = expression("-" * log[10]*"(Adjusted p-value)"),title = title) + ##title
    #labs(x = expression(log[2]*"(Fold Change)"),y = expression("-" * log[10]*"(Adjusted p-value)",title = title)) + ##no title 


    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      #plot.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      plot.background = element_blank(),    # Remove plot background
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Arial"),    # <== Title font (Arial, size 20)
      axis.title.x = element_text(size = 19, family = "Arial"),                                  # <== X-axis label font (Arial, size 16)
      axis.title.y = element_text(size = 19, family = "Arial"),                                  # <== Y-axis label font (Arial, size 16)
      axis.text.x = element_text(size = 18),  # Increase size of X-axis numbers
      axis.text.y = element_text(size = 18),  # Increase size of Y-axis numbers
      legend.position = "bottom",
      legend.title = element_blank(),
      text = element_text(family = "Arial")                                                     # <== Default text font (Arial)
    ) +

        
    scale_color_manual(values = c("Upper left" = downcol, "Middle left" = downcol, "Upper right" = upcol, "Middle right" = upcol, "Other" = 'gray'),
                      labels = c("Upper left" = downreg_label, "Upper right" = upreg_label),
                      breaks = c("Upper left", "Upper right")) +
    scale_alpha_manual(values = c("Upper left" = 1, "Middle left" = 0.2, "Upper right" = 1, "Middle right" = 0.2),
                      guide = "none") +
    geom_text_repel(data = label_genes, aes(label = GeneSymbol), 
                    box.padding = 0.5, point.padding = 0.1,
                    segment.size = 0.2, segment.color = "grey50",
                    size = 5,                # <== Increase the font size for gene labels here (size=5)
                    color = "black") +
    guides(color = guide_legend(override.aes = list(alpha = 1),label.theme = element_text(size = 14))) ###Downregulated / Upregulated fomnts
  

  tiff(paste0(outfile, ".volcano.tiff"), width = width, height = height, units = 'in', res = res)
  print(volcano_plot)
  dev.off()

}

main <- function() {
 
  parser <- ArgumentParser(description = 'Volcano Plot \n')

  parser$add_argument('-v_res', help = 'Output resolution, as an integer (ex: 150)', default = 300, type = "integer")
  parser$add_argument('-v_ht', help = 'Height of graph, as an integer (ex: 5)', default = 6, type = "integer")
  parser$add_argument('-v_w', help = 'Width of graph, as an integer (ex: 7)', default = 8, type = "integer")
  parser$add_argument('-v_o', help = 'Output file name', default = 'volcanoplot', type = "character")

  parser$add_argument('-v_d', help = 'Path to directory for DEGS file', required = TRUE, type = "character")

  parser$add_argument('-v_t', help = 'Graph title', default = "Volcano Plot", type = "character")
  parser$add_argument('-v_upcol', help = 'Color for upregulated points [default is red]', default = 'red', type = "character")
  parser$add_argument('-v_downcol', help = 'Color for downregulated points [default is blue]', default = 'blue', type = "character") 

  parser$add_argument('-v_pcut', help = 'Cutoff for p-value as a float value [default is 0.05]', default= 0.05, type = 'double') 
  parser$add_argument('-v_foldcut', help = 'Cutoff for log fold change as a float value [default is 0.263]', default = 0.263, type = 'double')
  parser$add_argument('-v_ptsz', help = 'Point size [default is 2]', default = 2, type = 'integer')

  parser$add_argument('-v_il', help = 'Optional desired genes to be labeled, name must match exactly to data sheet', nargs='+', type = "character")

  parser$add_argument('-v_gfp', help = 'Optional file path to data containing subset list of genes to be plotted', default = "None", type = "character")


  args <- parser$parse_args()
  make_volcano(
    res = args$v_res,
    height = args$v_ht,
    width = args$v_w,
    outfile = args$v_o,
    degsfile = args$v_d,
    title = args$v_t,
    upcol = args$v_upcol,
    downcol = args$v_downcol,
    pcut = args$v_pcut,
    foldcut = args$v_foldcut,
    size = args$v_ptsz,
    user_genes = args$v_il,
    subset = args$v_gfp
  )
}

main()
