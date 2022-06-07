library(dplyr)
library(ggplot2)
library(ggrepel)

d <- data.table::fread("Results_DGE_of_Treatment_vs_Control_shrinkage.csv") %>% as.data.frame() ## Load DGE data from L2FC shrinkage

## Select genes to show on the volcano plot
showgenes <- c("Fos","Met","Gata1","Ccr9") 
## Change font of label
font_size <- 4 
## Create a cap padj or pvalue.  The value is based on -log10(value).  You can use -log10(desired value) if preferred.  This cap is necessary if there are genes that just so high that it significantly shrink the plot
padj_cap_value <- 10 
## Similar to capping the significant value.  This cap significantly large FC when comparing to the other genes
l2fc_cap_value <- 3
## Create a minimum cut off for l2fc.  Can be set to 0 if you want all genes that is padj/pvalue significant
l2fc_cutoff_value <- 0.25
## Choose either padj to pvalue to use for graphing
signif_factor <- "padj"
## Select a low range color (Choices can be obtain here: https://r-charts.com/colors/)
low_color = "blue1"
## Select a high range color
high_color = "red1"
## Select mid range or near 0 L2FC
mid_color = "grey60"


processed_d <-d %>%
## create a column with capped -log10(padj) to contain extremely large value
  mutate(padj_cap = ifelse((-log10(get(signif_factor)))>padj_cap_value, padj_cap_value,-log10(get(signif_factor)))) %>% 
## Selecting cutoff FC and significant factor either pvalue or padj
  mutate(selected = ifelse((abs(log2FoldChange) > l2fc_cutoff_value)&(get(signif_factor) < .05), log2FoldChange,NA)) %>%
## Capping off FC
  mutate(selected = ifelse(selected < -l2fc_cap_value, -l2fc_cap_value, ifelse(selected > l2fc_cap_value, l2fc_cap_value, selected))) %>%
  as.data.frame()

ggplot(processed_d,aes(log2FoldChange,padj_cap))+
## Graping the geom_point of selected genes that are significant and coloring appropriately
  geom_point(aes(color = selected),size=1.5)+
  scale_color_gradient2(low = low_color,midpoint = 0, high = high_color, mid = mid_color,limits = c(-3,3), na.value="grey90")+
## Label for the axis.  Can be manually changed here.
  xlab("Fold Change (log2)")+
  ylab(paste(signif_factor,"(-log10)"))+
## Label and designate X axis (FC) based on range within given cap
  scale_x_continuous(breaks = c(-l2fc_cap_value:l2fc_cap_value), limits = c(-l2fc_cap_value,l2fc_cap_value))+
## Draw dashed horizontal and vertical lines to signify cut off from pvalue/padj and FC
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = -l2fc_cutoff_value, linetype = "dashed",alpha = ifelse(l2fc_cutoff_value > 0, 1, 0))+
  geom_vline(xintercept = l2fc_cutoff_value,linetype = "dashed",alpha = ifelse(l2fc_cutoff_value > 0, 1, 0))+
## Annotate the selected genes input from the showgenes onto the graph.
  geom_text_repel(data = subset(processed_d, symbol %in% showgenes),
                  aes(log2FoldChange,padj_cap, label = symbol),
                  size = font_size,
                  fontface = "bold",
                  nudge_y = .5,
                  segment.size = 1.5,
                  nudge_x=0,
                  min.segment.length = 0.2)+
## Minor thematic alteration.  Change for your desired parameter
  theme(panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(size=1, color = "black", fill = NA),
        axis.text = element_text(size=20, color = "black"),
        axis.title = element_text(size=20 , color = "black"),
        text = element_text(size = 10),
        aspect.ratio = 1/1)

ggsave("Treatment_vs_Control_VolcanoPlot.jpg",plot = p,device = "jpg", width = 1600, height = 1600, units = "px")
## Often time, the labels and segment generated from geom_text_repel are not in good position and manual curation such as nudging parameters and splitting up the geom_text_repel to each individual gene
## However, I found it more convenient to just spit it back as svg graphics and edit it manually with software like Inkscape or SVG editor and move the line and text around
## Alternatively, just omit the showgenes and create an empty volcano plot to manually add the gene in later

ggsave("Treatment_vs_Control_VolcanoPlot.svg",plot = p,device = "svg", width = 1600, height = 1600, units = "px") 

## dplyr
## Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2022). dplyr: A Grammar of Data Manipulation. R package version
## 1.0.8. https://CRAN.R-project.org/package=dplyr

## ggplot2
## H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

## ggrepel
## Kamil Slowikowski (2021). ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'. R package version 0.9.1.
## https://CRAN.R-project.org/package=ggrepel
