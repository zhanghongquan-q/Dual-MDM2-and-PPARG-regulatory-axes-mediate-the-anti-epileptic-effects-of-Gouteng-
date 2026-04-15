library(GEOquery)
library(limma)
library(affy)
library(Biobase)
library(annotate)
library(hugene10sttranscriptcluster.db)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(stringr)
library(tibble)
library(clusterProfiler)
library(mouse430a2.db)
library(ggplot2)
gset <- getGEO("GSE88992",destdir=".",AnnotGPL=FALSE,getGPL=FALSE)
class(gset)
gset[[1]]
rt <-pData(gset[[1]])
table(rt$title)
group_list <- ifelse(str_detect(rt$title," saline"),"Health","EP")
rt$group <-group_list
rt1<- rt %>% dplyr::select("group")
write.table(rt1,file = "group_list.txt",sep="\t",row.names = T,col.names=NA,quote =F)

exp<- exprs(gset[[1]])
exp_normalized <- normalizeBetweenArrays(exp)


ids <- toTable(mouse430a2SYMBOL)
length(unique(ids$symbol))

exp_normalized1 <- as.data.frame(exp_normalized)

exp_normalized1 <- exp_normalized1 %>%
  mutate(probe_id = rownames(exp_normalized1))

exp_normalized1 <- exp_normalized1 %>%
  inner_join(ids,by = "probe_id")

exp_normalized1 <- exp_normalized1 %>%
  dplyr::select(probe_id,symbol,everything())

exp_normalized1 <- exp_normalized1[!duplicated(exp_normalized1[["symbol"]]), ]

rownames(exp_normalized1) <- exp_normalized1$symbol

exp_normalized1 <- exp_normalized1[,-(1:2)]

identical(rownames(rt1),colnames(exp_normalized1))
group_list <- factor(rt1$group,levels = c("Health","EP"))
design <- model.matrix(~group_list)

fit <- lmFit(exp_normalized1,design)

fit2 <- eBayes(fit)
deg <- topTable(fit2,coef = 2,number = Inf)
DEG = na.omit(deg)

mart_data <- read.table(
  "mart_export.txt", 
  sep = "\t", 
  header = TRUE, 
  stringsAsFactors = FALSE
)

filtered_data <- data.frame(
  Gene.name = mart_data[["Gene.name"]],
  Human.gene.name = mart_data[["Human.gene.name"]],
  Human.homology.type = mart_data[["Human.homology.type"]],
  Human.orthology.confidence = mart_data[["Human.orthology.confidence..0.low..1.high."]]
)
filtered_data <- unique(filtered_data)
exp_normalized1$symbol <- rownames(exp_normalized1)
exp_normalized1 <- merge(
  exp_normalized1, 
  filtered_data, 
  by.x = "symbol", 
  by.y = "Gene.name", 
  all.x = TRUE
)

success_count <- sum(
  exp_normalized1$Human.gene.name != "F" & 
    !is.na(exp_normalized1$Human.gene.name))
fail_count <- nrow(exp_normalized1) - success_count
conversion_rate <- success_count / nrow(exp_normalized1) * 100
stats_data <- data.frame(
  Status = c("success", "fail"),
  Count = c(success_count, fail_count),
  Percentage = c(conversion_rate, 100 - conversion_rate)
)


p = ggplot(stats_data, aes(x = "", y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = sprintf("%.2f%%", round(Percentage, 2))), 
            position = position_stack(vjust = 0.5)) +
  theme_void() +
labs(title = "Mouse-to-Human Gene Homology Conversion Success Rate")
print(p)
ggsave("gene_homology_conversion_pie_chart.tiff", p, width = 10, height = 8, dpi = 600)
dev.off()

if ("Human.orthology.confidence" %in% colnames(exp_normalized1)) {
  exp_normalized1 <- exp_normalized1[order(exp_normalized1$symbol, -exp_normalized1$Human.orthology.confidence), ]  
  exp_normalized1 <- exp_normalized1[!duplicated(exp_normalized1$symbol), ]  
}

exp_normalized1$Human.gene.name <- gsub("[^A-Za-z0-9_-]", "", exp_normalized1$Human.gene.name)
