library(tidyverse)

expr_path <- "data/GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv"
meta_path <- "data/GSE60450_filtered_metadata.csv"

expr <- read_csv(expr_path, show_col_types = FALSE)
meta <- read_csv(meta_path, show_col_types = FALSE)

expr <- expr %>% rename(ensembl_id = 1, gene_symbol = 2)
meta <- meta %>% rename(sample_id = 1,
                        characteristics = 2,
                        immunophenotype = 3,
                        developmental_stage = 4)

dim(expr)
dim(meta)

gene_of_interest <- "Klf6"   

gene_df <- expr %>%
  filter(gene_symbol == gene_of_interest) %>%
  pivot_longer(cols = -c(ensembl_id, gene_symbol),
               names_to = "sample_id",
               values_to = "expression") %>%
  left_join(meta, by = "sample_id")

stopifnot(nrow(gene_df) > 0)

p <- ggplot(gene_df, aes(x = immunophenotype, y = log2(expression + 1))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
  labs(title = paste0(gene_of_interest, " expression by cell type"),
       x = "Cell type (immunophenotype)",
       y = "log2(CPM + 1)  (TMM normalized)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

dir.create("results", showWarnings = FALSE)
out_file <- file.path("results", paste0(gene_of_interest, "_expression_by_celltype.png"))
ggsave(out_file, plot = p, width = 7, height = 4.5, dpi = 300)

out_file

