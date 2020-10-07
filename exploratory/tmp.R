
library(magrittr)
library(ggplot2)
library(autonomics)

object <- dator:::data_sirt
object <- dator:::data_dp



### add_cluster
cluster_type = c("pval", "pbins", "bins", "custom")
contrasts = NULL
breaks = c(0.05, 0.10, 0.20)
sig_level = 0.05
k = 5


### add_ora
db = c("wiki", "kegg", "go_bp", "go_mf", "go_cc", "msigH", "msigC2", "msigC5", "custom")

cluster_type = c("pval", "pbins", "bins")
contrasts = NULL
sig_level = 0.05
p_cutoff = 0.05
q_cutoff = 0.3
simplify_go = FALSE
breaks = c(0.05, 0.1, 0.2)
k = 5
custom_db = NULL


### plot_ora
plot_type = NULL
cluster_type = NULL
contrasts = NULL
db = NULL
exclude = NULL
sig_level = 0.05
p_cutoff = 0.05
q_cutoff = 0.5



# presentation ------------------------------------------------------------

contrasts <- "Sirt_vs_Mock"
plot_type <- "barplot"
db <- "wiki"
cluster_type <- "bins"
p_cutoff <- 0.05


object <- dator:::data_sirt

object %<>% dator::add_ora(db = "wiki", contrasts = "Sirt_vs_Mock")

object %>% dator::plot_ora(db = "wiki", contrasts = "Sirt_vs_Mock", cluster_type = "bins", plot_type = "barplot", p_cutoff = 0.01)
object %>% dator::plot_ora(db = "wiki", contrasts = "Sirt_vs_Mock", cluster_type = "pbins", plot_type = "barplot", p_cutoff = 0.01)



object %<>% dator::add_ora(db = "kegg", contrasts = "Sirt_vs_Mock", cluster_type = "pbins", breaks = c(0.01, 0.05, 0.10))

object %>% dator::plot_ora(db = "kegg", contrasts = "Sirt_vs_Mock", cluster_type = "pbins", plot_type = "barplot")

object %>% dator::plot_ora(db = "kegg", contrasts = "Sirt_vs_Mock", cluster_type = "pbins", plot_type = "dotplot1")

object %<>% dator::add_ora(db = "kegg", contrasts = "Sirt_vs_Mock", cluster_type = "pbins", breaks = c(0.05, 0.1))
object %>% dator::plot_ora(db = "kegg", contrasts = "Sirt_vs_Mock", cluster_type = "pbins", plot_type = "dotplot2")

object %<>% dator::add_ora(db = "kegg", contrasts = "Sirt_vs_Mock", cluster_type = "bins", k = 8)
object %>% dator::plot_ora(db = "kegg", contrasts = "Sirt_vs_Mock", cluster_type = "bins", plot_type = "heatmap")


object %<>% dator::add_ora(db = "go_bp", contrasts = "Sirt_vs_Mock", cluster_type = "pval", sig_level = 0.05, simplify_go = FALSE)
object %>% dator::plot_ora(db = "go_bp", contrasts = "Sirt_vs_Mock", cluster_type = "pval", plot_type = "barplot")

object %<>% dator::add_ora(db = "go_bp", contrasts = "Sirt_vs_Mock", cluster_type = "pval", sig_level = 0.01, simplify_go = FALSE)
object %>% dator::plot_ora(db = "go_bp", contrasts = "Sirt_vs_Mock", cluster_type = "pval", plot_type = "barplot")

# custom cluster from heatmap ---------------------------------------------

plot_type = NULL
cluster_type = "custom"
contrasts = "Sirt_vs_Mock"
db = "kegg"
exclude = NULL
sig_level = 0.05
p_cutoff = 0.05
q_cutoff = 0.5

plot_heatmap(object, contrasts = "Sirt_vs_Mock", p_cutoff = 0.001, scale = TRUE, cluster_columns = FALSE,
             row_km = 3, name=" ", show_row_names = FALSE, border = TRUE,
             breaks = c(-1.5, 0, 1.5),
             colors = c('#046D68', "#E8E8E8", '#B13407')) -> hmap
hmap

object %<>% dator::add_ora(db = "kegg", contrasts = "Sirt_vs_Mock", cluster_type = "custom", custom_cluster = hmap@layout$custom_clust)

object %>% dator::plot_ora(db = "kegg", contrasts = "Sirt_vs_Mock", cluster_type = "custom", plot_type = "dotplot1")




# tests custom db ---------------------------------------------------------


autonomics.import::fdata(object) %>%
    dplyr::select(hydrophobicity, ENTREZID) %>%
    `<<-`(custom_db, .)
autonomics.import::fdata(object) %>%
    dplyr::select(pI, ENTREZID) %>%
    `<<-`(custom_db, .)
autonomics.import::fdata(object) %>%
    dplyr::select(hmoment, ENTREZID) %>%
    `<<-`(custom_db, .)
autonomics.import::fdata(object) %>%
    dplyr::select(instability, ENTREZID) %>%
    `<<-`(custom_db, .)


object %<>% add_ora(db = "custom", contrasts = "dp_vs_up", custom_db = custom_db)

object %>% plot_ora(db = "custom", contrasts = "dp_vs_up", plot_type = "heatmap", cluster_type = "bins", p_cutoff = 0.3, q_cutoff = 1)

object@metadata[["ora_custom"]][["dp_vs_up"]][["bins"]] %>%
    .[order(.$ID), ] %>%
    dplyr::select("Description", "logp", "clust_bins") %>%
    dplyr::mutate(logp = ifelse(logp >= 10, 10, logp)) %>%
    reshape2::dcast(as.formula(paste("Description", "clust_bins", sep="~")),
                    value.var = "logp") %>%
    replace(is.na(.), 0) %>%
    tibble::column_to_rownames(var = "Description") ->tmp


p <- tmp %>% pheatmap::pheatmap(cluster_cols = FALSE,
                                cluster_rows = FALSE,
                                #labels_col = label2,
                                angle_col = 0,
                                color = colorRampPalette(
                                    RColorBrewer::brewer.pal(n = 7, name = "PuBu"))(100),
                                silent = TRUE)
p


