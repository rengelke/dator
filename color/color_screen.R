
data(iris)
data(mtcars)

# RColorBrewer ------------------------------------------------------------

library(RColorBrewer)

RColorBrewer::brewer.pal.info %>% 
  .[.$category == "div", ] %>%
  tibble::rownames_to_column(var = "palette") %>%
  apply(., 1, function (x) {
    n <- as.numeric(x["maxcolors"]) - 2
    pal_name <- as.character(x["palette"])
    my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n, pal_name))(50)
    pheatmap::pheatmap(test, color = my_palette, main = pal_name)
    superheat::superheat(mtcars[sample(1:nrow(mtcars), 10), ],
                         scale = T, 
                         title = pal_name,
                         left.label.size = 0.5,
                         bottom.label.size = 0.15,
                         bottom.label.text.angle = 90,
                         bottom.label.text.alignment = "right",
                         bottom.label.col = "white",
                         row.dendrogram = TRUE,
                         col.dendrogram = TRUE,
                         grid.vline.col = "white",
                         grid.hline.col = "white",
                         heat.pal = my_palette)
  })

colorspace::hcl_palettes("diverging", n = 7, plot = F) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var = "palette") %>%
    apply(., 1, function (x) {
        n <- 7
        pal_name <- as.character(x["palette"])
        my_palette <- colorRampPalette(grDevices::hcl.colors(n, pal_name))(50)
        pheatmap::pheatmap(test, color = my_palette, main = pal_name)
        superheat::superheat(mtcars[sample(1:nrow(mtcars), 10), ],
                             scale = T, 
                             title = pal_name,
                             left.label.size = 0.5,
                             bottom.label.size = 0.15,
                             bottom.label.text.angle = 90,
                             bottom.label.text.alignment = "right",
                             bottom.label.col = "white",
                             row.dendrogram = TRUE,
                             col.dendrogram = TRUE,
                             grid.vline.col = "white",
                             grid.hline.col = "white",
                             heat.pal = my_palette)
    })
    


RColorBrewer::brewer.pal.info %>% 
  .[.$category == "seq", ] %>%
  tibble::rownames_to_column(var = "palette") %>%
  apply(., 1, function (x) {
    n <- as.numeric(x["maxcolors"]) - 2
    pal_name <- as.character(x["palette"])
    my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n, pal_name))(50)
    pheatmap::pheatmap(test, color = my_palette, main = pal_name)
    superheat::superheat(mtcars[sample(1:nrow(mtcars), 10), ],
                         scale = T, 
                         title = pal_name,
                         left.label.size = 0.5,
                         bottom.label.size = 0.15,
                         bottom.label.text.angle = 90,
                         bottom.label.text.alignment = "right",
                         bottom.label.col = "white",
                         row.dendrogram = TRUE,
                         col.dendrogram = TRUE,
                         grid.vline.col = "white",
                         grid.hline.col = "white",
                         heat.pal = my_palette)
  })




RColorBrewer::brewer.pal.info %>% 
  .[.$category == "qual", ] %>%
  tibble::rownames_to_column(var = "palette") %>%
  apply(., 1, function (x) {
    n <- as.numeric(x["maxcolors"]) 
    pal_name <- as.character(x["palette"])
    qual_palette <- RColorBrewer::brewer.pal(n, pal_name)
    ggpubr::ggscatterhist(
      iris, x = "Sepal.Length", y = "Sepal.Width",
      color = "Species", title = pal_name, size = 3, alpha = 0.6,
      palette = qual_palette,   #c("#00AFBB", "#E7B800", "#FC4E07"),
      margin.params = list(fill = "Species", color = "black", size = 0.2)
    )
  })




colorspace::choose_color()



# from jpr #1752A4 (blue)
my_palette <- c('#1952A5', '#8E3C0A', '#B4705B', '#128CB1', '#73A7C3', 
                '#BE9A78', '#D0B49B', '#98BED1', '#C0D5E0', '#E1CEC5')[c(1,4,5,8,9,10,7,6,3,2)]
c("#1752A4", '#1952A5') %>% scales::show_col()


