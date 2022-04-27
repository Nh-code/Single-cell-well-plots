
###################################################################
# File Name: xx.R
# Author: zhuzhiyong
# mail: zhuzhiyong@genomics.cn
#create time   : Wed 27 Apr 2022 11:32:13 PM CST
#last modified : Wed 27 Apr 2022 11:32:13 PM CST
#Update Count    : 0
#=============================================================

clutstDot <- function(object,features){
    library(dplyr);library(ggtree);library(cowplot)
    df <- DotPlot(object,features = features)$data
    mat <- df %>% dplyr::select(-pct.exp, -avg.exp.scaled) %>%  # drop unused columns to faciliate widening
        pivot_wider(names_from = features.plot, values_from = avg.exp ) %>% 
        data.frame() # make df as tibbles -> matrix annoying
    row.names(mat) <- mat$id  # put gene in `row`
    mat <- mat[,-1] #drop gene column as now in rows
    clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
    ddgram <- as.dendrogram(clust) # create dendrogram
    ##################################
    ggtree_plot <- ggtree::ggtree(ddgram,branch.length="none")+
        geom_tippoint(color="#FDAC4F", shape=8, size=3)
    df$id <-  factor(df$id , levels = clust$labels[clust$order])
    dotplot <- ggplot(df,aes(x=features.plot,y = id,size = pct.exp, color = avg.exp.scaled))+
        geom_point() + 
        scale_size("% detected", range = c(0,6)) +
        scale_y_discrete(position = "right") +
        scale_color_gradientn(colours = viridis::viridis(20),
                              guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                              name = "Average\nexpression") +
        cowplot::theme_cowplot() + 
        ylab("") + xlab("Markers") + theme_bw() +
        theme(
            axis.text.x = element_text(size=10, angle=0, hjust=0.5, color="black"),
            axis.text.y = element_text(size=12, color="black"),
            axis.title = element_text(size=14)
        )
    ggtree_plot_yset <- ggtree_plot + aplot::ylim2(dotplot)
    p <- plot_grid(ggtree_plot_yset,NULL,dotplot, nrow = 1, 
                   rel_widths = c(0.5,-0.025, 2), #rel_widths:合并绘图的三张图的相对大小
                   align = 'h')
    return(p)
}