net_top <- function(infile) {
  #
  # Script for creating network topology PCA plots, biplot, and dendrogram from PCLRC result.
  # 
  # Input:
  #       infile: rds file of the result of PCLRC per network
  #
  # Output:
  #       top_pca: ggplot, topology 2D PCA plot
  #       bi_plot: plot, biplot of topology measures
  #       s3d: plot, 3D PCA plot of topology
  #       dendr: plot, dendrogram of networks
  #
  
  source('TopDendro.R')
  install_github("kassambara/factoextra")
  install_github('vqv/ggbiplot')
  
  
  require(stats)
  require(viridis)
  require(scatterplot3d)
  require(tidyverse)
  require(ggsci)
  require(ggbiplot)
  require(devtools)
  require(factoextra)
  require(FactoMineR)
  
  
  # Functions
  TopCalc <- function(infile) {
    
    #  Function for reading RDS PCLRC result, creating an igraph of filtered median correlation matrix and 
    #  calculating topologies
    #
    #  infile: rds file
    #  net_top: dataframe containing network topology
    
    
    source('CalcNetTopology.R')
    require(igraph)
    require(CINNA)
    
    # read RDS file, containing PCLRC_res
    PCLRC_res <- readRDS(infile)
    
    # Select filtered median correlation matrix
    corrmat = PCLRC_res$MedianFilteredCorr
    
    # Create igraph object
    net = graph_from_adjacency_matrix(corrmat, mode = "undirected",  weighted = T, diag = F)
    
    # calculate topologies 
    net_top = CalcNetTopology(net, corrmat)
    
    return(net_top)
    
  }
  
  
  matchbind <- function(...){
    
    # Function for merging multiple dataframes by matching rownames.
    #
    # input: dataframes
    # output: merged dataframe, rows equalized and matched
    
    Reduce(function(x,y){
      cbind(x,y[match(row.names(x),row.names(y)),])}, list(...))
  }
  
  
  # Calculate topologies for all networks
  ctrl_top = TopCalc('control_net_res.rds')
  thyother_top = TopCalc('THYother_net_res.rds')
  thy_top = TopCalc('THY_net_res.rds')
  pso_top = TopCalc('PSO_net_res.rds')
  other_top = TopCalc('other_net_res.rds')
  ms_top = TopCalc('MS_net_res.rds')
  ibdother_top = TopCalc('IBDother_net_res.rds')
  ibd_top = TopCalc('IBD_net_res.rds')
  
  
  # Merge all topology networks and match rownames
  all_top = matchbind(ctrl_top, thyother_top, thy_top, pso_top, other_top, ms_top, ibdother_top, ibd_top)
  # Make columnnames pretty
  colnames(all_top) <- c('ctrl', 'THYother', 'THY', 'PSO', 'other', 'MS', 'IBDother', 'IBD')
  # Remove rows with constant values
  all_top = all_top[apply(all_top, 1, var, na.rm=TRUE) != 0,]
  # Transpose
  all_top = data.frame(t(all_top))
  
  # Run PCA
  top.pca <- prcomp(all_top, scale = T, center = T)
  # Bind data together with first 2 PCs
  pc_data = cbind(all_top, top.pca$x[,1:2])
  # Add type column
  pc_data$Type <- rownames(pc_data)
  
  
  # PC percentage calculations for x and y axis
  percentage <- round(top.pca$sdev / sum(top.pca$sdev) * 100, 2)
  percentage <- paste(colnames(top.pca$x), "(", paste( as.character(percentage), "%", ")", sep="") )
  
  cols = c('other'='#0DCEFF', 'IBDother'='#006E00', 'ctrl'='#DEAE65', 'THYother'='#FF32B8', 'MS'='#009E84','THY'='#AA88EF', 'IBD'='#9B9C0D', 'PSO'='#5883e8')
  
  
  # Plot PCA
  top_pca <- ggplot(pc_data,aes(PC1, PC2, color= Type, label = Type)) + 
    geom_point(size=3) + geom_text(aes(label = Type), hjust=0.5, vjust=1.5, check_overlap = T, size = 4) + 
    xlab(percentage[1]) + ylab(percentage[2]) +  theme(legend.position="none") +
    scale_color_manual(values = cols) +
    ggtitle('PCA plot of network topologies')
  
  
  ############# BIPLOT ##############
  
  biplot(top.pca)
  
  fviz_pca_biplot(top.pca, repel = T, col.var = 'contrib', habillage = 'none')
  
  # factominer to get contributions of variables
  res <- PCA(all_top)
  var <- get_pca_var(res)
  
  # assess variable contributions to PC, select top 15
  most.contrib = fviz_contrib(res, choice = 'var', axes = 1:2, top = 15)
  
  # correlation plot of most important / most contributing variables 
  fviz_pca_var(res, col.var = "contrib",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
  
  
  fviz_pca_biplot(res, repel = T, select.var = list(contrib = 10), col.var = 'contrib')
  
  
  bi_plot = ggbiplot(pcobj = top.pca, choices = c(1,2), obs.scale = 1, var.scale = 1)
  
  
  ######## 3D TOPOLOGY PCA ##########
  
  # bind first 3 PCs to topology df
  pca_3d = cbind(all_top, top.pca$x[,1:3])
  # add type column
  pca_3d$Type <- rownames(pca_3d)
  
  # set colors
  colcodes = c('#ff5c33', '#e6ac00', '#b2b266', '#59b300', '#33cccc', '#00ffff', '#9999ff', '#bf80ff', '#ff66ff')
  # colors <- colcodes[as.factor(pca_3d$Type)]
  # 
  # s3d <- scatterplot3d(pca_3d$PC1, pca_3d$PC2, pca_3d$PC3, pch=20, color = colors, cex.symbols = 1.8,
  #                      main = 'PCA plot of network topologies', 
  #                      xlab = 'PC1', ylab = 'PC2', zlab = 'PC3')
  
  
  colors <- cols[as.factor(pca_3d$Type)]
  
  s3d <- scatterplot3d(pca_3d$PC1, pca_3d$PC2, pca_3d$PC3, pch=20, color = cols[pca_3d$Type], cex.symbols = 2,
                       main = 'PCA plot of network topologies', 
                       xlab = 'PC1', ylab = 'PC2', zlab = 'PC3')
  # add labels
  s3d.coords <- s3d$xyz.convert(pca_3d$PC1, pca_3d$PC2, pca_3d$PC3)
  text(s3d.coords$x, s3d.coords$y, 
       labels=pca_3d$Type,
       pos=4, cex=0.75, pch=25)          
  
  
  ######## TOPOLOGY DENDROGRAM ##########
  dendr <- top_dendr(all_top)
  
  
}
