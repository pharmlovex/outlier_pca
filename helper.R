
# Load Library  -----------------------------------------------------------
library(AnnotationHub)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ggfortify)
library(plotly)
library(tidyr)

# Define a functions -------------------------------------------------------
# variance compute function 
rowVars <- function(
    x
){
  apply(x, 1, var, na.rm=T)
}

# Application function 
wrangle <- function(data_path) {
  # Read the data into R env
  df = read.delim(data_path, row.names = 1)
  # Covert the row name to gene symbol
  df$genes = mapIds(org.Hs.eg.db,
                    keys = row.names(df),
                    keytype = 'ENSEMBL',
                    column = 'SYMBOL')
  # Remove rows with no gene symbol
  df <- df %>% 
    drop_na(genes)
  # Remove duplicates from gene symbols
  df <- df[!duplicated(df$genes),]
  # Assign gene name as row names 
  row.names(df) <- df$genes
  # Remove the gene name column
  df <- subset(df, select= -genes)
  # Call the variance compute function
  df$var_exp = rowVars(df)
  df <- filter(df, var_exp > quantile(df$var_exp, 0.9))
  # Remove the mean columns 
  df <- subset(df, select= -var_exp)
  # Transpose the dataframe
  t(df)-> df2
  df2 = as.data.frame(df2)
  
  # Data scaling 
  # Bio-marker density plot before scaling 
  db <- ggplot(data = df2, aes(x=CCNL2)) + 
    geom_density(alpha= 0.5, fill = "#21c661") +
    theme_bw() +
    labs(title = 'Before Scaling the data') +
    theme(plot.title = element_text(hjust = 0.5,))
    
  
  # scale data
  df2_scaled = scale(df2)
  df2_scaled = as.data.frame(df2_scaled)
  # Bio-marker density plot after scaling
  da <- ggplot(data = df2_scaled, aes(x=CCNL2)) + 
    geom_density(alpha= 0.5, fill = "#21c661") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(title = 'After Scaling the data')
  # Perform principal components analysis on the scaled data. 
  
  res <- prcomp(df2_scaled) 
  res$rotation <- -1*res$rotation
  df_pcomp <- res$rotation
  
  # Variance explained by each pc 
  
  var_exp_df = data.frame(PC = paste0("PC",1:175),
                          Var_explained = res$sdev^2/sum(res$sdev^2))
  
  # Create a subset data of the first 10 PCs
  var_df = var_exp_df[1:10,]
  # Factorise the PC column
  var_df$PC = factor(var_df$PC,
                     levels = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7",
                                "PC8","PC9","PC10")
  )
  
  # Visualize the Variance explained by the first 10 PCs
  var_df %>% 
    ggplot(aes(x=PC, y= Var_explained *100, group = 1)) +
    geom_col(fill = "#00cae7")+
    geom_line()+
    labs(x = "Principal Components",
         y = "Variance Explained [%]",
         title = "Variance Explained by first 10 Principal Components")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5)) -> plot_var
  #Weight of bio-marker by PCs
  gene_wt = res$rotation[1:20,1:6]
  
  # First 20 PCs of the 10 samples 
  sample_pc = res$x[1:20,1:6] * -1
  #sample_pc = t(sample_pc) * -1
  
  # Create a subset data for plot
  df_sample = res$x[1:80,1:10]
  
  # PCA Plot
  ggplot(data = df_sample, aes(x=PC1, y=PC2))+
    geom_point(aes(color = "orange"), show.legend = FALSE)+
    geom_text_repel(aes(label= row.names(df_sample))) +
    labs(title = "Principal Components Plots") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))-> pca_plot
  
  # Joint probability density
  ggplot(data = df_sample, aes(x = PC1, y = PC2))+
    geom_density2d(linewidth = 1,
                   show.legend = FALSE)+
    geom_point(aes(color = "red"), show.legend = FALSE)+
    geom_text_repel(aes(label= row.names(df_sample)))+
    labs(title = "Probable Outlier Detection",
         y = "PC2 [10.9% variance explained]",
         x = "PC1 [51.2% variance explained]") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))  -> coutour_plot
  
  
  
  out <- list(db, da, plot_var, gene_wt, sample_pc, pca_plot, coutour_plot)
  names(out) <- c("dp_before_scale", "dp_after_scale", 
                  "variance_explained_plot", "gene_pc_df", "sample_pc_df", 
                  "pca_plot", "joint_prob_density_plot")
  return(out)
}






