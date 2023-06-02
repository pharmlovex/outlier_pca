# Define a functions -------------------------------------------------------
# variance compute function 
rowVars <- function(
    x
){
  apply(x, 1, var, na.rm=T)
}


  df = read.delim("GSE229904_mRNA_counts.tsv", row.names = 1)
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