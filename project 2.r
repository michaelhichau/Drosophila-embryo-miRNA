install.packages("pheatmap", dep=T)
library("pheatmap")

#loading data into data frames
miRNA_count <- read.table("Galaxy-[Join_two_Datasets_on_data_SRR069838_and_data_SRR069840].tabular", col.names = c('miRNA', 'embryo.6_10hr.count', 'miRNA.SRR069840', 'embryo.0_1hr.count'))
row.names(miRNA_count) <- miRNA_count$miRNA
miRNA_count <- miRNA_count[c('embryo.6_10hr.count', 'embryo.0_1hr.count')]




### Normalized data ###
# Note: Every 4 lines in the .fastqsanger file is a single read.
SRR069838_miRNA_total_reads <- length(readLines("Galaxy-[Clip_on_data_SRR069838].fastqsanger"))/4
SRR069840_miRNA_total_reads <- length(readLines("Galaxy-[Clip_on_data_SRR069840].fastqsanger"))/4

# Calculate the normalized read count in RPM
miRNA_count$embryo.6_10hr.normalized <- (miRNA_count$embryo.6_10hr.count)/(SRR069838_miRNA_total_reads/1000000)
miRNA_count$embryo.0_1hr.normalized <- (miRNA_count$embryo.0_1hr.count)/(SRR069840_miRNA_total_reads/1000000)




### Logged Normalized data ###
# remove data where both data sets has a read count lower than 10 RPM
miRNA_count <- miRNA_count[miRNA_count$embryo.6_10hr.normalized > 10 & miRNA_count$embryo.0_1hr.normalized > 10,]

# Logging the normalized read count
miRNA_count$embryo.6_10hr.normalized_log2 <- log2(miRNA_count$embryo.6_10hr.normalized)
#miRNA_count$embryo.6_10hr.normalized_log2[is.infinite(miRNA_count$embryo.6_10hr.normalized_log2)]<- 0 #avoid -inf from log
miRNA_count$embryo.0_1hr.normalized_log2 <- log2(miRNA_count$embryo.0_1hr.normalized)
#miRNA_count$embryo.0_1hr.normalized_log2[is.infinite(miRNA_count$embryo.0_1hr.normalized_log2)]<- 0 #avoid -inf from log




### Fold change and Logged Fold change ###
# Calculate the fold change in RPM and log2(RPM)
miRNA_count$embryo.0_1hr.fold <- miRNA_count$embryo.6_10hr.count/miRNA_count$embryo.0_1hr.count
miRNA_count$embryo.0_1hr.fold_log2 <- log2(miRNA_count$embryo.6_10hr.count/miRNA_count$embryo.0_1hr.count)





### plots ###



### plot 1.1 ###
# sorting order: fold_log2, descending 
sortby_fold_log2 <- miRNA_count[order(- miRNA_count$embryo.0_1hr.fold_log2),]
# data values: normalized_log2
plot1.1 <- sortby_fold_log2[c('embryo.0_1hr.normalized_log2', 'embryo.6_10hr.normalized_log2')]


miRNA_normalized_log2_matrix <- (data.matrix(plot1.1))
miRNA_heatmap <- pheatmap(miRNA_normalized_log2_matrix,
                          cluster_row = FALSE, 
                          cluster_cols = FALSE, 
                          fontsize_row = 5,
                          fontsize_col = 5,
                          border_color = "transparent",
                          angle_col = "0",
                          col = heat.colors(256),
                          scale="none",
)


### plot 1.2 ###
# sorting order: fold_log2, descending 
# data values: normalized_log2, fold_log2 >= 0
plot1.2 <- sortby_fold_log2[sortby_fold_log2$embryo.0_1hr.fold_log2 >= 0,][c('embryo.0_1hr.normalized_log2', 'embryo.6_10hr.normalized_log2')]

plot1.2_matrix <- (data.matrix(plot1.2))
miRNA_heatmap <- pheatmap(plot1.2_matrix, 
                         cluster_row = FALSE, 
                         cluster_cols = FALSE, 
                         fontsize_row = 5,
                         fontsize_col = 5,
                         border_color = "transparent",
                         angle_col = "0",
                         col = heat.colors(256),
                         scale="none",
)


### plot 1.3 ###
# sorting order: fold_log2, descending 
# data values: normalized_log2, fold_log2 <= 0
plot1.3 <- sortby_fold_log2[sortby_fold_log2$embryo.0_1hr.fold_log2 <= 0,][c('embryo.0_1hr.normalized_log2', 'embryo.6_10hr.normalized_log2')]

plot1.3_matrix <- (data.matrix(plot1.3))
miRNA_heatmap <- pheatmap(plot1.3_matrix, 
                         cluster_row = FALSE, 
                         cluster_cols = FALSE, 
                         fontsize_row = 5,
                         fontsize_col = 5,
                         border_color = "transparent",
                         angle_col = "0",
                         col = heat.colors(256),
                         scale="none",
)






### plot 2.1 ###
# sorting order: fold, descending 
sortby_fold <- miRNA_count[order(- miRNA_count$embryo.0_1hr.fold),]
# data values: normalized
plot2.1 <- sortby_fold[c('embryo.6_10hr.normalized', 'embryo.0_1hr.normalized')]


miRNA_normalized_matrix <- (data.matrix(plot2.1))
miRNA_heatmap <- pheatmap(miRNA_normalized_matrix,
                          cluster_row = FALSE, 
                          cluster_cols = FALSE, 
                          fontsize_row = 5,
                          fontsize_col = 3,
                          border_color = "transparent",
                          angle_col = "0",
                          col = heat.colors(256),
                          scale="none",
)


### plot 2.2 ###
# sorting order: fold, descending 
# data values: normalized, fold >= median
normalized_median = median(miRNA_count$embryo.6_10hr.normalized)
plot2.2 <- sortby_fold[sortby_fold$embryo.6_10hr.normalized >= normalized_median,][c('embryo.6_10hr.normalized', 'embryo.0_1hr.normalized')]

plot2.2_matrix <- (data.matrix(plot2.2))
miRNA_heatmap <- pheatmap(plot2.2_matrix, 
                          cluster_row = FALSE, 
                          cluster_cols = FALSE, 
                          fontsize_row = 5,
                          fontsize_col = 3,
                          border_color = "transparent",
                          angle_col = "0",
                          col = heat.colors(256),
                          scale="none",
)


### plot 2.3 ###
# sorting order: fold, descending 
# data values: normalized, fold <= 0
plot2.3 <- sortby_fold[sortby_fold$embryo.6_10hr.normalized <= normalized_median,][c('embryo.6_10hr.normalized', 'embryo.0_1hr.normalized')]

plot2.3_matrix <- (data.matrix(plot2.3))
miRNA_heatmap <- pheatmap(plot2.3_matrix, 
                          cluster_row = FALSE, 
                          cluster_cols = FALSE, 
                          fontsize_row = 5,
                          fontsize_col = 3,
                          border_color = "transparent",
                          angle_col = "0",
                          col = heat.colors(256),
                          scale="none",
)







#loading CHD.xlsx
#into data frames
#install.packages("xlsx") # First install the xlsx R package
#library(xlsx) # load the xlsx package
#CHD_xlsx <- read.xlsx("CHD.xlsx", sheetIndex = 1) # read in the data from sheet number 1

# Exercise 7
check_or_compare_10 <- function(x){
  if (x > 10) 
    print('Our variable is larger than 10')
  else if (x < 10)
    print('Our variable is smaller than 10')
  else if (x == 10) 
    print('Our variable is equal to 10')
}

our_var <- 9
check_or_compare_10(our_var)
our_var <- 10
check_or_compare_10(our_var)
our_var <- 11
check_or_compare_10(our_var)