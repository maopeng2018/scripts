source('ttest_normal.r');
options(warn=-1)
data <- read.table('data/F-CA_over_F-NC.txt',sep = '\t',header = TRUE);
rownames(data) <- data[, 1];
data <- data[, -1];
results <- TTest(data,2,2);  #### the 2 means two replicates for experiment and control samples

write.csv(results, file = 'Ttest_normal/F-CA_over_F-NC.csv', col.names = T, row.names = T);



