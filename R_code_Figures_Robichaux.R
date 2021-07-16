###############################################################################################
###############################################################################################
# This Script contains the source code for R for figures of the manuscript "Structure-based   #
# classification predicts drug response in EGFR mutant lung cancer" by Robichaux J et al.     #
# published in Nature 2021.                                                                   #
#                                                                                             #
# The code here provided was used for the creation of the following panels using the provided #
# source data that was published with the manuscript.                                         #
# - Figure 2A                                                                                 #
# - Figure 3E                                                                                 #
# - Supplementary Figure 2A                                                                   #
# - Supplementary Figure 3A-C                                                                 #
# - Supplementary Figure 5D                                                                   #
# - Supplementary Figure 6A                                                                   #
#                                                                                             #
# The Figures might have been adjusted afterwards for publishing using external software      #
# however without altering the representation of the data.                                    #
###############################################################################################
###############################################################################################

#load libraries used for the creation of the figures
library(tidyverse)
library(readxl)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(psych)

#####################
##### Figure 2A #####
#####################

#load data
data.2A <- read_excel("2021.6.3 Figure 2 Source Data.xlsx", sheet = "Panel A")

#calculate mean per treatment and mutation
data.2A.melt <- melt(data.2A, id.vars = c("Exon 1", "Exon 2", "Exon 3", "Group", "EGFR Mutation"))
data.2A.melt$variable <- gsub("\\..*", "", data.2A.melt$variable)
data.2A.mean <- data.2A.melt %>% group_by(`EGFR Mutation`, variable) %>% summarize(mean = mean(value, na.rm = TRUE))

data.2A.spread <- dcast(data.2A.mean, `EGFR Mutation` ~ variable)
rownames(data.2A.spread) <- data.2A.spread$`EGFR Mutation`
data.2A.spread <- data.2A.spread[, -1]

#create Heatmap

#prepare color function
col_fun_IC50 = colorRamp2(c(min(data.2A.spread, na.rm = TRUE), -0.5, 1), c("dodgerblue3", "white", "firebrick3"))

#add annotation heatmap with drug clusters
anno.drug <- cbind.data.frame(drug = colnames(data.2A.spread))
anno.drug$DrugClass <- ifelse(anno.drug$drug %in% c("Erlotinib", "Gefitinib", "AZD3759", "Sapatinib"), "Non covalent",
                              ifelse(anno.drug$drug %in% c("Afatinib", "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI"), "2nd Gen",
                                     ifelse(anno.drug$drug %in% c("Osimertinib", "Nazartinib", "Olmutinib", "Rociletinib", "Naquotinib", "Lazertinib"), "3rd Gen",
                                            ifelse(anno.drug$drug %in% c("CLN-081", "AZ5104", "Mobocertinib"), "Ex20ins Specific", 
                                                   ifelse(anno.drug$drug %in% c("CUDC-101"), "EGFR/HDAC",
                                                          ifelse(anno.drug$drug %in% c("Brigatinib", "AZD3463"), "ALK",
                                                                 ifelse(anno.drug$drug %in% c("Ruboxistaurin", "Midostaurin", "Sotrastaurin"), "PKC", NA)))))))


anno.drug.mat <- as.data.frame(anno.drug[, -1])
rownames(anno.drug.mat) <- anno.drug$drug
colnames(anno.drug.mat) <- "Drug Class"
anno.drug.mat <- as.matrix(anno.drug.mat)

#create data for EGFR groups
anno.groups <- cbind.data.frame(Mutation = rownames(data.2A.spread)) %>% left_join(select(data.2A, Group, `EGFR Mutation`), 
                                                                                  by = c("Mutation" = "EGFR Mutation"))

#add annotation for the exon
anno.exon <- cbind.data.frame(Mutation = rownames(data.2A.spread)) %>% left_join(select(data.2A, `EGFR Mutation`, `Exon 1`, `Exon 2`, `Exon 3`), 
                                                                                 by = c("Mutation" = "EGFR Mutation"))
anno.exon$`Exon 1` <- as.character(anno.exon$`Exon 1`)
anno.exon$`Exon 2` <- as.character(anno.exon$`Exon 2`)
anno.exon$`Exon 3` <- as.character(anno.exon$`Exon 3`)

#create heatmap annotation
hm.a <- HeatmapAnnotation(Group = anno.groups$Group, 
                          "Mutation 1" = anno.exon$`Exon 1`,
                          "Mutation 2" = anno.exon$`Exon 2`,
                          "Mutation 3" = anno.exon$`Exon 3`, col = list(Group = c("Classical-Like" = "darkorchid1", "PACC" = "royalblue1", 
                                                                           "T790M-like-3S" = "palegreen2", "Ex20ins-L" = "tomato2", "T790M-like-3R" = "mediumseagreen"),
                                                                 `Mutation 1` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4"),
                                                                 `Mutation 2` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4"),
                                                                 `Mutation 3` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4")), 
                          annotation_height = unit(c(6, 2, 2, 2), "mm"), gap = unit(c(1, 0, 0, 1), "mm"))

#create heatmap with data
hm <- Heatmap(t(data.2A.spread), column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8), col = col_fun_IC50, name = "Log (Mut/WT)", show_row_names = TRUE,
              row_order = c("Erlotinib", "Gefitinib", "AZD3759", "Sapatinib", "Afatinib", "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI",
                            "Osimertinib", "Nazartinib", "Olmutinib", "Rociletinib", "Naquotinib", "Lazertinib", "CLN-081", "AZ5104", "Mobocertinib"), top_annotation = hm.a, 
              clustering_distance_columns = "euclidean")

#create heatmap with drug classes for annotation
hm.drug <- Heatmap(anno.drug.mat, name = "Drug Class", show_row_names = TRUE, col = c("Non covalent" = "darkorchid4", "2nd Gen" = "blue3", "3rd Gen" = "green3",
                                                                                      "Ex20ins Specific" = "red3"))

#draw heatmap
draw(hm + hm.drug)

###############################################################################################


#####################
##### Figure 3E #####
#####################

#load data
data.3E <- read_excel("2021.6.3 Figure 3 Source Data.xlsx", sheet = "Panel E")

#calculate mean per treatment and mutation
data.3E.melt <- melt(data.3E, id.vars = c("Exon 1", "Exon 2", "Group", "EGFR Mutation"))
data.3E.melt$variable <- gsub("\\..*", "", data.3E.melt$variable)
data.3E.mean <- data.3E.melt %>% group_by(`EGFR Mutation`, variable) %>% summarize(mean = mean(value, na.rm = TRUE))

data.3E.spread <- dcast(data.3E.mean, `EGFR Mutation` ~ variable)
rownames(data.3E.spread) <- data.3E.spread$`EGFR Mutation`
data.3E.spread <- data.3E.spread[, -1]



#prepare color function
col_fun_IC50 = colorRamp2(c(min(data.3E.spread, na.rm = TRUE), -0.5, 1), c("dodgerblue3", "white", "firebrick3"))

#add annotation heatmap with drug clusters
anno.drug <- cbind.data.frame(drug = colnames(data.3E.spread))
anno.drug$DrugClass <- ifelse(anno.drug$drug %in% c("Erlotinib", "Gefitinib", "AZD3759", "Sapatinib"), "Non covalent",
                              ifelse(anno.drug$drug %in% c("Afatinib", "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI"), "2nd Gen",
                                     ifelse(anno.drug$drug %in% c("Osimertinib", "Nazartinib", "Olmutinib", "Rociletinib", "Naquotinib", "Lazertinib"), "3rd Gen",
                                            ifelse(anno.drug$drug %in% c("CLN-081", "AZ5104", "Mobocertinib"), "Ex20ins Specific", 
                                                   ifelse(anno.drug$drug %in% c("CUDC-101"), "EGFR/HDAC",
                                                          ifelse(anno.drug$drug %in% c("Brigatinib", "AZD3463"), "ALK",
                                                                 ifelse(anno.drug$drug %in% c("Ruboxistaurin", "Midostaurin", "Sotrastaurin"), "PKC", NA)))))))


anno.drug.mat <- as.data.frame(anno.drug[, -1])
rownames(anno.drug.mat) <- anno.drug$drug
colnames(anno.drug.mat) <- "Drug Class"
anno.drug.mat <- as.matrix(anno.drug.mat)

#create data for EGFR groups
anno.groups <- cbind.data.frame(Mutation = rownames(data.3E.spread)) %>% left_join(select(data.3E, Group, `EGFR Mutation`), 
                                                                                   by = c("Mutation" = "EGFR Mutation"))

#add annotation for the exon
anno.exon <- cbind.data.frame(Mutation = rownames(data.3E.spread)) %>% left_join(select(data.3E, `EGFR Mutation`, `Exon 1`, `Exon 2`), 
                                                                                 by = c("Mutation" = "EGFR Mutation"))
anno.exon$`Exon 1` <- as.character(anno.exon$`Exon 1`)
anno.exon$`Exon 2` <- as.character(anno.exon$`Exon 2`)

#create heatmap annotation
hm.a <- HeatmapAnnotation(Group = anno.groups$Group, 
                          "Mutation 1" = anno.exon$`Exon 1`,
                          "Mutation 2" = anno.exon$`Exon 2`, col = list(Group = c("Classical-Like" = "darkorchid1", "PACC" = "royalblue1"),
                                                                        `Mutation 1` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4"),
                                                                        `Mutation 2` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4")), 
                          annotation_height = unit(c(6, 2, 2), "mm"), gap = unit(c(1, 0, 0, 1), "mm"))

#create heatmap with data
hm <- Heatmap(t(data.3E.spread), column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8), col = col_fun_IC50, name = "Log (Mut/WT)", show_row_names = TRUE,
              row_order = c("Erlotinib", "Gefitinib", "AZD3759", "Sapatinib", "Afatinib", "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI",
                            "Osimertinib", "Nazartinib", "Olmutinib", "Rociletinib", "Naquotinib", "Lazertinib", "CLN-081", "AZ5104", "Mobocertinib"), top_annotation = hm.a, 
              clustering_distance_columns = "euclidean")

#create heatmap with drug classes for annotation
hm.drug <- Heatmap(anno.drug.mat, name = "Drug Class", show_row_names = TRUE, col = c("Non covalent" = "darkorchid4", "2nd Gen" = "blue3", "3rd Gen" = "green3",
                                                                                      "Ex20ins Specific" = "red3"))

#draw heatmap
draw(hm + hm.drug)


###############################################################################################


###################################
##### Supplementary Figure 2 #####
###################################

#load data
data.S2 <- read_excel("2021.6.3 Figure 2 Source Data.xlsx", sheet = "Panel A")

#calculate mean per treatment and mutation
data.S2.melt <- melt(data.S2, id.vars = c("Exon 1", "Exon 2", "Exon 3", "Group", "EGFR Mutation"))
data.S2.melt$variable <- gsub("\\..*", "", data.S2.melt$variable)
data.S2.mean <- data.S2.melt %>% group_by(`EGFR Mutation`, variable) %>% summarize(mean = mean(value, na.rm = TRUE))

data.S2.spread <- dcast(data.S2.mean, `EGFR Mutation` ~ variable)
rownames(data.S2.spread) <- data.S2.spread$`EGFR Mutation`
data.S2.spread <- data.S2.spread[, -1]

#create data for EGFR groups
anno.groups <- cbind.data.frame(Mutation = rownames(data.S2.spread)) %>% left_join(select(data.S2, Group, `EGFR Mutation`), 
                                                                                   by = c("Mutation" = "EGFR Mutation"))

#add annotation for the exon
anno.exon <- cbind.data.frame(Mutation = rownames(data.S2.spread)) %>% left_join(select(data.S2, `EGFR Mutation`, `Exon 1`, `Exon 2`, `Exon 3`), 
                                                                                 by = c("Mutation" = "EGFR Mutation"))
anno.exon$`Exon 1` <- as.character(anno.exon$`Exon 1`)
anno.exon$`Exon 2` <- as.character(anno.exon$`Exon 2`)
anno.exon$`Exon 3` <- as.character(anno.exon$`Exon 3`)


#create Heatmap

#prepare color function
col_fun_IC50 = colorRamp2(c(min(data.S2.spread, na.rm = TRUE), -0.5, 1), c("dodgerblue3", "white", "firebrick3"))

#add annotation heatmap with drug clusters
anno.drug <- cbind.data.frame(drug = colnames(data.S2.spread))
anno.drug$DrugClass <- ifelse(anno.drug$drug %in% c("Erlotinib", "Gefitinib", "AZD3759", "Sapatinib"), "Non covalent",
                              ifelse(anno.drug$drug %in% c("Afatinib", "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI"), "2nd Gen",
                                     ifelse(anno.drug$drug %in% c("Osimertinib", "Nazartinib", "Olmutinib", "Rociletinib", "Naquotinib", "Lazertinib"), "3rd Gen",
                                            ifelse(anno.drug$drug %in% c("CLN-081", "AZ5104", "Mobocertinib"), "Ex20ins Specific", 
                                                   ifelse(anno.drug$drug %in% c("CUDC-101"), "EGFR/HDAC",
                                                          ifelse(anno.drug$drug %in% c("Brigatinib", "AZD3463"), "ALK",
                                                                 ifelse(anno.drug$drug %in% c("Ruboxistaurin", "Midostaurin", "Sotrastaurin"), "PKC", NA)))))))


anno.drug.mat <- as.data.frame(anno.drug[, -1])
rownames(anno.drug.mat) <- anno.drug$drug
colnames(anno.drug.mat) <- "Drug Class"
anno.drug.mat <- as.matrix(anno.drug.mat)

#order data by Exon 1
data.S2.spread.ordered <- data.S2.spread %>% rownames_to_column(var = "Mutation") %>% left_join(anno.exon[, 1:4]) %>% arrange(`Exon 1`)
rownames(data.S2.spread.ordered) <- data.S2.spread.ordered$Mutation
data.S2.spread.ordered.m <- as.matrix(data.S2.spread.ordered[, -c(1, 20:22)])

anno.groups <- cbind.data.frame(Mutation = data.S2.spread.ordered$Mutation) %>% left_join(select(data.S2, Group, `EGFR Mutation`), by = c("Mutation" = "EGFR Mutation"))
anno.exon <- cbind.data.frame(Mutation = data.S2.spread.ordered$Mutation) %>% left_join(select(data.S2, `EGFR Mutation`, `Exon 1`, `Exon 2`, `Exon 3`), by = c("Mutation" = "EGFR Mutation"))
anno.exon$`Exon 1` <- as.character(anno.exon$`Exon 1`)
anno.exon$`Exon 2` <- as.character(anno.exon$`Exon 2`)
anno.exon$`Exon 3` <- as.character(anno.exon$`Exon 3`)





#create heatmap annotation
hm.a <- HeatmapAnnotation(Group = anno.groups$Group, 
                          "Mutation 1" = anno.exon$`Exon 1`,
                          "Mutation 2" = anno.exon$`Exon 2`,
                          "Mutation 3" = anno.exon$`Exon 3`, col = list(Group = c("Classical-Like" = "darkorchid1", "PACC" = "royalblue1", 
                                                                                  "T790M-like-3S" = "palegreen2", "Ex20ins-L" = "tomato2", "T790M-like-3R" = "mediumseagreen"),
                                                                        `Mutation 1` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4"),
                                                                        `Mutation 2` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4"),
                                                                        `Mutation 3` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4")), 
                          annotation_height = unit(c(6, 2, 2, 2), "mm"), gap = unit(c(1, 0, 0, 1), "mm"))

#create heatmap with data
hm <- Heatmap(t(data.S2.spread.ordered.m), column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8), col = col_fun_IC50, name = "Log (Mut/WT)", show_row_names = TRUE,
              row_order = c("Erlotinib", "Gefitinib", "AZD3759", "Sapatinib", "Afatinib", "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI",
                            "Osimertinib", "Nazartinib", "Olmutinib", "Rociletinib", "Naquotinib", "Lazertinib", "CLN-081", "AZ5104", "Mobocertinib"), top_annotation = hm.a, 
              cluster_columns = FALSE)

#create heatmap with drug classes for annotation
hm.drug <- Heatmap(anno.drug.mat, name = "Drug Class", show_row_names = TRUE, col = c("Non covalent" = "darkorchid4", "2nd Gen" = "blue3", "3rd Gen" = "green3",
                                                                                      "Ex20ins Specific" = "red3"))

#draw heatmap
draw(hm + hm.drug)


###################################
##### Supplementary Figure 2B #####
###################################

#the data is teh same for supplementary Figure 2A and 2B with the single difference that the 2B panel is ordered by structure-function
#groups while the 2A panel is ordered by exons. Thus the code is largely the same and only the changes are highlighted here


#order data by Group
data.S2.spread.ordered <- data.S2.spread %>% rownames_to_column(var = "Mutation") %>% left_join(anno.groups) %>% arrange(Group)
rownames(data.S2.spread.ordered) <- data.S2.spread.ordered$Mutation
data.S2.spread.ordered.m <- as.matrix(data.S2.spread.ordered[, -c(1, 20:22)])

anno.groups <- cbind.data.frame(Mutation = data.S2.spread.ordered$Mutation) %>% left_join(select(data.S2, Group, `EGFR Mutation`), by = c("Mutation" = "EGFR Mutation"))
anno.exon <- cbind.data.frame(Mutation = data.S2.spread.ordered$Mutation) %>% left_join(select(data.S2, `EGFR Mutation`, `Exon 1`, `Exon 2`, `Exon 3`), by = c("Mutation" = "EGFR Mutation"))
anno.exon$`Exon 1` <- as.character(anno.exon$`Exon 1`)
anno.exon$`Exon 2` <- as.character(anno.exon$`Exon 2`)
anno.exon$`Exon 3` <- as.character(anno.exon$`Exon 3`)





#create heatmap annotation
hm.a <- HeatmapAnnotation(Group = anno.groups$Group, 
                          "Mutation 1" = anno.exon$`Exon 1`,
                          "Mutation 2" = anno.exon$`Exon 2`,
                          "Mutation 3" = anno.exon$`Exon 3`, col = list(Group = c("Classical-Like" = "darkorchid1", "PACC" = "royalblue1", 
                                                                                  "T790M-like-3S" = "palegreen2", "Ex20ins-L" = "tomato2", "T790M-like-3R" = "mediumseagreen"),
                                                                        `Mutation 1` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4"),
                                                                        `Mutation 2` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4"),
                                                                        `Mutation 3` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4")), 
                          annotation_height = unit(c(6, 2, 2, 2), "mm"), gap = unit(c(1, 0, 0, 1), "mm"))

#create heatmap with data
hm <- Heatmap(t(data.S2.spread.ordered.m), column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8), col = col_fun_IC50, name = "Log (Mut/WT)", show_row_names = TRUE,
              row_order = c("Erlotinib", "Gefitinib", "AZD3759", "Sapatinib", "Afatinib", "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI",
                            "Osimertinib", "Nazartinib", "Olmutinib", "Rociletinib", "Naquotinib", "Lazertinib", "CLN-081", "AZ5104", "Mobocertinib"), top_annotation = hm.a, 
              cluster_columns = FALSE)

#create heatmap with drug classes for annotation
hm.drug <- Heatmap(anno.drug.mat, name = "Drug Class", show_row_names = TRUE, col = c("Non covalent" = "darkorchid4", "2nd Gen" = "blue3", "3rd Gen" = "green3",
                                                                                      "Ex20ins Specific" = "red3"))

#draw heatmap
draw(hm + hm.drug)


###############################################################################################


###################################
##### Supplementary Figure 3A #####
###################################

#the results of Figure S3A are based on the data in figure 2A and thus the same dataset is used

#load data
data.2A <- read_excel("2021.6.3 Figure 2 Source Data.xlsx", sheet = "Panel A")

#calculate mean per treatment and mutation
data.2A.melt <- melt(data.2A, id.vars = c("Exon 1", "Exon 2", "Exon 3", "Group", "EGFR Mutation"))
data.2A.melt$variable <- gsub("\\..*", "", data.2A.melt$variable)
data.2A.mean <- data.2A.melt %>% group_by(`EGFR Mutation`, variable) %>% summarize(mean = mean(value, na.rm = TRUE))

data.2A.spread <- dcast(data.2A.mean, `EGFR Mutation` ~ variable)
rownames(data.2A.spread) <- data.2A.spread$`EGFR Mutation`
data.2A.spread <- data.2A.spread[, -1]

#get group information for each data. Melt data frame
data.S3 <- data.2A.spread %>% rownames_to_column(var = "EGFR Mutation") %>% melt()

#add information with Group
data.S3 <- data.S3 %>% left_join(select(data.2A, `EGFR Mutation`, `Exon 1`))

##################################
# Prepare Exon based correlation #
##################################

#for exon 18
data.corr.18 <- filter(data.S3, `Exon 1` == 18)
data.corr.18.spread <- data.corr.18 %>% spread(`EGFR Mutation`, value)

results.exon18 <- c()

for(i in colnames(data.corr.18.spread)[-c(1:2)]){
  data.corr.18.spread$mean <- rowMeans(select(data.corr.18.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.18.spread[, i], data.corr.18.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.18.spread[, i], data.corr.18.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.exon18[[i]] <- res
  
}

results.exon18 <- do.call("rbind", results.exon18)

#for exon 19
data.corr.19 <- filter(data.S3, `Exon 1` == 19)
data.corr.19.spread <- data.corr.19 %>% spread(`EGFR Mutation`, value)


results.exon19 <- c()

for(i in colnames(data.corr.19.spread)[-c(1:2)]){
  data.corr.19.spread$mean <- rowMeans(select(data.corr.19.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.19.spread[, i], data.corr.19.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.19.spread[, i], data.corr.19.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.exon19[[i]] <- res
  
}

results.exon19 <- do.call("rbind", results.exon19)


#for exon 20
data.corr.20 <- filter(data.S3, `Exon 1` == 20)
data.corr.20.spread <- data.corr.20 %>% spread(`EGFR Mutation`, value)


results.exon20 <- c()

for(i in colnames(data.corr.20.spread)[-c(1:2)]){
  data.corr.20.spread$mean <- rowMeans(select(data.corr.20.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.20.spread[, i], data.corr.20.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.20.spread[, i], data.corr.20.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.exon20[[i]] <- res
  
}

results.exon20 <- do.call("rbind", results.exon20)



#for exon 21
data.corr.21 <- filter(data.S3, `Exon 1` == 21)
data.corr.21.spread <- data.corr.21 %>% spread(`EGFR Mutation`, value)


results.exon21 <- c()

for(i in colnames(data.corr.21.spread)[-c(1:2)]){
  data.corr.21.spread$mean <- rowMeans(select(data.corr.21.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.21.spread[, i], data.corr.21.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.21.spread[, i], data.corr.21.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.exon21[[i]] <- res
  
}

results.exon21 <- do.call("rbind", results.exon21)


#combine data
results.exon <- rbind.data.frame(results.exon18,
                                 results.exon19,
                                 results.exon20,
                                 results.exon21)

###################################
# Prepare Group based correlation #
###################################

#add information with Group
data.S3 <- data.2A.spread %>% rownames_to_column(var = "EGFR Mutation") %>% melt()
data.S3 <- data.S3 %>% left_join(select(data.2A, `EGFR Mutation`, Group))


#for classical
data.corr.classical <- filter(data.S3, Group == "Classical-Like")
data.corr.classical.spread <- data.corr.classical %>% spread(`EGFR Mutation`, value)


results.classical <- c()

for(i in colnames(data.corr.classical.spread)[-c(1:2)]){
  data.corr.classical.spread$mean <- rowMeans(select(data.corr.classical.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.classical.spread[, i], data.corr.classical.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.classical.spread[, i], data.corr.classical.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.classical[[i]] <- res
  
}


results.classical <- do.call("rbind", results.classical)


#for PACC
data.corr.PACC <- filter(data.S3, Group == "PACC")
data.corr.PACC.spread <- data.corr.PACC %>% spread(`EGFR Mutation`, value)


results.PACC <- c()

for(i in colnames(data.corr.PACC.spread)[-c(1:2)]){
  data.corr.PACC.spread$mean <- rowMeans(select(data.corr.PACC.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.PACC.spread[, i], data.corr.PACC.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.PACC.spread[, i], data.corr.PACC.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.PACC[[i]] <- res
  
}

results.PACC <- do.call("rbind", results.PACC)

#for Ex20
data.corr.Ex20 <- filter(data.S3, Group == "Ex20ins-L")
data.corr.Ex20.spread <- data.corr.Ex20 %>% spread(`EGFR Mutation`, value)


results.Ex20 <- c()

for(i in colnames(data.corr.Ex20.spread)[-c(1:2)]){
  data.corr.Ex20.spread$mean <- rowMeans(select(data.corr.Ex20.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.Ex20.spread[, i], data.corr.Ex20.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.Ex20.spread[, i], data.corr.Ex20.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.Ex20[[i]] <- res
  
}

results.Ex20 <- do.call("rbind", results.Ex20)

#for T790M-like
data.corr.T790M <- filter(data.S3, Group == "T790M-like-3S")
data.corr.T790M.spread <- data.corr.T790M %>% spread(`EGFR Mutation`, value)


results.T790M <- c()

for(i in colnames(data.corr.T790M.spread)[-c(1:2)]){
  data.corr.T790M.spread$mean <- rowMeans(select(data.corr.T790M.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.T790M.spread[, i], data.corr.T790M.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.T790M.spread[, i], data.corr.T790M.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.T790M[[i]] <- res
  
}

results.T790M <- do.call("rbind", results.T790M)

#for T790M-like-R
data.corr.T790MR <- filter(data.S3, Group == "T790M-like-3R")
data.corr.T790MR.spread <- data.corr.T790MR %>% spread(`EGFR Mutation`, value)


results.T790MR <- c()

for(i in colnames(data.corr.T790MR.spread)[-c(1:2)]){
  data.corr.T790MR.spread$mean <- rowMeans(select(data.corr.T790MR.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.T790MR.spread[, i], data.corr.T790MR.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.T790MR.spread[, i], data.corr.T790MR.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.T790MR[[i]] <- res
  
}

results.T790MR <- do.call("rbind", results.T790MR)



#Combine data
results.group <- rbind.data.frame(results.classical,
                                  results.PACC,
                                  results.Ex20,
                                  results.T790M,
                                  results.T790MR)



#combine Exon-based and group-based
results.combined <- results.exon %>% rownames_to_column(var = "Mutation") %>% inner_join(rownames_to_column(results.group, var = "Mutation"), 
                                                                                         by = c("Mutation" = "Mutation"))

colnames(results.combined) <- c("Mutation", "rho_exon", "p_exon", "rho_groups", "p_groups")


#Create Figure
results.barplot <- results.combined %>% mutate("rho_exon_bar" = -rho_exon, "delta" = rho_groups - rho_exon) %>% arrange(delta)
results.barplot$Mutation <- factor(results.barplot$Mutation, levels = results.barplot$Mutation)

barplot.compare <- ggplot(results.barplot, aes(x = rho_exon_bar, y = Mutation)) + geom_bar(stat = "identity", fill = "palegoldenrod") + 
  geom_bar(stat = "identity", aes(x = rho_groups), fill = "yellowgreen") + geom_bar(aes(x = delta), stat = "identity", fill = "gray80", col = "black", alpha = 0.5) +
  theme_bw() + geom_vline(xintercept = 0) + theme(axis.text.y = element_text(size = 6)) + xlim(-1, 1)

#plot
barplot.compare


###############################################################################################


###################################
##### Supplementary Figure 3B #####
###################################

#the results of Figure S3B are based on the data in Figure 2A and thus the same dataset is used.


#load data
data.2A <- read_excel("2021.6.3 Figure 2 Source Data.xlsx", sheet = "Panel A")

#calculate mean per treatment and mutation
data.2A.melt <- melt(data.2A, id.vars = c("Exon 1", "Exon 2", "Exon 3", "Group", "EGFR Mutation"))
data.2A.melt$variable <- gsub("\\..*", "", data.2A.melt$variable)
data.2A.mean <- data.2A.melt %>% group_by(`EGFR Mutation`, variable) %>% summarize(mean = mean(value, na.rm = TRUE))

data.2A.spread <- dcast(data.2A.mean, `EGFR Mutation` ~ variable)
rownames(data.2A.spread) <- data.2A.spread$`EGFR Mutation`
data.2A.spread <- data.2A.spread[, -1]

#get group information for each data. Melt data frame
data.S3B <- data.2A.spread %>% rownames_to_column(var = "EGFR Mutation") %>% melt()

#add information with Group
data.S3B <- data.S3B %>% left_join(select(data.2A, `EGFR Mutation`, `Exon 1`))

#filter out samples with T790M mutation
data.S3B <- data.S3B[!grepl(".*T790M.*", data.S3B$`EGFR Mutation`), ]

##################################
# Prepare Exon based correlation #
##################################

#for exon 18
data.corr.18 <- filter(data.S3B, `Exon 1` == 18)
data.corr.18.spread <- data.corr.18 %>% spread(`EGFR Mutation`, value)

results.exon18 <- c()

for(i in colnames(data.corr.18.spread)[-c(1:2)]){
  data.corr.18.spread$mean <- rowMeans(select(data.corr.18.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.18.spread[, i], data.corr.18.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.18.spread[, i], data.corr.18.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.exon18[[i]] <- res
  
}

results.exon18 <- do.call("rbind", results.exon18)

#for exon 19
data.corr.19 <- filter(data.S3B, `Exon 1` == 19)
data.corr.19.spread <- data.corr.19 %>% spread(`EGFR Mutation`, value)


results.exon19 <- c()

for(i in colnames(data.corr.19.spread)[-c(1:2)]){
  data.corr.19.spread$mean <- rowMeans(select(data.corr.19.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.19.spread[, i], data.corr.19.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.19.spread[, i], data.corr.19.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.exon19[[i]] <- res
  
}

results.exon19 <- do.call("rbind", results.exon19)


#for exon 20
data.corr.20 <- filter(data.S3B, `Exon 1` == 20)
data.corr.20.spread <- data.corr.20 %>% spread(`EGFR Mutation`, value)


results.exon20 <- c()

for(i in colnames(data.corr.20.spread)[-c(1:2)]){
  data.corr.20.spread$mean <- rowMeans(select(data.corr.20.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.20.spread[, i], data.corr.20.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.20.spread[, i], data.corr.20.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.exon20[[i]] <- res
  
}

results.exon20 <- do.call("rbind", results.exon20)



#for exon 21
data.corr.21 <- filter(data.S3B, `Exon 1` == 21)
data.corr.21.spread <- data.corr.21 %>% spread(`EGFR Mutation`, value)


results.exon21 <- c()

for(i in colnames(data.corr.21.spread)[-c(1:2)]){
  data.corr.21.spread$mean <- rowMeans(select(data.corr.21.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.21.spread[, i], data.corr.21.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.21.spread[, i], data.corr.21.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.exon21[[i]] <- res
  
}

results.exon21 <- do.call("rbind", results.exon21)


#combine data
results.exon <- rbind.data.frame(results.exon18,
                                 results.exon19,
                                 results.exon20,
                                 results.exon21)

###################################
# Prepare Group based correlation #
###################################

#add information with Group
data.S3B <- data.2A.spread %>% rownames_to_column(var = "EGFR Mutation") %>% melt()
data.S3B <- data.S3B %>% left_join(select(data.2A, `EGFR Mutation`, Group))
data.S3B <- data.S3B[!grepl(".*T790M.*", data.S3B$`EGFR Mutation`), ]

#for classical
data.corr.classical <- filter(data.S3B, Group == "Classical-Like")
data.corr.classical.spread <- data.corr.classical %>% spread(`EGFR Mutation`, value)


results.classical <- c()

for(i in colnames(data.corr.classical.spread)[-c(1:2)]){
  data.corr.classical.spread$mean <- rowMeans(select(data.corr.classical.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.classical.spread[, i], data.corr.classical.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.classical.spread[, i], data.corr.classical.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.classical[[i]] <- res
  
}


results.classical <- do.call("rbind", results.classical)


#for PACC
data.corr.PACC <- filter(data.S3B, Group == "PACC")
data.corr.PACC.spread <- data.corr.PACC %>% spread(`EGFR Mutation`, value)


results.PACC <- c()

for(i in colnames(data.corr.PACC.spread)[-c(1:2)]){
  data.corr.PACC.spread$mean <- rowMeans(select(data.corr.PACC.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.PACC.spread[, i], data.corr.PACC.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.PACC.spread[, i], data.corr.PACC.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.PACC[[i]] <- res
  
}

results.PACC <- do.call("rbind", results.PACC)

#for Ex20
data.corr.Ex20 <- filter(data.S3B, Group == "Ex20ins-L")
data.corr.Ex20.spread <- data.corr.Ex20 %>% spread(`EGFR Mutation`, value)


results.Ex20 <- c()

for(i in colnames(data.corr.Ex20.spread)[-c(1:2)]){
  data.corr.Ex20.spread$mean <- rowMeans(select(data.corr.Ex20.spread[, -c(1:2)], !i)) 
  res.r <- corr.test(data.corr.Ex20.spread[, i], data.corr.Ex20.spread$mean, method = "spearman")[["r"]]
  res.p <- corr.test(data.corr.Ex20.spread[, i], data.corr.Ex20.spread$mean, method = "spearman")[["p"]]
  res <- cbind.data.frame("rho" = res.r, "p-value" = res.p)
  results.Ex20[[i]] <- res
  
}

results.Ex20 <- do.call("rbind", results.Ex20)

#Combine data
results.group <- rbind.data.frame(results.classical,
                                  results.PACC,
                                  results.Ex20)



#combine Exon-based and group-based
results.combined <- results.exon %>% rownames_to_column(var = "Mutation") %>% inner_join(rownames_to_column(results.group, var = "Mutation"), 
                                                                                         by = c("Mutation" = "Mutation"))

colnames(results.combined) <- c("Mutation", "rho_exon", "p_exon", "rho_groups", "p_groups")


#Create Figure
results.barplot <- results.combined %>% mutate("rho_exon_bar" = -rho_exon, "delta" = rho_groups - rho_exon) %>% arrange(delta)
results.barplot$Mutation <- factor(results.barplot$Mutation, levels = results.barplot$Mutation)

barplot.compare <- ggplot(results.barplot, aes(x = rho_exon_bar, y = Mutation)) + geom_bar(stat = "identity", fill = "palegoldenrod") + 
  geom_bar(stat = "identity", aes(x = rho_groups), fill = "yellowgreen") + geom_bar(aes(x = delta), stat = "identity", fill = "gray80", col = "black", alpha = 0.5) +
  theme_bw() + geom_vline(xintercept = 0) + theme(axis.text.y = element_text(size = 6)) + xlim(-1, 1)

#plot
barplot.compare


###############################################################################################


###################################
##### Supplementary Figure 5D #####
###################################

data.5 <- read_excel("2021.6.3 Extended Figure 5 Source Data .xlsx", sheet = "Panel D")

#Prepare data
data.5.melt <- melt(data.5, id.vars = "EGFR Mutation")
data.5.melt$variable <- gsub("\\..*", "", data.5.melt$variable)

data.5.mean <- data.5.melt %>% group_by(`EGFR Mutation`, variable) %>% summarize(mean = mean(value, na.rm = TRUE))

#spread
data.5.spread <- dcast(data.5.mean, `EGFR Mutation` ~ variable)
rownames(data.5.spread) <- data.5.spread$`EGFR Mutation`
data.5.spread <- data.5.spread[, -1]

#create Heatmap

#prepare color function
col_fun_IC50 = circlize::colorRamp2(c(-4, -0.5, 1), c("dodgerblue3", "white", "firebrick3"))

#add annotation heatmap with drug clusters
anno.drug <- cbind.data.frame(drug = colnames(data.5.spread))
anno.drug$DrugClass <- ifelse(anno.drug$drug %in% c("Erlotinib", "Gefitinib", "AZD3759", "Sapatinib"), "Non covalent",
                              ifelse(anno.drug$drug %in% c("Afatinib", "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI"), "2nd Gen",
                                     ifelse(anno.drug$drug %in% c("Osimertinib", "Nazartinib", "Olmutinib", "Rocelitinib", "Naquotinib", "Lazartinib"), "3rd Gen",
                                            ifelse(anno.drug$drug %in% c("TAS6417", "AZ5104", "Mobocertinib", "CLN-081"), "Exon20 Spec", 
                                                   ifelse(anno.drug$drug %in% c("CUDC-101"), "EGFR/HDAC",
                                                          ifelse(anno.drug$drug %in% c("Brigatinib", "AZD3463"), "ALK",
                                                                 ifelse(anno.drug$drug %in% c("Ruboxistaurin", "Midostaurin", "Sotrastaurin"), "PKC", NA)))))))


anno.drug.mat <- as.data.frame(anno.drug[, -1])
rownames(anno.drug.mat) <- anno.drug$drug
colnames(anno.drug.mat) <- "Drug Class"
anno.drug.mat <- as.matrix(anno.drug.mat)


anno.groups <- cbind.data.frame(Mutation = rownames(data.5.spread))
anno.groups$Group <- gsub(c("*.773.*", "*.774.*"), "Far",  anno.groups$Mutation)
anno.groups$Group <- ifelse(str_detect(anno.groups$Mutation, c(".*773.*|.*774.*")), "Far", "Near")

#manually check for double mutations adn change if necessary
anno.groups$Group <- ifelse(anno.groups$Mutation == "D770insY/H773Y", "Near", anno.groups$Group)

hm.a <- HeatmapAnnotation(Group = anno.groups$Group, 
                          col = list(Group = c("Near" = "firebrick1", "Far" = "firebrick4")), 
                          annotation_height = unit(6, "mm"))

hm <- Heatmap(t(data.5.spread), column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8), col = col_fun_IC50, name = "Log (Mut/WT)", show_row_names = TRUE,
              row_order = c("Afatinib",  "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI", "CLN-081", "AZ5104", "Mobocertinib"), top_annotation = hm.a, 
              clustering_distance_columns = "euclidean")


hm.drug <- Heatmap(anno.drug.mat, name = "Drug Class", show_row_names = TRUE, col = c("2nd Gen" = "blue3", "Exon20 Spec" = "red3"))

#safe to disk
draw(hm + hm.drug)


###############################################################################################


###################################
##### Supplementary Figure 6A #####
###################################

data.S6 <- read_excel("2021.6.3 Extended Figure 6 Source Data .xlsx", sheet = "Panel A")

#calculate mean per treatment and mutation
data.S6.melt <- melt(data.S6, id.vars = c("Exon 1", "Exon 2", "Exon 3", "Group", "EGFR Mutation"))
data.S6.melt$variable <- gsub("\\..*", "", data.S6.melt$variable)
data.S6.mean <- data.S6.melt %>% group_by(`EGFR Mutation`, variable) %>% summarize(mean = mean(value, na.rm = TRUE))

data.S6.spread <- dcast(data.S6.mean, `EGFR Mutation` ~ variable)
rownames(data.S6.spread) <- data.S6.spread$`EGFR Mutation`
data.S6.spread <- data.S6.spread[, -1]

#create Heatmap

#prepare color function
col_fun_IC50 = colorRamp2(c(min(data.S6.spread, na.rm = TRUE), -0.5, 1), c("dodgerblue3", "white", "firebrick3"))

#add annotation heatmap with drug clusters
anno.drug <- cbind.data.frame(drug = colnames(data.S6.spread))
anno.drug$DrugClass <- ifelse(anno.drug$drug %in% c("Erlotinib", "Gefitinib", "AZD3759", "Sapatinib"), "Non covalent",
                              ifelse(anno.drug$drug %in% c("Afatinib", "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI"), "2nd Gen",
                                     ifelse(anno.drug$drug %in% c("Osimertinib", "Nazartinib", "Olmutinib", "Rociletinib", "Naquotinib", "Lazertinib"), "3rd Gen",
                                            ifelse(anno.drug$drug %in% c("CLN-081", "AZ5104", "Mobocertinib"), "Ex20ins Specific", 
                                                   ifelse(anno.drug$drug %in% c("CUDC-101"), "EGFR/HDAC",
                                                          ifelse(anno.drug$drug %in% c("Brigatinib", "AZD3463"), "ALK",
                                                                 ifelse(anno.drug$drug %in% c("Ruboxistaurin", "Midostaurin", "Sotrastaurin"), "PKC", NA)))))))


anno.drug.mat <- as.data.frame(anno.drug[, -1])
rownames(anno.drug.mat) <- anno.drug$drug
colnames(anno.drug.mat) <- "Drug Class"
anno.drug.mat <- as.matrix(anno.drug.mat)

#create data for EGFR groups
anno.groups <- cbind.data.frame(Mutation = rownames(data.S6.spread)) %>% left_join(select(data.S6, Group, `EGFR Mutation`), 
                                                                                   by = c("Mutation" = "EGFR Mutation"))

#add annotation for the exon
anno.exon <- cbind.data.frame(Mutation = rownames(data.S6.spread)) %>% left_join(select(data.S6, `EGFR Mutation`, `Exon 1`, `Exon 2`, `Exon 3`), 
                                                                                 by = c("Mutation" = "EGFR Mutation"))
anno.exon$`Exon 1` <- as.character(anno.exon$`Exon 1`)
anno.exon$`Exon 2` <- as.character(anno.exon$`Exon 2`)
anno.exon$`Exon 3` <- as.character(anno.exon$`Exon 3`)

#create heatmap annotation
hm.a <- HeatmapAnnotation(Group = anno.groups$Group, 
                          "Mutation 1" = anno.exon$`Exon 1`,
                          "Mutation 2" = anno.exon$`Exon 2`,
                          "Mutation 3" = anno.exon$`Exon 3`, col = list(Group = c("Classical-Like" = "darkorchid1", "PACC" = "royalblue1", 
                                                                                  "T790M-like-3S" = "palegreen2", "Ex20ins-L" = "tomato2", "T790M-like-3R" = "mediumseagreen"),
                                                                        `Mutation 1` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4"),
                                                                        `Mutation 2` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4"),
                                                                        `Mutation 3` = c("18" = "lightgoldenrod1", "19" = "navy", "20" = "orange2", "21" = "burlywood4")), 
                          annotation_height = unit(c(6, 2, 2, 2), "mm"), gap = unit(c(1, 0, 0, 1), "mm"))

#create heatmap with data
hm <- Heatmap(t(data.S6.spread), column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8), col = col_fun_IC50, name = "Log (Mut/WT)", show_row_names = TRUE,
              row_order = c("Erlotinib", "Gefitinib", "AZD3759", "Sapatinib", "Afatinib", "Dacomitinib", "Neratinib", "Poziotinib", "Tarlox-TKI",
                            "Osimertinib", "Nazartinib", "Olmutinib", "Rociletinib", "Naquotinib", "Lazertinib", "CLN-081", "AZ5104", "Mobocertinib",
                            "CUDC-101", "Brigatinib", "AZD3463", "Ruboxistaurin", "Midostaurin", "Sotrastaurin"), top_annotation = hm.a, 
              clustering_distance_columns = "euclidean")

#create heatmap with drug classes for annotation
hm.drug <- Heatmap(anno.drug.mat, name = "Drug Class", show_row_names = TRUE, col = c("Non covalent" = "darkorchid4", "2nd Gen" = "blue3", "3rd Gen" = "green3",
                                                                                      "Ex20ins Specific" = "red3", "EGFR/HDAC" = "grey75", "ALK" = "darkorange3", "PKC" = "cyan3"))

#draw heatmap
draw(hm + hm.drug)

###############################################################################################

#Session Info
R version 4.0.4 (2021-02-15)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 17763)

Matrix products: default

locale:
  [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] psych_2.1.3          circlize_0.4.12      ComplexHeatmap_2.6.2 reshape2_1.4.4       readxl_1.3.1         forcats_0.5.1       
[7] stringr_1.4.0        dplyr_1.0.4          purrr_0.3.4          readr_1.4.0          tidyr_1.1.2          tibble_3.1.0        
[13] ggplot2_3.3.3        tidyverse_1.3.0     

loaded via a namespace (and not attached):
  [1] nlme_3.1-152        matrixStats_0.58.0  fs_1.5.0            lubridate_1.7.9.2   RColorBrewer_1.1-2  httr_1.4.2         
[7] ggplus_0.1          tools_4.0.4         backports_1.2.1     utf8_1.1.4          R6_2.5.0            DBI_1.1.1          
[13] BiocGenerics_0.36.1 colorspace_2.0-0    GetoptLong_1.0.5    withr_2.4.2         tidyselect_1.1.0    mnormt_2.0.2       
[19] curl_4.3            compiler_4.0.4      cli_2.4.0           rvest_1.0.0         Cairo_1.5-12.2      xml2_1.3.2         
[25] labeling_0.4.2      scales_1.1.1        digest_0.6.27       foreign_0.8-81      rio_0.5.26          pkgconfig_2.0.3    
[31] dbplyr_2.1.1        rlang_0.4.10        GlobalOptions_0.1.2 rstudioapi_0.13     farver_2.0.3        shape_1.4.5        
[37] generics_0.1.0      jsonlite_1.7.2      zip_2.1.1           car_3.0-10          magrittr_2.0.1      Rcpp_1.0.6         
[43] munsell_0.5.0       S4Vectors_0.28.1    fansi_0.4.2         abind_1.4-5         lifecycle_1.0.0     stringi_1.5.3      
[49] yaml_2.2.1          carData_3.0-4       plyr_1.8.6          parallel_4.0.4      crayon_1.4.1        lattice_0.20-41    
[55] haven_2.3.1         hms_1.0.0           magick_2.6.0        tmvnsim_1.0-2       pillar_1.6.0        ggpubr_0.4.0       
[61] rjson_0.2.20        ggsignif_0.6.1      stats4_4.0.4        reprex_2.0.0        glue_1.4.2          data.table_1.13.6  
[67] modelr_0.1.8        png_0.1-7           vctrs_0.3.6         cellranger_1.1.0    gtable_0.3.0        clue_0.3-58        
[73] assertthat_0.2.1    xfun_0.21           openxlsx_4.2.3      broom_0.7.6         rstatix_0.7.0       tinytex_0.31       
[79] IRanges_2.24.1      cluster_2.1.1       ellipsis_0.3.1    
