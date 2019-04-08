.libPaths()
setwd("D:/R_WORKHOME/GEO")

#BiocManager::install(c("genefu","org.Hs.eg.db","hgu133plus2.db" ),ask = F,update = F)
#Error in loadNamespace(name) : 不存在叫‘BiocManager’这个名字的程辑包
#所以升级R
#############################updateR################################################
install.packages("installr")
library(installr)
updater()
################################################################################

rm(list = ls())
options()$repos
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos
options()$BioC_mirror




# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocGenerics", version = "3.8")
if(!require("KEGG.db")) BiocManager::install("KEGG.db",ask = F,update = F)
if(!require(c("GSEABase","GSVA","clusterProfiler" ))) BiocManager::install(c("GSEABase","GSVA","clusterProfiler" ),ask = F,update = F)
if(!require(c("GEOquery","limma","impute" )))BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
if(!require(c("genefu","org.Hs.eg.db","hgu133plus2.db" ))) BiocManager::install(c("genefu","org.Hs.eg.db","hgu133plus2.db" ),ask = F,update = F)


######test 2019-03-31#########
# source("https://bioconductor.org/biocLite.R")
# biocLite('clusterProfiler')
# install.packages("githubinstall")
# library("BiocInstaller")
# biocLite('BiocGenerics')
# biocLite('KEGG.db')
#anzuang KEGG.db ,GO.db baocuo
#* installing *source* package 'KEGG.db' ...
#** R
#** inst
#** byte-compile and prepare package for lazy loading
#Error : package 'BiocGenerics' 0.24.0 is loaded, but >= 0.27.1 is required by 'Biobase'
#ERROR: lazy loading failed for package 'KEGG.db'
#* removing 'D:/Program Files/R/R-3.5.3/library/KEGG.db'
#In R CMD INSTALL

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BiocGenerics", version = "3.8")
# 查看r包安装版本   packageVersion("BiocGenerics")
#终于成功安装"BiocGenerics", version = "3.8"
############

# source("https://bioconductor.org/biocLite.R")
# library('BiocInstaller')
# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# BiocInstaller::biocLite("GEOquery")
# BiocInstaller::biocLite(c("limma"))
# BiocInstaller::biocLite(c("impute"))

#options()$repos
#install.packages('WGCNA')
#install.packages(c("FactoMineR", "factoextra"))
#install.packages(c("ggplot2", "pheatmap","ggpubr"))
options()$repos
if(!require("WGCNA")) install.packages('WGCNA')
if(!require(c("FactoMineR", "factoextra"))) install.packages(c("FactoMineR", "factoextra"))
if(!require(c("ggplot2", "pheatmap","ggpubr"))) install.packages(c("ggplot2", "pheatmap","ggpubr"))

library("FactoMineR")
library("factoextra")

library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(genefu)
library(ggplot2)
library(ggpubr)
library(hgu133plus2.db)
library(limma)
library(org.Hs.eg.db)
library(pheatmap)

# 只有这个报错
# library(clusterProfiler)
# Error: package or namespace load failed for ‘clusterProfiler’ in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]):
#   不存在叫‘DO.db’这个名字的程辑包
library("BiocInstaller")
biocLite('DO.db')
# Warning message:
#   'biocLite' is deprecated.
# Use 'BiocManager::install' instead.
# See help("Deprecated")
if(!require("DO.db")) BiocManager::install("DO.db",ask = F,update = F)

#finnally done

