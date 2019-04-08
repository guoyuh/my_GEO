## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2018-12-20 15:43:52
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-12-20  First version
###
### ---------------


rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata')
table(group_list)
# 每次都要检测数据
dat[1:4,1:4]

#检查矩阵中管家基因GAPDH ，ACTB的基因表达量是否较高
# boxplot(dat['GAPDH',])
# boxplot(dat[,1])

exprSet=dat
head(exprSet)
library(reshape2)
#exprSet_L=melt(exprSet,id.vars=c("probe")) #宽矩阵转换为长矩阵 
exprSet_L=melt(exprSet)
colnames(exprSet_L)=c('probe','sample','value')
exprSet_L$group=rep(group_list,each=nrow(exprSet))

# the up command is equal the flow command
#exprSet_L$group=c(rep('Control',18834),rep('Control',18834),rep('Control',18834),
#rep('Vemurafenib',18834),rep('Vemurafenib',18834),rep('Vemurafenib',18834))
#or
#exprSet_L$group=c(rep(group_list,18834))
####样本分组笔记#################################
############### 样本分组
##group_list = c(rep('not_TN',ncol(not_TN_expr)),rep('TN', ncol(TN_expr)))   ## 1:265
## rep('not_TN',ncol(not_TN_expr))：连续打和矩阵列数not_TN_expr一样数量的non_TN
## rep('TN', ncol(TN_expr))：连续打和矩阵列数TN_expr一样数量的TN



head(exprSet_L)
### ggplot2 
library(ggplot2)
p=ggplot(exprSet_L,
         aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)
p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_violin()
print(p)
p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 50)+facet_wrap(~sample, nrow = 2)
#bins 分组数量 ref：http://www.sohu.com/a/134994054_683794 （手把手教你使用ggplot2进行数据分布探索 ）
print(p)
p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 2)
print(p)
p=ggplot(exprSet_L,aes(value,col=group))+geom_density() 
print(p)

#
p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
p=p+theme_set(theme_set(theme_bw(base_size=20)))
p=p+theme(text=element_text(face='bold'),
          axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
#axis.text.x  坐标轴标题旋转
print(p)


## hclust 
colnames(exprSet)=paste(group_list,1:6,sep='')
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
plot(hclust(dist(t(exprSet))))
hc=hclust(dist(t(exprSet)))
par(mar=c(5,5,5,10)) 
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)






## 下面是画PCA的必须操作，需要看说明书。
dat=t(dat)#画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换
dat=as.data.frame(dat)#将matrix转换为data.frame
dat=cbind(dat,group_list) #cbind横向追加，即将分组信息追加到最后一列
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 
# The variable group_list (index = 54676) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('all_samples_PCA.png')




rm(list = ls())  ## 魔幻操作，一键清空~
load(file = 'step1-output.Rdata') #此步为一个小插曲，即计算一下从第一行开是计算每一行的sd值，知道最后一行所需要的时间
dat[1:4,1:4] 

cg=names(tail(sort(apply(dat,1,sd)),1000))#apply按行（'1'是按行取，'2'是按列取）取每一行的方差，从小到大排序，取最大的1000个
library(pheatmap)
pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #对那些提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵

####???????????????有必要对表达矩阵转置两次吗，等同于不转置？###################
n=t(scale(t(dat[cg,]))) # 'scale'可以对log-ratio数值进行归一化
#head(n_tmp)
#boxplot(n_tmp[,'GSM1052615'])

head(n)
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
#show_colnames=T ,show_colnames=T 分别显示样本的聚类 和基因的聚类


ac=data.frame(g=group_list)
rownames(ac)=colnames(n) #把ac的行名给到n的列名，即对每一个探针标记上分组信息

pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = 'heatmap_top1000_sd.png')
#annotation_col=ac  添加样本的分组标签
