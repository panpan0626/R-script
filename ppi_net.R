#加载包
library(dplyr)
library(igraph)
library(RColorBrewer)
library(readxl)

#创建绘图文件夹
#if (dir.exists("ppi")){
  #list.files("ppi")
#}else{
  #dir.create("ppi")
#}
newnot <- read.table("ppi_nodeFoldchange.xls",header = T, sep = "\t",na.strings = "")#读取数据表格
data2 <- read.table("ppi_nodes.xls",header = T, sep = "\t")
num<- c(data2$node1,data2$node2)
dd<- data.frame(table(num))
cc<- c()
for (w in 1:nrow(newnot)) {
  cc[w]<-if(newnot$Accession[w]%in%dd$num){dd[which(dd$num%in%newnot$Accession[w]),2]}else{0}
}
newnot$degree<- cc
data1<- newnot

data2 <- read.table("ppi_nodes.xls",header = T, sep = "\t")
colnames(data2) <-c("Accession","Accession2","score")#数据表格重命名列
colnames(data1) <-c("Accession","Gene_name","FC","degree")
#用蛋白名填补缺失的基因名
data1$Gene_name[which(data1$Gene_name=="  ")]=data1$Accession[which(data1$Gene_name=="  ")]
data1 <- data1%>%distinct(Accession,.keep_all = T)#根据蛋白名删除重复项
da_p <- data1[order(-data1$degree),]%>%head(25)#将data1按degree降序排列后取连接度前25的数据
da_p<- da_p[which(!da_p$degree==0),]
da_p2 <- data1[order(-data1$degree),]%>%head(25)#再将data1按degree降序排列后取连接度前25的数据并命名为Accession用于第二次匹配
da_p2<- da_p2[which(!da_p2$degree==0),]

colnames(da_p2)<- c("Accession2","Gene_name","FC","degree")
da_p$notes <- ifelse(da_p$FC>1,"up","down")#增加上下调的标签
dd<-merge(da_p,data2,by="Accession")#将da_p与data2按照Accession列连接为dd
ddd<-merge(da_p2,dd,by="Accession2")#将da_p2与dd连接为ddd

#以蛋白名绘图
da_pr_n <- c("Accession2","Accession","score")
dd_pro <- ddd[,colnames(ddd)%in%da_pr_n]#筛选列
na.omit(dd_pro)%>%as_tibble()->dd_pro#删去含有NA的列
nodes_p <- da_p[,colnames(da_p)%in%c("Accession","FC","notes")]
colnames(nodes_p)<-c("label","FC","location")
nodes_p=nodes_p[order(nodes_p$FC),]
DP<-dd_pro%>%rename(from=Accession2,to=Accession,weight=score)
net_DP<-graph_from_data_frame(d=DP,vertices=nodes_p,directed=TRUE)
###计算节点的度并做标准化处理
deg_p <- as.matrix(scale(degree(net_DP,mode="all"), center = F, scale = T))
#设置节点颜色
#判断是否有FC值，若没有则节点颜色为灰色，若有则上调为红色，下调为绿色，差异倍数越大，颜色越深
#
if (is.na(nodes_p$FC[1])) {V(net_DP)$color <- "gray50"}else{
  if("down"%in%nodes_p$location&'up'%in%nodes_p$location){col_up=c("#FF0000","#8B2323")
col_down=c("#228B22","#00FF00")
col_down1 <- colorRampPalette(col_down)(table(nodes_p$location)[1])
col_down2 <- colorRampPalette(col_up)(table(nodes_p$location)[2])
V(net_DP)$color <- c(col_down1,col_down2)}else{
  if("down"%in%nodes_p$location&!'up'%in%nodes_p$location){col_down=c("#228B22","#00FF00")
  col_down1 <- colorRampPalette(col_down)(table(nodes_p$location))
  V(net_DP)$color <- c(col_down1)}else{if(!"down"%in%nodes_p$location&'up'%in%nodes_p$location){col_up=c("#FF0000","#8B2323")
  col_down2 <- colorRampPalette(col_up)(table(nodes_p$location))
  V(net_DP)$color <- c(col_down2)}}}}

par(mfrow=c(1,1), mar=c(1,1,1,1))
png(file="ppi_query.png",width = 2400,height = 2400,res = 300)
plot(net_DP,layout=do.call("layout_in_circle", list(net_DP)),
     vertex.frame.color=F,
     vertex.size=(deg_p+2)*2,
     vertex.label=V(net_DP)$name,
     vertex.label.dist=1.2,
     vertex.label.cex=0.8,
     vertex.label.color="black",
     vertex.label.family="sans",
     vertex.label.degree=0408,
     edge.arrow.mode=0,
     edge.curved=0,
     edge.color="gray")
dev.off()
pdf(file = "ppi_query.pdf",width =12,height = 12)
plot(net_DP,layout=do.call("layout_in_circle", list(net_DP)),
     vertex.frame.color=F,
     vertex.size=deg_p*8,
     edge.arrow.mode=0,
     vertex.label.dist=1.2,
     vertex.label.cex=0.8,
     vertex.label.color="black",
     vertex.label.family="sans",
     vertex.label.degree=0408,
     edge.curved=0,
     edge.color="gray")
dev.off()

#以基因名绘图
da_ge_n <- c("Gene_name.x","Gene_name.y","score")
dd_ge <- ddd[,colnames(ddd)%in%da_ge_n]#按筛选列
na.omit(dd_ge)%>%as_tibble()->dd_ge
nodes_g <- da_p[,colnames(da_p)%in%c("Gene_name","FC","notes")]
colnames(nodes_g)<-c("label","FC","location")
nodes_g=nodes_g[order(nodes_g$FC),]
DG<-dd_ge%>%rename(from=Gene_name.x,to=Gene_name.y,weight=score)
net_DG<-graph_from_data_frame(d=DG,vertices=nodes_g,directed=TRUE)
###计算节点的度
deg_g <- as.matrix(scale(degree(net_DG,mode="all"), center = F, scale = T))
#设置节点颜色
#判断是否有FC值，若没有则节点颜色为灰色，若有则上调为红色，下调为绿色，差异倍数越大，颜色越深
if (is.na(nodes_g$FC[1])) {V(net_DG)$color <- "gray50"}else{
  if("down"%in%nodes_g$location&'up'%in%nodes_g$location){col_up=c("#FF0000","#8B2323")
  col_down=c("#228B22","#00FF00")
  col_down3 <- colorRampPalette(col_down)(table(nodes_g$location)[1])
  col_down4 <- colorRampPalette(col_up)(table(nodes_g$location)[2])
  V(net_DG)$color <- c(col_down3,col_down4)}else{
    if("down"%in%nodes_g$location&!'up'%in%nodes_g$location){col_down=c("#228B22","#00FF00")
    col_down3 <- colorRampPalette(col_down)(table(nodes_g$location))
    V(net_DG)$color <- c(col_down3)}else{if(!"down"%in%nodes_g$location&'up'%in%nodes_g$location){col_up=c("#FF0000","#8B2323")
    col_down4 <- colorRampPalette(col_up)(table(nodes_g$location))
    V(net_DG)$color <- c(col_down4)}}}}
par(mfrow=c(1,1), mar=c(1,1,1,1))
png(file="ppi_gene.png",width = 2400,height = 2400,res = 300)
plot(net_DG,layout=do.call("layout_in_circle", list(net_DG)),
     vertex.size=(deg_p+2)*2,
     vertex.frame.color=F,
     edge.arrow.mode=0,
     edge.width=1,
     vertex.label.dist=1.2,
     vertex.label.cex=0.8,
     vertex.label.color="black",
     vertex.label.family="sans",
     vertex.label.degree=0408,
     edge.curved=0)
dev.off()
pdf(file = "ppi_gene.pdf",width =12,height = 12)
plot(net_DG,layout=do.call("layout_in_circle", list(net_DG)),
     vertex.size=(deg_p+2)*2,
     vertex.frame.color=F,
     edge.arrow.mode=0,
     edge.width=1,
     vertex.label.dist=1.2,
     vertex.label.cex=0.8,
     vertex.label.color="black",
     vertex.label.family="sans",
     vertex.label.degree=0408,
     edge.curved=0)
dev.off()

