savepath='ttest_result/'
createdir(savepath)
filname='ttest_result_A6801_A8001_A6823_A8023_vs_oldcon_annotation_differential.xls'

datfile=readdata(filname)%>%unique()
#datfile$`p-adjust`<- p.adjust(datfile$p.value, method = "hochberg")
datfile$ID<- paste0(datfile$Assay,' ',datfile$OlinkID)
pr=datfile$ID
gr=readdata("Sample_Group_1.xlsx",sheet = 3)
cn=c('ID',gr$Sample)#设置需要读取的列名并根据列名筛选列
dd=datfile[datfile$ID%in%pr,colnames(datfile)%in%cn]
rownames(dd)=(dd[,'ID'])
dd=dd[,-which(colnames(dd)=='ID')]
as.character(unique(gr$Group))->b
#根据sample表中设置好的颜色确定箱线图箱子的颜色并命名
col2=unique(gr[,colnames(gr)%in%c("Group","color")])
col=col2$color
names(col)=col2$样品分组
for (i in 1:nrow(dd)) {
  #提取每个蛋白的数据并添加分组信息
  d1<- as.data.frame(t(dd[i,]))
  for (j in 1:nrow(d1)) {
    d1$group[j]=gr$Group[which(rownames(d1)[j]==gr$Sample)]
  }
  d1$Group <- factor(d1$group,levels = b)
  name=strsplit(colnames(d1)[1],' ')[[1]][1]
  
  pval=format(datfile$p.value[which(datfile$ID==colnames(d1)[1])],digits = 4)
  text=paste0('ANOVA','\n','p=',pval)
  p<- ggplot(d1,aes(x=Group,y=d1[,1],fill=Group))+
    geom_violin(alpha = 1,              # 透明度
                trim = F,               # 是否修剪尾巴，即将数据控制到真实的数据范围内
                scale = "area",         # 如果“area”(默认)，所有小提琴都有相同的面积(在修剪尾巴之前)。如果是“count”，区域与观测的数量成比例。如果是“width”，所有的小提琴都有相同的最大宽度。
                 color='white')+
    geom_boxplot(width=0.1,position = position_dodge(0.5))+
    scale_fill_manual(limits=b,
                      values=col)+ #颜色
    annotate("text",label=text,x=1.5,y=max(d1[,1])+1.5,colour="black",size=4)+ #添加文本
    ylim(min(d1[,1])-0.5, max(d1[,1])+3)+
    #scale_y_continuous(limits = c(min(d1[,1])-0.5, max(d1[,1])+3), breaks = seq(0,round(max(d1[,1])+3,0),2))+
    labs(x="Group",y="NPX",title = colnames(d1)[1])+ # 添加标题，x轴，y轴内容
    theme(title = element_text(size=20,color='black',face = 'bold',family = "sans"),
          axis.title.y = element_text(size = 18,family = "sans",face = 'bold',color = "black",vjust = 4,hjust = 0.5),
          legend.title = element_text(colour = 'black',face = 'bold',size = 18,family = 'sans'),
          legend.text = element_text(color="black",size = 15,family = "sans",face = 'bold'), # 设置图例标签文字
          axis.text.x = element_text(size = 15,color = "black",family = "sans",face = 'bold',hjust = 1,angle = 45),# 修改X轴上字体大小，颜色
          axis.text.y = element_text(size = 15,color = "black",family = "sans",face = 'bold',hjust = 0.5),
          axis.line = element_line(colour = 'black',size = 1),
          panel.grid.major =element_blank(),
          panel.background = element_blank(),
          panel.border=element_blank(),
          plot.margin=unit(c(2,2,2,5), "cm"),
          aspect.ratio = 3/2)
  ggplotsave(plot = p,savepath=savepath,
             mapname = paste0(name,'_violin'),
             width = width,
             height = height,
             imagetype = imagetype)
  
}
  
