library(openxlsx)
library(ggplot2)
library(ggrepel)


data=read.xlsx("KD-vs-NC.gene.xlsx", sheet=1)#读取工作表格
data$change = ifelse(data$`p-value` < 0.05 & data$FoldChange >= 1.2|data$`p-value` < 0.05&data$FoldChange <= (1/1.2),
                     ifelse(data$FoldChange>= 1.2 ,"Up","Down"),"Stable")
a=table(data$change)
data$`log2(FoldChange)`<-as.numeric(log2(data$FoldChange))
data$FoldChange<-as.numeric(data$FoldChange)
b <- c("ENTPD1","PRUNE1","NT5E")
data$label <- ifelse(data$gene_id %in%b,data$gene_id, "")

ggplot(data, aes(x = `log2(FoldChange)`, y = -log10(`p-value`),colour=change,label=gene_id))+
  geom_point(alpha=1, size=2,shape=16) +#设置颜色深浅以及散点的大小
  scale_color_manual(values=c("blue", "gray50","red"))+#设置点的颜色
  geom_vline(xintercept=c(-log2(1.2),log2(1.2)),lty=2,col="black",lwd=0.8)+ #根据x轴的数据做辅助线，辅助线的位置为xintercept，lty=样式,col=颜色,lwd=粗细
  annotate("text",label="FC=1/1.2",x=-2,y=55,colour="black",size=4)+ #添加文本
  annotate("text",label="FC=1.2",x=2,y=55,colour="black",size=4)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.8)+ # 辅助线y=-log10(0.05)
  annotate("text",label="P-value=0.05",x=8,y=0.7,colour="black",size=4,family="sans")+
  labs(x="log2(FC)",
       y="-log10(P-value)",size=15,family="sans")+# 设置坐标轴名称
  theme_bw()+
  ggtitle("Volcano Plot")+
  theme(legend.position="right",
        text=element_text(size=12,family="sans"),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.background = element_rect(fill = NA,color = "gray",linetype = 1),
        axis.text.y = element_text(colour="black",size=12),
        axis.text.x = element_text(colour="black",size=12),
        plot.title = element_text(hjust = 0.5,size = 20),
        aspect.ratio = 1.1,#更改图形比例
        plot.margin=unit(c(2,2,2,2), "lines"))+
  scale_y_continuous(breaks = seq(0,70,by=10))+
  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 2))+
  scale_colour_manual(values=c("blue", "gray50","red"),
                     labels=c(paste0("down-regulated","(",a[1],")"),
                              paste0("none-Significan","(",a[2],")"),
                              paste0("up-regulated","(",a[3],")")))+
  geom_text_repel(aes(`log2(FoldChange)`,-log10(`p-value`),label = `label`),
                  size=5,arrow = arrow(length = unit(0.02, "npc")),
                  box.padding = 1,segment.size=0.75,
                  show.legend = F,nudge_x =ifelse(data$`log2(FoldChange)`<0,-2,2),
                  color="black",
                  family="sans",
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 1000000))#添加基因名
ggsave("volcano-2.png",width = 10, height =8)
ggsave("volcano-2.pdf",width = 10, height =8)

