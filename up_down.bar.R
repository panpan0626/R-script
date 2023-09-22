library(openxlsx)
library(readxl)
library(stringr)
library(ggplot2)
library(reshape2)

data=read.xlsx("sample_information.xlsx",sheet = "比较组信息")#读取比较组信息
a=gsub(pattern = "/", replacement = "_", data[,1])
#dataR为差异表格中上下调数量
dataR=data.frame()
b <- dir(path = "./",pattern = ".xlsx")
b1 <- grep("筛选结果",b,value = T)
d <- if (b1=="差异蛋白筛选结果.xlsx") {"Number Of Proteins"}else{
  if (b1=="差异位点筛选结果.xlsx") {"Number Of Sites"}else{"Number Of Peptides"}
}
for (i in 1:length(a)) {
  diff=read_excel(b1,sheet=a[i])%>%as.data.frame()
  dataR[i,1] = a[i]
  dataR[i,2] = length(which(!as.numeric(diff$FC)<1))
  dataR[i,3] = length(which(as.numeric(diff$FC)<1))
}
colnames(dataR)= c("group","up","down")
rownames(dataR)=dataR$group
if(length(a)>10){
  x<- ceiling(length(a)*0.1)
  y<- ceiling(length(a)/x)
  a[(length(a)+1):(x*y)]<- NA
  dd<- matrix(a, nrow = y, ncol = x)
}else{
  dd<- matrix(a)
}
for(k in 1:ncol(dd)) {
  class<- dd[,k]
  data1<- dataR[class,]%>%na.omit()
  df = melt(data1,id= "group")
  df$variable<-factor(df$variable,levels = c("up","down"))
  df$group<-factor(df$group,levels = c(class))
  if(length(a)==1) {
  p=ggplot(df,aes(x=group,y=value,fill=variable))+
      geom_bar( stat="identity",position = position_dodge(0.3),width = 0.3)+
      labs(x = "", y = d)+
      geom_text(aes(label=value), size=5,
                position = position_dodge(width = 0.3), vjust=-0.5)+
      scale_fill_manual(values = c("red","blue"))+
      #scale_y_continuous(expand = c(0,0),limits = c(0, e),breaks = seq(0,e,e1))+
      theme(text=element_text(family="sans"),
            panel.grid = element_blank(),
            panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
            panel.background = element_blank(),
            axis.text.x = element_text(colour="black",hjust =1,angle = 45,size = ifelse(nchar(max(as.character(df$group)))>25,(ifelse(lnchar(max(as.character(df$group)))>30,9,11)),15)),
            axis.text.y = element_text(colour="black",size = 16),
            axis.title.y = element_text(margin = margin(0,1,0,1,"cm"),size = 20),
            axis.ticks.y = element_line(color = "black",size = 1),
            legend.text = element_text(size = 14),
            plot.margin=unit(c(2,2,2,2), "lines"),
            aspect.ratio = 3/2)+
      guides(fill=guide_legend(title=""))
    ggsave("foldchange_bars.png",p,width =11, height =8.5,units = "in")
    ggsave("foldchange_bars.pdf",p,width =11, height =8.5,units = "in")
  }else{
  p=ggplot(df,aes(x=group,y=value,fill=variable))+
      geom_bar( stat="identity",position = position_dodge(0.5),width = 0.5)+
      labs(x = "", y = d)+
      geom_text(aes(label=value), size=4,
                position = position_dodge(width = 0.5), vjust=-0.5)+
      scale_fill_manual(values = c("red","blue"))+
      #scale_y_continuous(expand = c(0,0),limits = c(0, e),breaks = seq(0,e,e1))+
      theme(panel.grid = element_blank(),
            panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
            panel.background = element_blank(),
            axis.text.x = element_text(colour="black",hjust =1,angle = 45,size = ifelse(nchar(max(as.character(df$group)))>25,(ifelse(lnchar(max(as.character(df$group)))>30,9,11)),15)),
            axis.text.y = element_text(colour="black",size = 16),
            axis.title.y = element_text(margin = margin(0,1,0,1,"cm"),size = 20),
            axis.ticks.y = element_line(color = "black",size = 1),
            legend.text = element_text(size = 14),
            plot.margin=unit(c(2,2,2,2), "lines"))+
            #aspect.ratio = 4/5)+
      guides(fill=guide_legend(title=""))
  ggsave(paste0("foldchange_bars-",k,".png"),p,width =11, height =8.5,units = "in")
  ggsave(paste0("foldchange_bars-",k,".pdf"),p,width =11, height =8.5,units = "in")
  }
}
