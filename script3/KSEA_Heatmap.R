#!/opt/conda/bin/Rscript
#' 激酶活性多比较组时绘制热图
#'
#' @param savepath 保存路径
#' @param m.cutoff 底物数量
#' @param stats p值
#' @param imagetype 图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 分辨率
#' @param fontfamily 字体
#' @param ... 
#' @export
KSEA_pheatmap<-function(savepath = "./KSEA_Heatmap/", inputpath = './KSEA_Activity/',m.cutoff=2,stats='p.value' ,imagetype = c("pdf", "png") , 
                        height = 12, width = 10, dpi = 300, fontfamily = "sans", ...){
  
  pacman::p_load(dplyr,stringr,Hmisc,readxl,openxlsx,KSEAapp,ggplot2,optparse,pheatmap,tidyr,lmbio)
  
  con<- dir(path = inputpath, pattern = "*")
  createdir(savepath)
  if(length(con)>1){
    KS<- list()
    for (g in 1:length(con)) {
      KS[[g]]<-read.xlsx(paste0(inputpath,con[g],"/KSEA Kinase Scores.xlsx"))
    }
    sample.labels<- con
    
    filter.m = function(dataset, m.cutoff) {
      filtered = dataset[(dataset$m >= m.cutoff), ]
      return(filtered)
    }
    KS.m = lapply(KS, function(...) filter.m(...,m.cutoff))
    for (i in 1:length(KS.m)) {
      names = colnames(KS.m[[i]])[c(2:7)]
      colnames(KS.m[[i]])[c(2:7)] = paste(names, i,sep = ".")
    }
    master = Reduce(function(...) merge(..., by = "Kinase.Gene",all = F), KS.m)
    if(!nrow(master)==0){
      row.names(master) = master$Kinase.Gene
      columns = as.character(colnames(master))
      merged.scores = as.matrix(master[, grep("z.score", columns)])
      colnames(merged.scores) = sample.labels
      merged.stats = as.matrix(master[, grep(stats, columns)])
      data_mark=merged.stats
      for(i in 1:nrow(merged.stats)){
        for(j in 1:ncol(merged.stats)){
          if (merged.stats[i,j]<0.0001){data_mark[i,j]="****"}else{
            if (merged.stats[i,j]<0.001){data_mark[i,j]="***"}else{
              if (merged.stats[i,j]<0.01){data_mark[i,j]="**"}else{
                if (merged.stats[i,j]<0.05){data_mark[i,j]="*"}else{
                  data_mark[i,j]=""}}}}}}
      bk <- c(seq(floor(min(merged.scores)),-0.1,by=0.001),seq(0,floor(max(merged.scores)),by=0.001))
      pheatmap::pheatmap(merged.scores,display_numbers=data_mark,
                   cluster_rows = F,
                   cluster_cols = F,
                   clustering_method = "complete",
                   fontfamily= fontfamily,
                   fontsize_number=14,
                   number_format ="%.3f",
                   number_color = "black",
                   show_colnames = T,
                   show_rownames = T,
                   angle_col=90,
                   color = c(colorRampPalette(colors = c("blue","white"))(length(which(bk<0))),colorRampPalette(colors = c("white","red"))(length(which(bk>=0)))),
                   legend_breaks=seq(floor(min(merged.scores)),floor(max(merged.scores)),1),
                   breaks=bk,
                   cellwidth =  ifelse(ncol(merged.scores)>3,400/ncol(merged.scores),50),
                   cellheight = ifelse(nrow(merged.scores)>3,ifelse(nrow(merged.scores)<7,150/nrow(merged.scores),350/nrow(merged.scores)),40),
                   scale = "none",
                   fontsize_col =ifelse(ncol(merged.scores)>3,ifelse(ncol(merged.scores)<7,80/ncol(merged.scores),50/ncol(merged.scores)),20),
                   fontsize_row =ifelse(nrow(merged.scores)>3,ifelse(nrow(merged.scores)<7,150/nrow(merged.scores),350/nrow(merged.scores)),40))
      #ggsave("KSEA_Heatmap/KSEA.Merged.Heatmap.pdf", width=ncol(merged.scores)*2, height=12,dpi = 300,p)
      #ggsave("KSEA_Heatmap/KSEA.Merged.Heatmap.png", width=ncol(merged.scores)*2, height=12,dpi = 300,p)
      plotsave(savepath=savepath,
                 mapname = "KSEA.Merged.Heatmap",
                 width = width,
                 height = height,
                 imagetype = imagetype,
                 family=fontfamily,
               units = "in")
      ms<-cbind(as.data.frame(rownames(merged.scores)),as.data.frame(merged.scores))
      colnames(ms)=c("kinase",paste0(sample.labels,"-z.scores"))
      write.xlsx(ms,paste0(savepath,"KSEA热图数据.xlsx"))
    }else{
      savetxt(data = "无交集激酶，不生成交集激酶热图。",
              filename = paste0(savepath,"/说明.txt"),append = T)
    }
    
    ###并集蛋白
    KS.m_nrow = lapply(KS.m, function(...)nrow(...))%>%unlist()
    if(0%in%KS.m_nrow){
      KS.m2 = KS.m[lapply(KS.m,nrow)>0]
      master2 = Reduce(function(...) merge(..., by = "Kinase.Gene",all = T), KS.m2)
      sample.labels2 = sample.labels[-(which(KS.m_nrow==0))]
      savetxt(data = "该项目有比较组底物数量均<2，因此筛选后激酶数量为0，绘图时删去该比较组。",
              filename = paste0(savepath,"/说明.txt"),append = T)
    }else{
      master2 = Reduce(function(...) merge(..., by = "Kinase.Gene",all = T), KS.m)
      sample.labels2 = sample.labels
    }
    row.names(master2) = master2$Kinase.Gene
    columns = as.character(colnames(master2))
    merged.scores2 = as.matrix(master2[, grep("z.score", columns)])
    colnames(merged.scores2) = sample.labels2
    merged.stats2 = as.matrix(master2[, grep(stats, columns)])
    data_mark2=merged.stats2
    for(i in 1:nrow(merged.stats2)){
      for(j in 1:ncol(merged.stats2)){
        if (is.na(merged.stats2[i,j])){data_mark2[i,j]=""}else{
          if (merged.stats2[i,j]<0.0001){data_mark2[i,j]="****"}else{
            if (merged.stats2[i,j]<0.001){data_mark2[i,j]="***"}else{
              if (merged.stats2[i,j]<0.01){data_mark2[i,j]="**"}else{
                if (merged.stats2[i,j]<0.05){data_mark2[i,j]="*"}else{
                  data_mark2[i,j]=""}}}}}}}
    bk2 <- c(seq(floor(min(merged.scores2,na.rm = T)),-0.1,by=0.001),seq(0,floor(max(merged.scores2,na.rm = T)),by=0.001))
    pheatmap::pheatmap(merged.scores2,display_numbers=data_mark2,
                 cluster_rows = F,
                 cluster_cols = F,
                 clustering_method = "complete",
                 fontfamily= 'sans',
                 fontsize_number=ifelse(nrow(merged.scores2)>25,680/nrow(merged.scores2),ifelse(nrow(merged.scores2)<10,250/nrow(merged.scores2),360/nrow(merged.scores2))),
                 number_format ="%.3f",
                 number_color = "black",
                 show_colnames = T,
                 show_rownames = T,
                 angle_col=90,
                 color = c(colorRampPalette(colors = c("blue","white"))(length(which(bk2<0))),colorRampPalette(colors = c("white","red"))(length(which(bk2>=0)))),
                 legend_breaks=seq(floor(min(merged.scores2,na.rm = T)),floor(max(merged.scores2,na.rm = T)),1),
                 breaks=bk2,
                 cellwidth =  ifelse(ncol(merged.scores2)>3,400/ncol(merged.scores2),50),
                 cellheight = ifelse(nrow(merged.scores2)>25,680/nrow(merged.scores2),350/nrow(merged.scores2)),
                 scale = "none",
                 fontsize_col =ifelse(ncol(merged.scores2)>3,ifelse(ncol(merged.scores2)<7,80/ncol(merged.scores2),50/ncol(merged.scores2)),20),
                 fontsize_row =ifelse(nrow(merged.scores2)>25,580/nrow(merged.scores2),ifelse(nrow(merged.scores2)<10,150/nrow(merged.scores2),260/nrow(merged.scores2))))
    plotsave(savepath=savepath,
             mapname = "KSEA.Merged_union.Heatmap",
             width = width,
             height = height,
             imagetype = imagetype,
             family=fontfamily,
             units = "in")
    ms2<-cbind(as.data.frame(rownames(merged.scores2)),as.data.frame(merged.scores2))
    colnames(ms2)=c("kinase",paste0(sample.labels2,"-z.scores"))
    write.xlsx(ms2,paste0(savepath,"KSEA热图数据-union.xlsx"))
    
    write.xlsx(KS, file = paste0(savepath,"phos_kinase-source.xlsx"),sheetName=sample.labels)
  }else{
    savetxt(data = "该项目为单个比较组，不进行热图绘制",
            filename = paste0(savepath,"/说明.txt"),append = T)
    return()
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 10, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 14, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-s","--savepath",default = "./KSEA_Heatmap/", help = "绘图存放分析文件夹路径，默认./KSEA_Heatmap/")
  parser$add_argument("-ip","--inputpath",default = "./KSEA_Activity/", help = "绘图数据导入路径，默认./KSEA_Activity/")
  parser$add_argument("-m","--m.cutoff",default=2,help="激酶对应底物的最小数量,默认2")
  parser$add_argument("-v","--stats",default='p.value',help="富集展示方法,默认'p.value'")
  args <- parser$parse_args()
  KSEA_pheatmap <- do.call(KSEA_pheatmap,args = args)
}
