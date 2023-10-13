#!/opt/conda/bin/Rscript
#' 复制联合分析所需文件
#'
#' @param ptmfile 修饰项目文件路径
#' @param profile 蛋白项目文件路径
#' @param path 运行结果保存路径
#' @export
copy_files<- function(ptmfile = './',profile= './',path = 'rawdata/'){
  pacman::p_load(lmbio,dplyr)
  
  grptm<- list.files(path = ptmfile,pattern = "Sample_Group",recursive = T)
  grpro<- list.files(path = profile,pattern = "Sample_Group",recursive = T)
  group<- list()
  for (i in 1:2) {
    group_pro<- readdata(paste0(profile,'/',grpro),sheet = i)
    group_ptm<- readdata(paste0(ptmfile,'/',grptm),sheet = i)
    group[[i]]<- merge(group_pro,group_ptm)
  }
  names(group)<- c('样品','比较组')
  savexlsx(group,paste0(path,'/Sample_Group.xlsx'))
  
  prodir <- list.files(path = profile,pattern = "^表达矩阵|差异表达矩阵\\(未筛选\\)",recursive = T)
  ptmdir <- list.files(path = ptmfile,pattern = "^表达矩阵|差异表达矩阵\\(未筛选\\)|classtype",recursive = T)
  
  for (r in prodir) {
    namem = strsplit(r,'/')[[1]][length(strsplit(r,'/')[[1]])]%>%gsub('.xlsx','-pro.xlsx',.)
    base::file.copy(from = paste0(profile,"/",r), to = paste0(path,namem))
  }
  for (t in ptmdir) {
    namem = strsplit(t,'/')[[1]][length(strsplit(t,'/')[[1]])]%>%gsub('.xlsx','-ptm.xlsx',.)
    base::file.copy(from = paste0(ptmfile,"/",t), to = paste0(path,namem))
  }
  ##FC筛选条件
  ratio_pro = readxlsx(paste0(profile,'/',
                              list.files(path = profile,pattern = "overview_information",recursive = T)),
                       sheet = '筛选条件')
  ratio_ptm = readxlsx(paste0(ptmfile,'/',
                              list.files(path = ptmfile,pattern = "overview_information",recursive = T)),
                       sheet = '筛选条件')
  colnames(ratio_pro)<- paste0(colnames(ratio_pro),'_pro')
  colnames(ratio_ptm)<- paste0(colnames(ratio_ptm),'_ptm')
  ratio = cbind(ratio_pro,ratio_ptm)
  savexlsx(ratio,paste0(path,'/ratio.xlsx'))
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-pt","--ptmfile",default = './', help = "修饰项目文件路径，例如：/data/hstore3/Projects/202306/ZLM2023010890")
  parser$add_argument("-pr","--profile",default = './', help = "蛋白项目文件路径，例如：/data/hstore3/Projects/202306/ZLM2023010891")
  parser$add_argument("-ip","--path",default = 'rawdata/' ,help = "运行结果保存路径,默认：'rawdata/")
  args <- parser$parse_args()
  copy_files <- do.call(copy_files,args = args)
}
