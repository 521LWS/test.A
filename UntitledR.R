package_name  <- "GEOmirror"
if (!requireNamespace(package_name , quietly = TRUE))
BiocManager::install(package_name )



remotes::install_github("jmzeng1314/GEOmirror")
#三个包同时加载
library(AnnoProbe)
library(GEOmirror)
library(GEOquery)
#下载获取GSE39582数据 
gset=AnnoProbe::geoChina('GSE39582')































