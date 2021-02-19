## heatmap
cpdb_heatmap("test_out/means.txt", "test_out/pvalues.txt", "test_meta.txt", 'test_out/tmp/')

cpdb_heatmap <- function(mean_file, pvalues_file, meta_file, outputRoot=""){ #用不到 mean_file
  # https://github.com/MISAKA-DaYu/CellPhoneDB.plots/blob/main/cpdb.R
  if(0){ #测试
    mean_file="test_out/means.txt"
    pvalues_file="test_out/pvalues.txt"
    
    meta_file="test_meta.txt"
    
    dir= paste0(getwd(), "/tmp/"); dir
  }else{
    dir=outputRoot;
  }
  
  library(pheatmap)

  count_filename <- paste0(dir,'heatmap_count.pdf')
  log_filename <- paste0(dir,'heatmap_log_count.pdf')
  count_network_filename <- paste0(dir,'count_network.txt')
  interaction_count_filename <- paste0(dir,'interactions_count.txt')
  
  count_network_separator <- '|'
  interaction_count_separator <- '|'
  #######   Network
  show_rownames = T
  show_colnames = T
  scale="none"
  cluster_cols = T
  border_color='white'
  cluster_rows = T
  fontsize_row=12
  fontsize_col = 12 
  main = ''
  treeheight_row=0
  family='Arial'
  treeheight_col = 0
  col1 = "dodgerblue4"
  col2 = 'peachpuff'
  col3 = 'deeppink4'
  meta_sep='\t'
  pvalues_sep='\t'
  pvalue=0.05
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1)){
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  all_count = all_count[-1,]
  write.table(all_count, count_network_filename, sep=count_network_separator, quote=F, row.names = F)
  
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1)){
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  if (any(count1)>0){
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])
    saveRDS(count_matrix,"heatmap_count.rds")
    saveRDS(log(count_matrix+1),"heatmap_log_count.rds")
    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    write.table(all_sum, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=F)
    
    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
    
    pheatmap::pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
                       border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                       main = main, cellwidth = 30, cellheight = 30, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)
    
    pheatmap::pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
                       border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                       main = main, cellwidth = 30, cellheight = 30, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}




########################
## dotplot
########################
cpdb_dotplot("test_out/means.txt", "test_out/pvalues.txt", "test_meta.txt", 'test_out/tmp/')

cpdb_dotplot <- function(mean_file, pvalues_file, meta_file, outputRoot=""){ #用不到 meta_file
  library(ggplot2)
  library(reshape)
  library(RColorBrewer)
  library(grid)
  library(Cairo)
  #setwd(dir)
  dir=outputRoot
  
  mfile <- read.csv(mean_file,sep="\t",check.names = F)
  rownames(mfile) <- mfile$interacting_pair
  pairs <- colnames(mfile)[-c(1:11)]
  # pairs <- pairs[c(grep("Macrophage",pairs),grep("DC",pairs))]
  ## step1: filter pfile, mfile
  # 筛选出有mean值的interacting_pair
  mfile[is.na(mfile)]=0
  mfile <- mfile[,as.character(pairs)]
  msum <- rowSums(mfile)
  #class(msum)
  msum <- as.data.frame(msum[msum >0])
  colnames(msum) <- 'mean_sum'
  msum$interacting_pair <- rownames(msum)
  # msum <- as.data.frame(msum)
  m_fi <- mfile[rownames(msum),]
  # m_fi[1:5,1:10]
  
  ## p-value prepare
  pfile <- read.csv( pvalues_file, sep='\t',check.names=F)
  rownames(pfile) <- pfile$interacting_pair
  pfile<-pfile[,colnames(m_fi)]
  # 根据m_fi筛选行，再根据p值筛选列（cluster）
  pfile <- pfile[rownames(m_fi),]
  dim(pfile)
  # 筛选出含有p值<0.05的cluster_pair
  l <- ''
  for (i in colnames(pfile)){
    if (min(pfile[,i])<=0.05){
      #print(i)
      l <- c(l,i)
    }
  }
  class(l)
  l <- l[-1]
  l <- as.character(l)
  pfile <- pfile[,l]
  dim(pfile)
  # 筛选含有p值<0.01的受体配体对
  r <- ""
  for (i in rownames(pfile)){
    if (min(pfile[i,])<=0.01){
      r <- c(r,i)
    }
  }
  r <- r[-1]
  length(r)
  r <- as.character(r)
  pfile <- pfile[r,]
  dim(pfile)
  pfile[1:5,1:5]
  m_fi<- m_fi[r,]
  write.csv(m_fi,paste0(dir,'means_filtered.csv'))
  write.csv(pfile,paste0(dir,'pvalue_filtered.csv'))
  ## 取p值之和最小的前50个受体配体对
  pfile$psum <- rowSums(pfile)
  pfile <- pfile[order(pfile$psum, decreasing = F),]
  pfileraw <- pfile
  pfile <- pfile[,-ncol(pfile)]
  pfile <- pfile[1:50,]
  r_new <- rownames(pfile)
  m_fi <- m_fi[r_new,]
  m_fi <- m_fi[,l]
  
  facets<-unique(unlist(strsplit(colnames(pfile),split='\\|')))
  selfid<-paste(facets,facets,sep="|")
  
  m.stat <- as.data.frame(colSums(m_fi))
  #head(m.stat)
  colnames(m.stat) <- 'sum_means'
  m.stat$interaction_cluster <- rownames(m.stat)
  m.stat <- m.stat[order(m.stat$sum_means,decreasing = TRUE),]
  write.csv(m.stat,paste0(dir,'means_stat_by_clusters.top50.csv'),row.names = FALSE)
  
  interacting_pair <- "interacting_pair"
  interacting_cluster <- "interacting_cluster"
  p_value <- "p_value"
  mean <- "mean"
  r <- rownames(pfile)
  
  for (x in l){
    #print(x)
    interacting_pair <- c(interacting_pair,r)
    interacting_cluster <- c(interacting_cluster,rep(x,time=length(r)))
    p_value <- c(p_value,pfile[,x])
    mean <- c(mean, m_fi[,x])
  }
  dat <- data.frame(interacting_pair,interacting_cluster,p_value,mean)
  #head(dat)
  dat <- dat[-1,]
  dat$p_value <- as.character(dat$p_value)
  dat$mean <- as.numeric(dat$mean)
  class(dat$mean)
  dat$log10mean <- log10(dat$mean)
  dat[,'-log10pvalue'] <- -log10(as.numeric(dat$p_value))
  max(dat[,'-log10pvalue'])
  min(dat[,'-log10pvalue'])
  dat$new_p[dat['-log10pvalue']>=0 & dat['-log10pvalue']<1] ='0'
  dat$new_p[dat['-log10pvalue']>=1 & dat['-log10pvalue']<2] ='1'
  dat$new_p[dat['-log10pvalue']>=2 & dat['-log10pvalue']<3] ='2'
  dat$new_p[dat['-log10pvalue']>=3] ='>=3'
  dat$new_p <- factor(dat$new_p, levels=c('0','1','2','>=3'))
  
  dat$p_group[dat$p_value>=0 & dat$p_value < 0.001] <- "0.001"
  dat$p_group[dat$p_value>=0.001 & dat$p_value < 0.01] <- "0.01"
  dat$p_group[dat$p_value>=0.01 & dat$p_value <= 0.05] <- "0.05"
  dat$p_group[dat$p_value>0.05] <- "NotSig"
  dat$p_group <- factor(dat$p_group, levels = c("NotSig", '0.05', '0.01','0.001'))
  # rownames(dat) <- dat$interacting_pair
  write.csv(dat,paste0(dir,'DotPlot.top50.csv'))
  #width=length(unique(dat$interacting_pair))*10
  #height=length(unique(dat$interacting_cluster))*10
  selfDat<-dat[dat$interacting_cluster %in% selfid,drop=F,]
  heightself=length(unique(selfDat$interacting_cluster))*66
  widthself=length(unique(selfDat$interacting_pair))*60
  #dat<-dat[!colnames(dat) %in% selfid, drop = F,]
  # dat<-dat[!dat$interacting_cluster %in% selfid,drop=F,]
  heightall=length(unique(dat$interacting_cluster))*66
  widthall=length(unique(dat$interacting_pair))*60
  # write.csv(selfDat,paste0(dir,'/dot_selfDot.csv'))
  mytheme <- theme_bw() + theme(plot.title = element_blank(),
                                panel.grid= element_blank(),
                                axis.title = element_blank(),
                                axis.ticks.length = unit(0.3, "cm"),
                                axis.text=element_text(family="ArialMT", size = 10),
                                axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size=10),
                                legend.title=element_text(family="ArialMT", size=10),
                                legend.text=element_text(family="ArialMT", size=10),
                                legend.background=element_rect(fill="transparent"),
                                plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
  mycol <- brewer.pal(11,"RdYlBu")
  
  pdf(file=paste0(dir,"dotAll_plt.top50.pdf"),width=15,height=15)
  print(ggplot(dat, aes(x=interacting_pair, y=interacting_cluster, colour = log10mean, size=p_group)) + geom_point() +
          scale_color_gradientn(name='log10(mean)', colors=rev(mycol)) +
          coord_flip()+
          mytheme)
  dev.off()
  
  ## 取p值之和最小的前20个受体配体对
  pfile <- pfile[1:20,]
  r_new <- rownames(pfile)
  m_fi <- m_fi[r_new,]
  m_fi <- m_fi[,l]
  
  facets<-unique(unlist(strsplit(colnames(pfile),split='\\|')))
  selfid<-paste(facets,facets,sep="|")
  
  m.stat <- as.data.frame(colSums(m_fi))
  #head(m.stat)
  colnames(m.stat) <- 'sum_means'
  m.stat$interaction_cluster <- rownames(m.stat)
  m.stat <- m.stat[order(m.stat$sum_means,decreasing = TRUE),]
  write.csv(m.stat,paste0(dir,'means_stat_by_clusters.top20.csv'),row.names = FALSE)
  
  interacting_pair <- "interacting_pair"
  interacting_cluster <- "interacting_cluster"
  p_value <- "p_value"
  mean <- "mean"
  r <- rownames(pfile)
  
  for (x in l){
    #print(x)
    interacting_pair <- c(interacting_pair,r)
    interacting_cluster <- c(interacting_cluster,rep(x,time=length(r)))
    p_value <- c(p_value,pfile[,x])
    mean <- c(mean, m_fi[,x])
  }
  dat <- data.frame(interacting_pair,interacting_cluster,p_value,mean)
  #head(dat)
  dat <- dat[-1,]
  dat$p_value <- as.character(dat$p_value)
  dat$mean <- as.numeric(dat$mean)
  class(dat$mean)
  dat$log10mean <- log10(dat$mean)
  dat[,'-log10pvalue'] <- -log10(as.numeric(dat$p_value))
  max(dat[,'-log10pvalue'])
  min(dat[,'-log10pvalue'])
  dat$new_p[dat['-log10pvalue']>=0 & dat['-log10pvalue']<1] ='0'
  dat$new_p[dat['-log10pvalue']>=1 & dat['-log10pvalue']<2] ='1'
  dat$new_p[dat['-log10pvalue']>=2 & dat['-log10pvalue']<3] ='2'
  dat$new_p[dat['-log10pvalue']>=3] ='>=3'
  dat$new_p <- factor(dat$new_p, levels=c('0','1','2','>=3'))
  
  dat$p_group[dat$p_value>=0 & dat$p_value < 0.001] <- "0.001"
  dat$p_group[dat$p_value>=0.001 & dat$p_value < 0.01] <- "0.01"
  dat$p_group[dat$p_value>=0.01 & dat$p_value <= 0.05] <- "0.05"
  dat$p_group[dat$p_value>0.05] <- "NotSig"
  dat$p_group <- factor(dat$p_group, levels = c("NotSig", '0.05', '0.01','0.001'))
  # rownames(dat) <- dat$interacting_pair
  write.csv(dat,paste0(dir,'DotPlot.top20.csv'))
  #width=length(unique(dat$interacting_pair))*10
  #height=length(unique(dat$interacting_cluster))*10
  selfDat<-dat[dat$interacting_cluster %in% selfid,drop=F,]
  heightself=length(unique(selfDat$interacting_cluster))*66
  widthself=length(unique(selfDat$interacting_pair))*60
  #dat<-dat[!colnames(dat) %in% selfid, drop = F,]
  # dat<-dat[!dat$interacting_cluster %in% selfid,drop=F,]
  heightall=length(unique(dat$interacting_cluster))*66
  widthall=length(unique(dat$interacting_pair))*60
  write.csv(selfDat,paste0(dir,'dot_selfDot.csv'))
  
  mytheme <- theme_bw() + theme(plot.title = element_blank(),
                                panel.grid= element_blank(),
                                axis.title = element_blank(),
                                axis.ticks.length = unit(0.3, "cm"),
                                axis.text=element_text(family="ArialMT", size = 10),
                                axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size=10),
                                legend.title=element_text(family="ArialMT", size=10),
                                legend.text=element_text(family="ArialMT", size=10),
                                legend.background=element_rect(fill="transparent"),
                                plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
  mycol <- brewer.pal(11,"RdYlBu")
  p=ggplot(dat, aes(x=interacting_pair, y=interacting_cluster, colour = log10mean, size=p_group)) + geom_point() +
          scale_color_gradientn(name='log10(mean)', colors=rev(mycol)) +
          coord_flip()+
          mytheme
  pdf(file=paste0(dir,"dotAll_plt.top20.pdf"),width=8,height=8)
  print(p)
  dev.off()
  return(p)
}

#
ggsave(file="out/R_dotAll_plt.top20-.pdf", p1, width=45, height=7)
#
pdf(file="out/R_dotAll_plt.top20-.pdf", width=75, height=7)
p1
dev.off()
