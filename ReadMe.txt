1.准备
(1)download data
https://support.10xgenomics.com/single-cell-gene-expression/datasets

(2) 教程
https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247496154&idx=3&sn=d3cfaa4a5b18235e0192619f64641635
SCENIC: https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247488383&idx=1&sn=7b8504ed4449df3a707d1c83ec0b0a7a

https://data.humancellatlas.org/

(3) 目录结构
scripts_.R_or_Py
backup/
 |-a1/
 |-a2/
 |-a3/

可视化实例:
R_plot_base/image_heatmap_nature2020.R.ipynb


(4)
init code:

setwd("/data/wangjl/scScripts/")
getwd()

subDir="backup/a2/"
if( ! dir.exists( subDir ) ){
    dir.create( subDir )
}
outputRoot=paste0( getwd(),"/", subDir)
outputRoot

Sys.time() #"2021-02-20 10:52:38 CST"



2. 开始
PBMC:
Single vs Dual Indexing Demonstration (v3.1 Chemistry)
Cell Ranger 4.0.0
(1) 10k Peripheral blood mononuclear cells (PBMCs) from a healthy donor, Single Indexed
filtered_feature_bc_matrix:
https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_filtered_feature_bc_matrix.tar.gz
store: ./backup/


## 基本分析：聚类、分群
a1: Seurat;
a2: Seurat recluster; redo v2;
a3: Seurat 可视化
a4: Scanpy – Single-Cell Analysis in Python


## 高级分析：轨迹分析、细胞通信分析、转录因子分析
b1: monocle;
b12: slingShot;


c1: GO & KEGG;
	cell cluster 之间的两两比较，也能很好的发现区分度很高的marker。
	options(repr.plot.width=12.5, repr.plot.height=3.5) #在jupyter中控制图像的大小

(2) adv: 包括多样本整合、转录因子分析、细胞通讯分析、基因集变异分析和更全面的基因集富集分析
b2: adv_cell_commu


(3) SCENIC: 
SCENIC转录因子分析结果的解读:
https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247497665&idx=1&sn=74ac0e87b9689d5df7c0208e1c1dc0ac


(4) 发育过程分析
https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html



