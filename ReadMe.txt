1.准备
(1)download data
https://support.10xgenomics.com/single-cell-gene-expression/datasets

(2) 教程
https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247496154&idx=3&sn=d3cfaa4a5b18235e0192619f64641635



2. 开始
PBMC:
Single vs Dual Indexing Demonstration (v3.1 Chemistry)
Cell Ranger 4.0.0
(1) 10k Peripheral blood mononuclear cells (PBMCs) from a healthy donor, Single Indexed
filtered_feature_bc_matrix:
https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_filtered_feature_bc_matrix.tar.gz
store: ./backup/

a1: Seurat;
a2: Seurat recluster;
a3: Seurat 可视化

b1: monocle;

c1: GO & KEGG;
	cell cluster 之间的两两比较，也能很好的发现区分度很高的marker。
	options(repr.plot.width=12.5, repr.plot.height=3.5) #在jupyter中控制图像的大小

(2) adv: 包括多样本整合、转录因子分析、细胞通讯分析、基因集变异分析和更全面的基因集富集分析


