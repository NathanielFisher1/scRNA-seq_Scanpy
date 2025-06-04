#all necessary modules imported
import pandas as pd #for data clustering
import numpy as np #for math and matrixes
import matplotlib.pyplot as plt #for plotting
import seaborn as sns #for plotting
sns.set() 
import scipy.stats as stats #for math/statistics
from sklearn.decomposition import PCA #for PCA
from sklearn import manifold, datasets # for TSNA
import scanpy as sc #for single cell methods

# link to data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130646

################################################################################################
### Step 1: Read in the data, format it, and normalize it. #####################################
################################################################################################

#reading in data: 4 different files (only using the human scRNAseq data from the paper)
#downloaded and unzipped files in linux
all_files = ['GSM3746212_Muscle_1_Counts.csv', 'GSM3746213_Muscle_2_Counts.csv', 'GSM3746214_Muscle_3_Counts.csv', 'GSM3746215_Muscle_4_Counts.csv']

#set empty list to append the files tranformed to dataframes to
df_list = []

#read in files and combine them into one dataframe
for filename in all_files:
    df = pd.read_csv(filename, header= 0, index_col= 0)
    df_list.append(df)
final_df = pd.concat(df_list, axis = 1, join = 'outer')

#transform into Anndata object for scanpy
adata = sc.AnnData(final_df.T) 

#Normalization: they were not explicit in how they normalized their data so I used the 
#defualt method in scanpy (make it so that total counts in each column = 10000 and then log1p tranform dataframe)

#normalize data to sum/column = 10000
sc.pp.normalize_total(adata, target_sum=1e4) 
#log1p transformation
sc.pp.log1p(adata) 
#identifying highly variable genes based on seurat method with cutoff of 3000 
# (not sure what they did bc it was not clear, so this is my own choice)
sc.pp.highly_variable_genes(adata, flavor = 'seurat', n_top_genes = 3000) 
#plotting highly variable genes of normalized (left) and not normalized (right) data (see below)
sc.pl.highly_variable_genes(adata) 


################################################################################################
### Step 2: Perform PCA ########################################################################
################################################################################################

# perform PCA (they take only first 10 PCs as stated in the paper)
sc.tl.pca(adata) 
#plot PCA variance ratios for each PC (just to get a sense of the distribution) 
sc.pl.pca_variance_ratio(adata, log=True)

################################################################################################
### Step 3: Perform Clustering and T-SNE  ######################################################
################################################################################################

# this neighbors function does clustering 
# (chose 10 PCs based on paper and arbitratrily chose 10 neighbors, could do more or less)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10) 
# this function performs TSNE 
# (chose 10 PCs, but other parameters are somewhat arbitrary, but generally 1000 should be good to run T-SNE to stable state)
sc.tl.tsne(adata, n_pcs = 10, use_rep="X_pca", perplexity = 25, random_state = 0, learning_rate = 1000)
# leiden does clustering based on neighbors matrix calculations from sc.pp.neighbors
sc.tl.leiden(adata, key_added = "clusters", resolution = 0.3)
# plot the T-SNE result and clusters
sc.pl.tsne(adata, color = "clusters",  cmap="tab20")

# result: 11 clusters as in the paper
# next step: do some further analysis to recreate papers in the figure


#####################################################################################################################
### Step 4: We need to identify top expressed marker genes for each cluster to identify cell type for each group ####
#####################################################################################################################

# exploratory: we can plot (25) top differentially expressed genes for each group compared to the rest
sc.tl.rank_genes_groups(adata, 'clusters', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# we can also look at the top expression of each marker gene that they identified in the paper
#list of top marker genes for each of the 11 groups that were identified in the paper
markers = ['VWF', 'LUM', 'APOC1','RGS5','PRG4','CCL4','DARC','MYH11','CD3D','LST1','MS4A1']
#function to plot the T-SNE result with color scale as expression of each marker gene
sc.pl.tsne(adata, color=markers)

# we can also look at each cluster individually for each gene
for i in range(0,11):
    print("Cluster",i)
    sc.pl.tsne(adata[adata.obs.clusters.isin(["{}".format(i)]),:], color=markers, legend_loc='right margin', groups = ['clusters'])


##Analyzing Results  

# Based on their biological information and classification in the supplement for expected highest expressed genes in each cell type.
# We can conldue the cell type based on these plots:

# Cluster 0: Endothelial Cells

# Cluster 1: Lum + FAP Cells

# Cluster 2: Satellite Cells

# Cluster 3: Pericytes

# Cluster 4: Endothelial Cells

# Cluster 5: PRG4 + FAP Cells

# Cluster 6: NK Cells

# Cluster 7: PCV Endothelial Cells

# Cluster 8: Smooth Muscle Cells

# Cluster 9: T and B Cells

# Cluster 10: Myeloid Cells

#assigning clusters based on gene expression and what they classified in the paper
#they combined two clusters of endothelial cells in the paper
old_to_new = {
'0': 'Endothelial Cells',
'1': 'LUM+ FAP Cells',
'2': 'Satellite Cells',
'3': 'Pericytes',
'4': 'Endothelial Cells',
'5': 'PRG4 + FAP Cells',
'6': 'NK Cells',
'7': 'PCV Endothelial Cells',
'8': 'Smooth Muscle Cells',
'9': 'B Cells', #actually T cells and B cells but will change it after
'10':'Myeloid Cells'
}
adata.obs['annotation'] = (
adata.obs['clusters']
.map(old_to_new)
.astype('category')
)

#subset Anndata object for just T and B cells (labeld B Cells)
k = adata[adata.obs["annotation"] == "B Cells"] 
#according to the paper T cells are those which express CD3D while B cells do not 
# so subsetting by this below and returning list of genes that are T cells
gene1_cells = list(k[(k[:, ['CD3D']].X > 0),].obs.index)
#loop through anndata object and assign identified genes as B cells to new obs object in Anndata 'adata'
for gene in list(adata.obs.index):
    if gene in gene1_cells:
        adata.obs.loc[gene, "final_clusters"] = "T Cells"
    else:
        adata.obs.loc[gene, "final_clusters"] = adata[[gene], ["CD3D"]].obs["annotation"][0]


################################################################################################
### Step 5: Make Final Plots and Some Exploratory Plots  #######################################
################################################################################################

#plotting graph with new labels
with plt.rc_context({'figure.figsize':(15,10)}):
    sc.pl.tsne(adata, color="final_clusters", legend_loc='on data', title='', frameon=True, save='.pdf', palette = 'tab20',
    outline_color = ('black','green'), na_color='black', colorbar_loc='right')

#variation with connections between nodes
with plt.rc_context({'figure.figsize':(15,10)}):
    sc.pl.tsne(adata, color="final_clusters", legend_loc='on data', frameon=False, save='.pdf', palette = 'tab20',
    edges = True, title = "T-SNE Clustering by Cell Types")

# exploratory, try to categorize by top expressed gene
#read in dictionary of marker genes for each cell type (as provided by the paper)
marker_genes_dict= {'Endothelial Cells' : ['VWF','AQP1','ITGA6','A2M','IFI27','GNG11','CD74','CLDN5','CD36','BTNL9','ITGA1','FABP4','RNASE1','EGFL7','GPR116','AC011526.1','CDH5','HLA-B','RGCC','TMEM88','HLA-DRB1','ABLIM3','HLA-E','CAV1','CD93','TM4SF1','IGF2','RP11-782C8.5','PODXL','MGLL','B2M','EMCN','CA4','FABP5','FLT1','WARS','TMSB10','KIF25','BST2','RP11-417J8.6','SLC9A3R2','SHANK3','ESAM','TCF4','MLLT4','F8','HLA-A','CD300LG','PTPRB','ENG'],
'LUM+ FAP Cells' : ['COL15A1','APOD','LUM','DCN','ADH1B','MYOC','C3','CFD','ABCA8','PTGDS','GSN','SMOC2','DPT','PDGFRL','CXCL14','COL6A3','C1S','FBLN2','SERPINF1','MGP','FBLN1','CFH','COL1A2','PCOLCE','CILP','SCN7A','ABCA10','IGF1','LTBP4','COL6A2','FBLN5','IGFBP6','LAMA2','COL6A1','PDGFRA','ANGPTL1','SFRP4','PLAC9','C1R','LRP1','MMP2','CHRDL1','PODN','HSD11B1','FN1','ISLR','MFAP5','CCDC80','SERPING1','PGF'],
'Satellite Cells' : ['CXCL14','ALDH1A1','SPRY1','MEG3','CD63','APOC1','APOE','MYF5','TRDN','MYF6','MEST','MT1X','TUBB2B','CHRNA1','SPATS2L','FOSB','STAC3','FOS','PAX7','CRYAB','CSRP2','EGR1','IGFBP5','DLK1','JUNB','JUN','CADM2','PON2','PDLIM4','RASD1','PPP1R14B','SGCA','RXRG','CNN3','DAG1','DES','TPD52L1','ZFP36','DUSP1','HMGN2','NACA','CIRBP','PDLIM3','GNB2L1','MYC','HSP90AB1','RSL1D1','CYR61','EIF3L','CLU'],
'Pericytes' : ['TIMP3','COL3A1','IGFBP7','MALAT1','NDRG2','RGS5','NDUFA4L2','HIGD1B','ABCC9','NOTCH3','COX4I2','ASPN','KCNJ8','PDGFRB','MFGE8','STEAP4','CPE','MYO1B','TNFRSF21','CPM','EPS8','TESC','MYL9','ACTA2','FAM162B','NR2F2','SPARC','ITM2C','CYGB','PPP1R14A','CALD1','GJA4','EGFLAM','CALM2','C20orf27','ARHGDIB','SEPT7','SEPT4','ASAH1','THBS4','MCAM','FRMD3','ITGB1','ADIRF','MT2A','FAM213A','ISYNA1','TAGLN','TPM2','SELENBP1'],
'PRG4 + FAP Cells' : ['DCN','CFD','C1S','FBLN2','COL1A2','IGFBP6','PLAC9','MMP2','CHRDL1','MFAP5','CCDC80','NOVA1','TNXB','COL1A1','MFAP4','PI16','MGST1','VCAN','GPNMB','EFEMP1','ACKR3','FSTL1','TIMP2','IGFBP5','FBN1','PRG4','PCOLCE2','LINC01133','DEFB1','PLA2G2A','CD55','GFPT2','CD248','SEMA3C','UAP1','NTM','ADAMTS5','LOXL1','SCARA5','SPRR2B','FNDC1','CREB5','ALDH1A3','PTGIS','PROCR','FABP3','TRIO','TGFBR3','COL14A1','C17orf58'],
'NK Cells' : ['CCL4','NKG7','GNLY','GZMA','CCL3','CCL5','KLRB1','CTSW','FCER1G','HCST','KLRD1','GZMB','CST7','KLRF1','FCGR3A','GZMM','TYROBP','CORO1A','PTPRCAP','CD7','ITGB2','HOPX','PRF1','FGFBP2','IL2RB','GZMH','LCP1','CD247','PTPRC','PLAC8','DOK2','SAMD3','CD69','XCL2','PYHIN1','KLRC1','ALOX5AP','PLEK','CD53','MYOM2','CLIC3','LCK','IL2RG','C1orf162','RAC2','BIN2','UCP2','GPR65','C12orf75','TBC1D10C'],
'PCV Endothelial Cells' : ['AQP1','CD74','AC011526.1','HLA-DRB1','HLA-DRA','ELTD1','TM4SF18','HLA-DPA1','ECSCR','HLA-DMA','TMSB4X','HLA-DPB1','GIMAP4','FKBP1A','CLU','EEF1A1','DARC','CCL14','PLVAP','PLAT','HLA-DQA1','RAMP3','PERP','NOSTRIN','HLA-DQB1','JAM2','NUAK1','DUSP23','PRCP','CRHBP','PKP4','PIM3','SNHG7','TSPAN7','NPDC1','LIFR','FXYD5','NRN1','C10orf128','LMCD1','PLK2','LYST','TGFBR2','IL33','ETS2','TSHZ2','KCTD12','ERG','SLC2A3','LRRFIP1'],
'B Cells' : ['ETS1','GMFG','TSC22D3','JUNB','HNRNPA1','LDHB','EEF1A1','NPM1','EEF1B2','TMEM66','PABPC1','UBA52','ARHGDIB','SH3BGRL3','HCST','CORO1A','PTPRCAP','PTPRC','PLAC8','CD69','LCK','IL2RG','TBC1D10C','LIMD2','FYB','LAPTM5','CD48','EVL','CD37','ABRACL','CDC42SE1','CD52','C1orf56','STK4','BTG1','PFDN5','LTB','CD3D','IL7R','CD3E','SELL','CD3G','CD2','CXCR4','TRAF3IP3','ISG20','COTL1','RNASET2','NOSIP','COX7C'],
'T Cells' : ['ETS1','GMFG','TSC22D3','JUNB','HNRNPA1','LDHB','EEF1A1','NPM1','EEF1B2','TMEM66','PABPC1','UBA52','ARHGDIB','SH3BGRL3','HCST','CORO1A','PTPRCAP','PTPRC','PLAC8','CD69','LCK','IL2RG','TBC1D10C','LIMD2','FYB','LAPTM5','CD48','EVL','CD37','ABRACL','CDC42SE1','CD52','C1orf56','STK4','BTG1','PFDN5','LTB','CD3D','IL7R','CD3E','SELL','CD3G','CD2','CXCR4','TRAF3IP3','ISG20','COTL1','RNASET2','NOSIP','COX7C'],
'Smooth Muscle Cells' : ['COX4I2','MYL9','ACTA2','SEPT4','MCAM','ISYNA1','TAGLN','TPM2','SNCG','CRIP1','NTRK2','RERGL','MYH11','PLN','CASQ2','ACTG2','PTN','PHLDA2','ENTPD3','LBH','RCAN2','MTHFD2','CNN1','SORBS2','MYLK','MAP3K7CL','CCDC3','FRZB','RP11-332H18.4','SLC38A11','ITGA8','FBXL13','KCNMB1','PPP1R12B','GUCY1A3','BGN','SLIT3','NET1','RGS16','SOD3','LMOD1','PTPLA','EDNRA','GPRC5C','PIP5K1B','KLHL23','CKB','TBX2','RBPMS2','PRSS23'],
'Myeloid Cells' : ['FCER1G','TYROBP','C1orf162','CD48','COTL1','S100A8','LYZ','AIF1','S100A9','LST1','FCN1','CLEC7A','SERPINA1','IGSF6','CSTA','MS4A6A','IL1B','S100A12','HCK','SPI1','RP11-290F20.3','CFP','IFI30','MS4A7','MPEG1','LILRA5','C5AR1','RGS18','FPR1','CSF2RA','BCL2A1','FCGR1A','CYBB','AMICA1','IL8','RNASE6','RGS2','PILRA','APOBEC3A','CTSS','FAM26F','LY86','STXBP2','TNFSF13B','CD68','CD83','MNDA','PLAUR','TBXAS1','CD14']
}

#swtich to pandas dataframe since it is easier to parse
j = adata.to_df()
#get column name (gene) for each row maximum value
genes = j.idxmax(axis = 1)
#make dictionary into list that can be used in for loop
marker_list = sum(list(marker_genes_dict.values()),[])
#remove genes that are not marker genes 
for gene,ky,index in zip(genes,adata.obs.final_clusters,range(len(genes))):
    if gene not in marker_genes_dict[ky]:
        genes[index] = None

#assign new obs column for top genes
adata.obs["top_gene"]=genes

#plot top genes
with plt.rc_context({'figure.figsize':(10,5)}):
    sc.pl.tsne(adata, color="top_gene", legend_loc='right margin', title='', frameon=True, save='.pdf', palette = 'tab10')


# Analysis

# Clearly most of the top expressed genes are not marker genes, which kind of makes sense.  
# I need to do a more sensitive anaysis, to not just look at the top expressed gene, 
# but all the genes above some threshold that are some marker genes.  But that is highly complicated, 
# so I am just going to go with the first plot as my final plot.  Also getting a color legend and a 
# cluster label is very difficult to do with the standard scanpy commands.
    
#exploratory: top 3 differentially expressed for each cluster
sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, values_to_plot='logfoldchanges', min_logfoldchange=3, vmax=7, vmin=-7, cmap='bwr')

#recreating plot in Figure 2C
# the gene PERGL was not in the dataset even though it is in their plot so I used PERP instead...
# need to define order of cells
order = ['B Cells', 'Myeloid Cells',  'T Cells', 'Smooth Muscle Cells', 'PCV Endothelial Cells', 'NK Cells',
         'PRG4 + FAP Cells', 'Pericytes',
       'Satellite Cells', 'LUM+ FAP Cells', 'Endothelial Cells']
# plotting same markers as in figure in paper
markers2 = ['FABP4', 'IFI27', 'VWF', 'APOD', 'LUM', 'DCN','APOC1','APOE','MYF5','RGS5','NDUFA4L2','HIGD1B','FBN1',
'PRG4','PCOLCE2','CCL4','NKG7','GNLY','DARC','CCL14','PLVAP','PERP','MYH11','PLN','CD3D','IL7R','LTB','S100A8','LYZ','AIF1','MS4A1','TCL1A', 'FCRL1']
sc.pl.dotplot(adata, markers2, groupby='final_clusters', dendrogram=False, cmap = "Blues", 
mean_only_expressed = False,
categories_order = order,
use_raw = False,
standard_scale = 'var',
smallest_dot = 0)

#plotting figure in 2b from paper

#plot dataframe by cell type and color by sample.  customize titles and legend
sns.set_style("white")
ax = sns.histplot(data=dff, y="cell_type", hue="sample", multiple="stack", palette = sns.color_palette("tab10"))
ax.set(xlabel = 'Cell Count', ylabel = "")
plt.legend(title='Sample', loc='upper right', labels=['Sample 1', 'Sample 2','Sample 3','Sample 4'])
plt.show()