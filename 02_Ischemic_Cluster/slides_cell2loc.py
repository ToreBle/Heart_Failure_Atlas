import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import gc
os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'
import cell2location
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
from cell2location.models import RegressionModel

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

sp_data_folder = /rwthfs/rz/cluster/home/noco0013/projects/mi_hearts/visium/

sample_name = sys.argv[1]

def read_and_qc(sample_name, path=sp_data_folder):
    r""" This function reads-in the data for all 10X spatial experiment into the anndata object 
    and calculates QC metrics.

    :param sample_name: Name of the sample
    :param path: path to data
    """
    adata = sc.read_visium(path + str(sample_name) + '/outs/',
                           count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.obs['sample'] = sample_name
    adata.var_names_make_unique(join='-')
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['ENSEMBL']
    adata.var.drop(columns='ENSEMBL', inplace=True)

    # Calculate QC metrics
    from scipy.sparse import csr_matrix
    adata.X = adata.X.toarray()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = csr_matrix(adata.X)
    adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    #adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var_names]
    adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'

    return adata

adata = read_and_qc(sample_name, path=sp_data_folder)

# mitochondria-encoded (MT) genes should be removed for spatial mapping
adata.obsm['mt'] = adata[:, adata.var['mt'].values].X.toarray()
adata = adata[:, ~adata.var['mt'].values]

with mpl.rc_context({'figure.figsize': [6,7],
                     'axes.facecolor': 'white'}):
        sc.pl.spatial(adata, img_key = "hires", cmap='magma',
                  library_id=list(adata.uns['spatial'].keys())[0],
                  color=['total_counts', 'n_genes_by_counts'], size=1,show=False, return_fig=True,
                  save="QC_"+sample_name+".pdf")

results_folder = '/hpcwork/noco0013/projects/novo_postdoc/data/MI_c2l/results/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures/'
run_name = f'{results_folder}/cell2location_map'

#For loading the model later
adata_file = f"{ref_run_name}sc.h5ad"
reference = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", reference)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in reference.varm.keys():
    inf_aver = reference.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in reference.uns['mod']['factor_names']]].copy()
else:
    inf_aver = reference.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in reference.uns['mod']['factor_names']]].copy()
inf_aver.columns = reference.uns['mod']['factor_names']

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata.var_names, inf_aver.index)
adata = adata[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata, batch_key="sample")

corr = inf_aver.corr()

# Generate a mask for the upper triangle
mask = np.triu(np.ones_like(corr, dtype=bool))

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(11, 9))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(corr, mask=mask, cmap=cmap, center=0.5,
            square=True, linewidths=.5, cbar_kws={'label': 'Pearson correlation - r', "shrink": .5}, )
plt.savefig("Pearson_correlation_reference_"+sample_name+".pdf", bbox_inches='tight')


# create and train the model
mod = cell2location.models.Cell2location(
    adata, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=5,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=200
)
mod.view_anndata_setup()

mod.train(max_epochs=15000,
          batch_size=None,
          train_size=1,
          use_gpu=True,
         )

# Export the estimated cell abundance (summary of the posterior distribution).
adata = mod.export_posterior(
    adata, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}/mod_"+sample_name, overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp_"+sample_name+".h5ad"
adata.write(adata_file)


adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
cell_types = list(adata.obs.columns[(adata.obs.columns.get_loc("Adipocytes")):adata.obs.columns.get_loc("VSMC")+1])

# add 5% quantile, representing confident cell abundance
adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
cell_types = list(adata.obs.columns[(adata.obs.columns.get_loc("Adipocytes")):adata.obs.columns.get_loc("VSMC")+1])

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'white',
                     "font.size" : 16,
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata, cmap='viridis',
                  # show first 8 cell types
                  color=cell_types,
                  ncols=1, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99', save=f"{run_name}/cell_type_deconvolution_"+sample_name+".pdf"
                 )
out = sc.get.obs_df(adata, keys=cell_types)
out.to_csv(f"{run_name}/c2l_meta_"+sample_name+".csv")