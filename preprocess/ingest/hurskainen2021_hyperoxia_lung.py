
import os.path
import pandas as pd
import scanpy as sc
import sys
sys.path.insert(1, '/home/carsten/alvira_bioinformatics/lungsc_ck')
from ck_functions import anndata_ingest
source = '/home/carsten/alvira_bioinformatics/data/external_datasets/hurskainen2021_hyperoxia_lung/raw'
geo = 'GSE151974'

if __name__ == '__main__':
    metadata = pd.read_csv(os.path.join(source, "GSE151974_cell_metadata_postfilter.csv.gz"),
                           compression='gzip',
                           header=0,
                           index_col= 0,
                           sep=',',
                          )

    adata = sc.read_csv(os.path.join(source,"GSE151974_raw_umi_matrix_postfilter.csv.gz"),
                           first_column_names=True,
                           delimiter=',',
                        ).T
    adata.obs = metadata
    adata.obs.rename(columns = {'Age':'Timepoint',
                      'Oxygen': 'Treatment',
                      'CellType': 'cellSubtype'},
                     inplace = True)

    lineage_dict = {'Col13a1+ fibroblast':'mesenchymal',
                     'Mesothelial':'mesenchymal',
                     'Neut 1': 'immune',
                     'AT2 1': 'epithelial',
                     'B cell 1': 'immune',
                     'AT1': 'epithelial',
                     'AT2 2': 'epithelial',
                     'NK cell': 'immune',
                     'Cap-a':'endothelial',
                     'SMC':'mesenchymal',
                     'Alv Mf': 'immune',
                     'Col14a1+ fibroblast':'mesenchymal',
                     'Ciliated': 'epithelial',
                     'Int Mf': 'immune',
                     'Art':'endothelial',
                     'Cap':'endothelial',
                     'DC1': 'immune',
                     'CD8 T cell 1': 'immune',
                     'Lymph': 'immune',
                     'Pericyte 1':'mesenchymal',
                     'Neut 2': 'immune',
                     'Mono': 'immune',
                     'Club': 'epithelial',
                     'ILC2': 'immune',
                     'gd T cell': 'immune',
                     'CD4 T cell 1': 'immune',
                     'Pericyte 2':'mesenchymal',
                     'CD4 T cell 2': 'immune',
                     'Vein':'endothelial',
                     'DC2': 'immune',
                     'CD8 T cell 2': 'immune',
                     'B cell 2': 'immune',
                     'Mast Ba2': 'immune',
                     'Myofibroblast':'mesenchymal'
                    }
    adata.obs['lineage'] = [lineage_dict[x] for x in adata.obs['cellSubtype']]
    anndata_ingest(adata, 'hurskainen2021', geo = geo)







