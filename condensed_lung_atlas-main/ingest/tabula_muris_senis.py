import sys
import scanpy as sc
sys.path.insert(1, '/home/carsten/alvira_bioinformatics/condensed_atlas/ingest')
from ingest_funtctions import anndata_ingest
facs_source = '/home/carsten/alvira_bioinformatics/data/external_datasets/tabula_muris_senis/raw_data/tabula-muris-senis-facs-official-raw-obj.h5ad'
drop_source = '/home/carsten/alvira_bioinformatics/data/external_datasets/tabula_muris_senis/raw_data/tabula-muris-senis-droplet-official-raw-obj.h5ad'
data_source = 'https://figshare.com/articles/dataset/Processed_files_to_use_with_scanpy_/8273102'
lin_dict = {'Adventitial Fibroblast':('mesenchymal', 'Adventitial Fibroblast'),
              'Airway Smooth Muscle': ('mesenchymal', 'Airway Smooth Muscle') ,
              'Alveolar Epithelial Type 2':('epithelial','Alveolar Epithelial Type 2'),
              'Alveolar Fibroblast': ('mesenchymal', 'Alveolar Fibroblast'),
              'Alveolar Macrophage': ('immune','Alveolar Macrophage'),
              'Artery': ('endothelial', 'Artery'),
              'B': ('immune', 'B'),
              'Basophil': ('immune','Basophil'),
              'CD4+ T': ('immune','T'),
              'CD8+ T': ('immune','T'),
              'Capillary': ('endothelial', 'General Capillary'),
              'Capillary Aerocyte': ('endothelial', 'Aerocyte'),
              'Ccr7+ Dendritic': ('immune', 'Dendritic'),
              'Ciliated':('epithelial', 'Ciliated'),
              'Classical Monocyte': ('immune','Classical Monocyte'),
              'Club':('epithelial', 'Club'),
              'Intermediate Monocyte': ('immune', 'Intermediate Monocyte'),
              'Interstitial Macrophage': ('immune', 'Interstitial Macrophage'),
              'Ly6g5b+ T': ('immune', 'T'),
              'Lympatic': ('endothelial','Lymphatic'),
              'Myeloid Dendritic Type 1':( 'immune', 'Myeloid Dendritic'),
              'Myeloid Dendritic Type 2': ('immune', 'Myeloid Dendritic'),
              'Myofibroblast':('mesenchymal', 'Myofibroblast'),
              'Natural Killer': ('immune', 'Natural Killer'),
              'Natural Killer T':( 'immune', 'Natural Killer T'),
              'Neuroendocrine':('epithelial', 'Neuroendocrine'),
              'Neutrophil':('immune', 'Neutrophil'),
              'Nonclassical Monocyte': ('immune', 'Nonclassical Monocyte'),
              'Pericyte':('mesenchymal', 'Pericyte'),
              'Plasma': ('immune', 'Plasma'),
              'Plasmacytoid Dendritic': ('immune', 'Plasmacytoid Dendritic'),
              'Proliferating Alveolar Macrophage': ('immune', 'Alveolar Macrophage'),
              'Proliferating Classical Monocyte': ('immune', 'Classical Monocyte'),
              'Proliferating Dendritic': ('immune', 'Dendritic'),
              'Proliferating NK': ('immune', 'Natural Killer'),
              'Proliferating T': ('immune', 'T'),
              'Regulatory T': ('immune', 'Regularory T'),
              'Vein': ('endothelial', 'Vein'),
              'Zbtb32+ B': ('immune', 'B'),
             'Alox5+ Lymphocyte':('immune', 'Lymphocyte'),
             'Alveolar Epithelial Type 1':('epithelial',  'Alveolar Epithelial Type 1'),
             'Alveolar Epithelial Type 2 Cell':('epithelial','Alveolar Epithelial Type 2'),
             'B Cell':('immune','B'),
             'Basal':('epithelial','Basal'),
             'Capillary Type 1 Cell':('endothelial', 'General Capillary'),
             'Capillary Type 2 Cell':('endothelial', 'Aerocyte'),
             'Ciliated Cell':('epithelial', 'Ciliated'),
             'Club Cell':('epithelial', 'Club'),
             'Dendritic Cell':('immune', 'Dendritic'),
             'Lymphatic Cell':('endothelial', 'Lymphatic'),
             'Natural Killer T Cell':('immune', 'Natural Killer T'),
             'Proliferating Immune':('immune', 'Proliferating Immune'),
             'T cell':('immune', 'T'),
             }

if __name__ == '__main__':
    for data in [('facs',facs_source), ('drop',drop_source)]:
        adata = sc.read(data[1])
        adata_lung = adata[adata.obs['tissue'] == 'Lung', :].copy()
        #standardizing obs categories, cleaning up so that categories are just raw annotations
        adata_lung.obs['cellSubtype'] = adata_lung.obs['free_annotation']
        adata_lung.obs['Mousename'] = adata_lung.obs['mouse.id']
        adata_lung.obs['Gender'] = adata_lung.obs['sex']
        adata_lung.obs['Gender'] = ['F' if x == 'female' else 'M' for x in adata_lung.obs['Gender']]
        adata_lung.obs['Timepoint'] = adata_lung.obs['age']
        adata_lung.obs['Treatment'] = 'normal'
        adata_lung.obs['lineage'] = [lin_dict[x][0] for x in adata_lung.obs['cellSubtype']]
        adata_lung.obs['lineage'] = [lin_dict[x][0] for x in adata_lung.obs['cellSubtype']]

        #obs that were renamed
        del adata_lung.obs['sex']
        del adata_lung.obs['mouse.id']
        del adata_lung.obs['age']
        del adata_lung.obs['n_genes']

        adata_lung.uns['data_source'] = data_source

        anndata_ingest(adata_lung, f'TMS2020_{data[0]}')







