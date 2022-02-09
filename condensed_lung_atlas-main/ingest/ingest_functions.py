import scanpy as sc
import pandas as pd
import numpy as np
import os
import GEOparse as geo


def anndata_geo_metadata(adata, geo_id):
    '''
    attaches metadata pulled from Gene Omnibus to anndata objects
    Args:
        adata: anndata object to be annotated
        geo: str of GEO accession number e.g "GSEXXXXX"
    Output:
    anndata object with metadata info appended to .uns of anndata object'''
    geo_id = geo_id
    gse = geo.get_GEO(geo=geo_id) #object methods called in uns_dict
    #read in geo_metadata_dict from file
    uns_ls = ["title",
              "geo_accession",
              "last_update_date",
              "contributor",
              "overall_design",
              "citation"
                ]
    for x in uns_ls:
        try:
            adata.uns[x] = gse.get_metadata_attribute(x)
        except:
            continue
    os.remove(f'{geo_id}_family.soft.gz')
    return
def anndata_var_names_standardize(adata):
    gene_aliases = pd.read_csv('/home/carsten/alvira_bioinformatics/data/marker_gene_lists/mouse_gene.info',
                               sep='\t')
    alias_dict = {}
    for x in gene_aliases.index:
        key = gene_aliases.at[x, 'Symbol']
        if key in alias_dict.keys():
            alias_dict[key] = alias_dict[key] + gene_aliases.at[x, 'Synonyms'].split('|')
        else:
            alias_dict[key] = gene_aliases.at[x, 'Synonyms'].split('|')

    new_var_names = []
    for index, x in enumerate(adata.var_names):
        if x in alias_dict.keys():
            new_var_names.append(x)
        else:
            for k, v in alias_dict.items():
                if x in v:
                    if k in adata.var_names:
                         new_var_names.append(f'{k}_{x}')
                         print(f'{k} from {x}_k_exists')
                         break
                    else:
                        new_var_names.append(k)
                        print(f'{k} from {x}')
                        break
                else:
                    continue
            if index != len(new_var_names) - 1:
                print(f'{x} is unique')
                new_var_names.append(x)
    adata.var_names = new_var_names
    return
def anndata_ingest(adata, name,output ='/home/carsten/alvira_bioinformatics/data/external_datasets/ingested',geo = None,var_index = 'mgi_symbol',toc_fn = '/home/carsten/alvira_bioinformatics/data/external_datasets/ingested/anndata_ingest_toc.csv'):
    adata.uns['name'] = name
    fp = os.path.join(output,f"{name.strip()}.gz.h5ad")
    if geo != None:
        anndata_geo_metadata(adata, geo)
    ## metadata grab
    annot = sc.queries.biomart_annotations('mmusculus',
                                           ["ensembl_gene_id",
                                            'mgi_symbol',
                                            "chromosome_name",
                                            'gene_biotype',
                                            'description'
                                            ])
    print(annot.columns)
    for column in annot.columns:
        print(column)
        if column == var_index:
            continue
        else:
            annot[column] = annot.groupby([var_index])[column].transform(lambda x: '_'.join(x))
    annot.drop_duplicates(subset= var_index,
                          inplace = True)
    annot.set_index(var_index,
                    inplace = True)
    adata.var[annot.columns] = annot

    ## scanpy builtins
    sc.pp.calculate_qc_metrics(adata,
                               percent_top = None,
                               inplace = True,
                               log1p =False) #Does it in natural log which is non-intuative
    sc.pp.normalize_total(adata,
                          target_sum = 1e6,
                          key_added = 'coverage'
                          )
    sc.pp.log1p(adata,
                base = 2)
    # if var_index == 'mgi_symbol':
    #     anndata_var_names_standardize(adata)

    adata.write_h5ad(fp, compression = 'gzip')
    toc_dict = {'name':name,
                'cells': len(adata.obs_names),
                'log10_median_counts': np.log10(adata.obs['total_counts'].median()),
                'lineages': sorted(adata.obs['lineage'].unique().tolist()),
                'timepoints': sorted(adata.obs['Timepoint'].unique().tolist()),
                'treatments' : sorted(adata.obs['Treatment'].unique().tolist()),
                'file': fp
                }
    ## Add entry to table of context
    if os.path.exists(toc_fn) == True:
        toc = pd.read_csv(toc_fn,
                    header = 0,
                    index_col = 0,
                    )
        new = pd.DataFrame(pd.Series(toc_dict)).T
        new.set_index('name', inplace=True)
        toc = toc.append(new)
        toc = toc[~toc.index.duplicated(keep='last')]
        toc.to_csv(toc_fn)
    else:
        toc = pd.DataFrame(pd.Series(toc_dict)).T
        toc.set_index('name', inplace=True)
        toc.to_csv(toc_fn)

    return