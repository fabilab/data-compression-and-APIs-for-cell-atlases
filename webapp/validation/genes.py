# vim: fdm=indent
'''
author:     Fabio Zanini
date:       30/04/22
content:    Validate gene names et al.
'''
import pandas as pd

from models import gene_order


genes_idx = gene_order.index
gene_matrix = (pd.Series(gene_order.index)
                 .str
                 .split('', expand=True)
                 .iloc[:, 1:-1]
                 .fillna('')
                 .values
                 .astype('U1'))
gene_maxlen = gene_matrix.shape[1]


def convert_numbers_in_gene_name(gene):
    '''Convert numbers in gene name into digits'''
    from text_recognition.assets import numbers

    # It's typical for gene names to end with a number (e.g. Car4)
    endswithdigit = False
    for i in range(len(gene)):
        tail = gene[len(gene) - 1 - i:]
        if tail.isdigit():
            endswithdigit = True
            continue
        break

    if endswithdigit:
        # No gene name is purely a number, this should be fine
        tail = tail[1:]
        ndigit = int(tail)
    else:
        # Check if we can convert the end to a digit
        for ntext, ndigit in numbers[::-1]:
            if gene.endswith(ntext):
                gene = gene[:-len(ntext)]+str(ndigit)
                endswithdigit = True
                break

    # If there are no convertible-to-digits at the end, we are done
    # FIXME: be more flexible than this, of course
    if not endswithdigit:
        return gene

    # If we found or converted an end digit, we should look for
    # internal digit-likes, e.g. Col1a1
    sfx = len(str(ndigit)) + 1
    for ntext, ndigit in numbers[::-1]:
        if gene[:-sfx].endswith(ntext):
            gene = gene[:-sfx-len(ntext)] + str(ndigit) + gene[-sfx:]
            break

    print(gene)
    return gene


def validate_correct_gene(gene, species='mouse'):
    '''Validate and correct misspellings for a single gene name'''

    # Murine genes are capitalized
    # Human/monkey genes are upper
    if species == 'mouse':
        gene = gene.capitalize()
    else:
        gene = gene.upper()

    if gene in genes_idx:
        return gene

    gene = convert_numbers_in_gene_name(gene)

    # Not found in whitelist, try to correct
    gene_array = list(gene)[:gene_maxlen]
    hamming = (gene_array != gene_matrix[:, :len(gene_array)]).sum(axis=1)

    # If the uncorrected gene is shorter (e.g. Col1) it can be a perfect match
    # for multiple whitelist genes (e.g. Col1a1, Col1a2), then ask for confirmation
    idx = (hamming == 0).nonzero()[0]
    # TODO: build mismatch scoring table based on natural English (e.g. p-t)
    if len(idx) == 0:
        idx = (hamming == 1).nonzero()[0]
    # If there is only one (partial) perfect match, take it. If multiple
    # perfect matches, ask. If no perfect and one imperfect match, take it.
    # If multiple imperfect or no imperfect, ask.
    if len(idx) == 1:
        gene_closest = genes_idx[idx[0]]
        print(gene, gene_closest)
        return gene_closest
    else:
        # At least one gene was not understood, ask for a written confirmation
        return None


def validate_correct_genestr(genestr, species='mouse'):
    '''Validate gene names and correct misspellings if possible'''

    # TODO: check punctuation more accurately
    genes = genestr.strip(' ').replace('.', ',').replace(';', ',').split(',')

    # Validate
    genesv = []
    for gene in genes:
        genev = validate_correct_gene(gene, species=species)
        if genev is None:
            return None
        genesv.append(genev)

    genestr = ','.join(genesv)
    return genestr
