# Compressed atlas of the murine lung
This web application exemplifies the idea of a "compressed cell atlas", i.e. a nonredundant distillation of one or more single cell omics data sets (three in this case).

The architecture of the compressed atlas is the following:
- A RESTful APIs to request the compressed data (e.g. `/data/gene_names=Car4`).
- A set of interative plots (mostly heatmaps or variations on the theme, e.g. dot plots) to visuaise the compressed data.
- A voice control system enabling a natural language interface to the visualisations.

At this time, this application is pre-alpha, so the API changes all the time. If you are interested in how it works, write me an email at fabio _DOT_ zanini _AT_ unsw _DOT_ edu _DOT_ au.
