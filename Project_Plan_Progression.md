#### Friday 18/02/2022

- Fix axis name of the heat map (it was hidden by the cell type names atm) ✅
- Add a rest button on the page, once clicked, the heatmap automatically changed from 'selected genes' back to the default plot. ✅
- Add a toolbar at the top left, with 3 button : 'default', 'Log10', 'Hierachical Clustering'. Enable the user to visualise the plot with different normailisation method and unsupervise learning. ✅



#### 05/03/2022

- Create a new page that displays 3 heat map at different timepoints:
  -  (data, dataset,dataset_timepoint) 
- Try to move everything to react (using a new git branch)
- Start Intrim Report:
  - 2-3 pages of Literature Review (single cell paper given by Givanna)
  - 1 page of Deployment (AWS or Google)
  - screen shot of current web portal 
  - Moke diagrams of future pages (what we also want to achieve)
  - Schematic view of the web application (how different entities connect with each other)
  - Screen shot of my codes



#### 11/03/2022

- Create a new page: heatmap by timepoints ✅
- adding a ? Info button next to each heatmap, display a small user menu when clicked.
- Enable gene search for the heatmap by timepoints.✅



#### 10/04/2022

- Create new dataset for the 3rd heatmap (combined heatmap with all cell types over different timepoints)✅
- Complete most part of the Interim report by 14/04/2022✅



#### 16/05/2022

- Heatmap by celltype (1st page): A equation is used to auto adjust the cell's size for the heatmap. For the scatter plot, since it looks a bit messy, I applied a log10 function to all values.✅
- Heatmap by dataset: after using the jinja template, the timepoint tooltip doesn't work anymore, need to be fix. size of the heatmaps are also needed to be adjusted.
- 3rd heatmap need to be done before Thursday. ✅


#### 22/05/2022
- Speed up the webpage (read_files function and heatmap loading time)✅



#### June to Sep Plan (14 weeks)

Week1: find a tutorial for building python package to handle compressed data 

(access the data / visualise the data)

* (high priority) access the data: functions that enable access of specific gene / celltype / timepoint. -----> send request to the API ----> JSON --->return as dataframe
* push it to github
* Test the function
* User guide

https://realpython.com/python-modules-packages/

https://www.tutorialsteacher.com/python/python-package

https://www.freecodecamp.org/news/build-your-first-python-package/

fixing small features: 

auto-completion of the seach bar (or case-insensitive) , gives ncbi or geneCards (GoTerm) link gene name.

Week2: Make a new front end page (show heatmaps of marker gene)

Week3: start to build interface for reading h5 file 

- import ying_package
- Data = ying_package.load_compressed_cell_atlas("filename.h5")

Week4: replace heatmap with dot plot, size of dot ~ expression level of genes?

###### ** Week5: build package + learn Carsten's algo (goal of the compression, how and why, challenges vs CELLXGENE)

Week6: improve frontend (responsive page: size of the page change when the user's device is mobile)

Week7: build package + learn Carsten's algo + presentation prep

Week8: improve frontend + presentation prep (demo of website and ppt).  ------> slides ready

Week9: presentation prep (demo of website and ppt)

##### * week 10: Thesis B presentation 



#### Plan by priority:

##### High:

build package + learn Carsten's algo

Presentations (making slides and script)

- goal of the compression, how and why, challenges vs CELLXGENE
- show people how to import, run and use the package ----> should show the returned data
- Demo the webportal

Q:

The needs for data compression

Long term goal (continue?)

- Applied to different data (different organs)

If you could start the project, what would you improve

why gene expression in these celltypes ?

why scRNA-seq but not bulk RNA seq ?

Target users 

-  Biologist: can not code, but understand the biology and interested in the result, 
- Bioinformaticians: who can code, but want to look at a simplied version of the data



Plan for thesis C:

- pushing package to pipy
- polishing / user friendly
- documentation of the package
- Deployment



Skills/Gaining:

- making plans/ Multi-taskings / project management
- Coding (font,backend,data-analysis)





#### Marker Genes page: (23/Jun/2022)

##### Algorithm Validation (Keyi's suggestion)

1) Randomly pick 5 celltypes from the list

2) plot the marker genes heatmap

3) pick the genes that is highly expressed (red color)

4) find papers to confirm whether they are real marker genes for that cell type

- if yes, algorithm works
- if not, possibily be a novel marker gene, need more experiment to validate.

5) why single cell ---> discover new cell types and new genes



suggestion from Afrina, the bioinformatician:

1) Search paper for all the house keeping genes in Neonatal Lung of mice

2) Remove all the house keeping genes from the original dataset, to reduce the dataset size in the first place.

3) When apply the current algorithm, also consider the P-value. (some genes has high expression, but may not be significant)

--> this page help discover not just possible novel marker genes, but also novel house keeping genes.

--> visualisation of data, and provide a new way for exploring dataset.





- Adding an new tab on the top: order by: expression level of the selected cell (highest to lowest)

##### House keeping genes



##### marker genes from Fabio's paper:

https://www.biorxiv.org/content/10.1101/2021.04.27.441649v1



