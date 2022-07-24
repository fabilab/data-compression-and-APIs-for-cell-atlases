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

##### * Current Issue:

![Screen Shot 2022-06-27 at 9.01.15 pm](/Users/yingxu/Desktop/Screen Shot 2022-06-27 at 9.01.15 pm.png)

It is validated that the current algorithm does successfully selected the marker genes. However,



##### Algorithm Validation (Keyi's suggestion)

1) Randomly pick 5 celltypes from the list

2) plot the marker genes heatmap

3) pick the genes that is highly expressed (red color)

4) find papers to confirm whether they are real marker genes for that cell type

- if yes, algorithm works
- if not, possibily be a novel marker gene, need more experiment to validate.

5) why single cell ---> discover new cell types and new genes



##### suggestion from Afrina, the bioinformatician:

1) Search paper for all the house keeping genes in Neonatal Lung of mice

2) Remove all the house keeping genes from the original dataset, to reduce the dataset size in the first place.

3) When apply the current algorithm, also consider the P-value. (some genes has high expression, but may not be significant)

--> this page help discover not just possible novel marker genes, but also novel house keeping genes.

--> visualisation of data, and provide a new way for exploring dataset.



- Adding an new tab on the top: order by: expression level of the selected cell (highest to lowest)





##### marker genes from Fabio's paper:

https://www.biorxiv.org/content/10.1101/2021.04.27.441649v1







Package:

- Go to where setup.py is, cd .. (back one step)
- pip install cell_atlas_portal_2022 (for new user), pip install -e cell_atlas_portal_2022 (for update)
- Go to the directory that contains the data file (h5)
- python3
- from package import class
  - from cell_atlas_portal_2022 import h5Reader
- h5_file = h5Reader('./condensed_lung_atlas_in_cpm.h5')



#### Week7: Improve webpage usability and adding user guide and package documentation to the page

- Combine all heatmap display into one html page: data explore

- Switching tabs to generate different type of heatmap.



#### Week8 Presentation prepare

- create a new page: All data 

  - (provide an option for the user to look at the entire dataset as the table)
  - user should execute this page with caution, since the dataset is large
  - provide sorting button to the page (sort by gene names (alphabetically), sort by expression level within a cell type)

- Create a new page:

  - Explain the compression alogirthm

  



background 2-3 slides

- goal of the compression, how and why, 
- the need for visualisation
- the need for a simple to use interface (can show some example of currently available once but don't have too much detail here, the main point is to show that the current ones are slow and a bit difficult to use)

design (DO NOT TALK ABOUT USER ACTION, you do that in demo)

- what did you use for frontend
  - jquery 
    - lightweight, fast simple, no dependencies
  - Plotly js
    - powerful plotting package 
- what did you use for backend
  - flask api 
  - Pandas for data manipulation 
  - did you use scanpy? 
- what types of plots and why
  - scatter plot
    - investigate linear relationship between expression values in cell types of two genes
  - heatmap
    - shows averaged expr values of multiple genes in multiple cell types simultaneously 
    - Color => expression value
    - support Clustering of cell types
  - bubble map
    - just like heatmap
    - shows one more dimension: the proportion of cells expressing a specific gene
  - Why?
    - yes, there are more complicated ways to visualise sc rna data like t-SNE plots, but that may not be acheievable due to compressed data, which is done for easy sharing of experiment results
    - this portal aims to support early investigation of data / target users may not have prior experience with bioinformatics  -> therefore we use simple diagram
- what questions are the plots trying to answer
  - I want to know the expression information (level, proportion) of a list of genes in different cell types  (by cell type)
  - Given a gene, what is its expression level in different time points (development progression of this gene)
  - What are the outlier genes for a given cell types threshold (mention later in 'future plan section' saying that you will improve this later by doing ...) 

demo

- Demo the portal
- demo the package
  - how to install
  - how to import
  - mention that you can can use your own compressed dataset
  - show the user guide html

future plan

- your future plan
  - Improvements (mainly for marker gene page I guess)
  - documentation (better documentation)
  - put to public domain? both portal and package
  - test the portal using different cell atlas from other organs 
- reflection
  - compare to your industry training experience
  - compare to your uni study experience
  - Not just to get things down, but to think of better or smarter way to achieve it. 
    - e.g: speed up the program, minimise duplicate things. 
  - Learn to use different tools
  - How does the result make sense biologically (
    - e.g: when doing the marker genes page, I select the outlier genes of each celltype and called them 'marker genes', sometimes we can't trust the result, we must validate it using different resources. e.g: % of the genes in the list are actual marker genes that has been proofed by literature) ---> am I doing this correctly? 
  - Always seek for other people (professional's help)
  - At the begining, I was always relied on my supervisor and aks him what I should do next. If he told me to think of it myself, I would feel sad and don't know what to do. But now I realise that a research project should be something that what I want to do, and what I want to achieve through it. It is more like an opportunity...blabla
  - Always start things from simple and make it work. 





QR code (scan and use) -----> for poster in Thesis C



Cell type oder (change)



dataset(link to paper)

 p7 -1, p7-2 (timepoint unified), hover over can not be displayed

...Rik : not name () --> go to ncbi ---> unknot ---> find its functional (most of the naming is based on its functionality and families)



case insensitive change is good 

mouse: capitalise

Human: 



Celltype: link (or description)



marker genes (rank) --> gives a filter (top 10/top30 high expression).



Unknot genes ()RIK ----> obvious marker ---> meaning of research

H5--->AnnData (convertor)



dot plot title



install package css



Click on questions and redirected to the corresponded page.





























