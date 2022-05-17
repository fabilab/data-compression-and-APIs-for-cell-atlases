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
- Enable gene search for the heatmap by timepoints.



#### 10/04/2022

- Create new dataset for the 3rd heatmap (combined heatmap with all cell types over different timepoints)
- Complete most part of the Interim report by 14/04/2022✅



#### 16/05/2022

- Heatmap by celltype (1st page): A equation is used to auto adjust the cell's size for the heatmap. For the scatter plot, since it looks a bit messy, I applied a log10 function to all values.✅
- Heatmap by dataset: after using the jinja template, the timepoint tooltip doesn't work anymore, need to be fix. size of the heatmaps are also needed to be adjusted.
- 3rd heatmap need to be done before Thursday. 