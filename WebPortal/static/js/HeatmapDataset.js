var dataForPlots = {};

function HeatmapDataset(result, html_element_id,dataset_name) {
    let useLog = dataForPlots['useLog'];
    
    
    
    if (!result) {
        alert("Error:Input gene name is invalid, please make sure you type in the corrent gene names")
    } else {
        // x-axis: celltypes
        let x_axis = Object.keys(result[Object.keys(result)[0]]);
        // y-axis: genes
        let y_axis = Object.keys(result);
        let data_content = [];
        for (var i = 0; i < Object.keys(result).length; i++) {
            cell_type = Object.keys(result)[i] // get the cell_type name as a string
            all_gene_expression = result[cell_type]         // find it from the dictionary as a key
            
            data_content.push(Object.values(all_gene_expression))
        }
        let ncelltypes = x_axis.length;
        let ngenes = y_axis.length;
        let heatmap_width = 1300;
        let heatmap_height = 270 + 41 * ngenes;
        var data = [
            {
                z: data_content,
                x: x_axis,
                y: y_axis,
                type: 'heatmap',
                hoverongaps: false
            }
            ];
        var layout = {
            title: 'belongs to dataset: '+ dataset_name,
            xaxis: {
                title: '<b>Celltypes<b>',
                automargin: true,
                tickangle: 45,
            },
            yaxis: {
                title: '<b>Timepoints<b>',
                automargin: true,
            },
            with: heatmap_width,
            height: heatmap_height,
        };
            
        Plotly.newPlot(document.getElementById(html_element_id), data,layout); 
    };
} 
