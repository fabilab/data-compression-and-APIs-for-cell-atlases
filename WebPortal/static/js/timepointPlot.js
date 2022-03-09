function HeatMapTimepoint(result, html_element_id,dataset_name) {
    if (!result) {
        alert("Error:Input gene name is invalid, please make sure you type in the corrent gene names")
    } else {
        // x-axis: 5 genes of interest
        let x_axis = Object.keys(result[Object.keys(result)[0]]);
        // y-axis:41 cell types
        let y_axis = Object.keys(result);
        let data_content = [];
        for (var i = 0; i < Object.keys(result).length; i++) {
            cell_type = Object.keys(result)[i] // get the cell_type name as a string
            all_gene_expression = result[cell_type]         // find it from the dictionary as a key
            
            data_content.push(Object.values(all_gene_expression))
        }
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
            title: 'dataset: '+ dataset_name,
            xaxis: {
                title: 'Cell types',
                automargin: true,
                tickangle: 45,
            },
            yaxis: {
                title: 'Timepoint',
                automargin: true,
            },
        };
            
        Plotly.newPlot(document.getElementById(html_element_id), data,layout); 
    };
} 
