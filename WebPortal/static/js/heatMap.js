function HeatMap(result, html_element_id) {
        if (!result) {
            alert("Error:Input gene name is invalid, please make sure you type in the corrent gene names")
        } else {
            // x-axis: cell types
            let x_axis = Object.keys(result[Object.keys(result)[0]]);
            // y-axis: genes
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
                title: 'Heatmap of gene expression level in selected cell types',
                xaxis: {
                    title: '<b>Cell types<b>',
                    automargin: true,
                    tickangle: 45,
                },
                yaxis: {
                    title: '<b>Genes<b>',
                    automargin: true
                },
                with: 700,
                height: 700,
            };
                
            Plotly.newPlot(document.getElementById(html_element_id), data,layout); 
        };
    } 
