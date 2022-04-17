function HeatMap(result, html_element_id) {
        if (!result) {
            alert("Error:Input gene name is invalid, please make sure you type in the corrent gene names")
        } else {
            // x-axis: genes of interest
            var x_axis = Object.keys(result[Object.keys(result)[0]]);
            var ngenes =  Object.keys(result)[0].length;
            var graph_height = 300 + 44 * ngenes;
            // y-axis:41 cell types
            var y_axis = Object.keys(result);
            var data_content = [];
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
                autosize: true,
                height: graph_height,
                title: 'Heatmap of gene expression level in selected cell types',
                xaxis: {
                    title: 'Cell types',
                    automargin: true,
                    tickangle: 45,
                },
                yaxis: {
                    title: 'Genes',
                    automargin: true,
                    scaleanchor: 'x',
                    scaleratio: 1,
                    autorange: "reversed",
                },
            };
                
            Plotly.newPlot(
                document.getElementById(html_element_id),
                data,
                layout,
            ); 
        };
    } 
