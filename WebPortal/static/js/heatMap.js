// Plot heatmap by celltype as a callback for the AJAX request
function HeatMap(result, html_element_id) {
        if (!result) {
            alert("Error:Input gene name is invalid, please make sure you type in the corrent gene names")
        } else {
            // x-axis: genes of interest
            var x_axis = Object.keys(result[Object.keys(result)[0]]);
            var y_axis = Object.keys(result);
            var ngenes =  y_axis.length;
            var graph_width = 1300;
            var graph_height = 370 + 26 * ngenes;
            // y-axis:41 cell types
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
                    hoverongaps: false,
                    colorscale: 'Reds',
                }
                ];
            var layout = {
                autosize: true,
                width: graph_width,
                height: graph_height,
                title: 'Heatmap of gene expression level in selected cell types',
                xaxis: {
                    title: 'Cell types',
                    automargin: true,
                    tickangle: 60,
                    scaleanchor: 'y',
                    scaleratio: 1,
                },
                yaxis: {
                    title: 'Genes',
                    automargin: true,
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
