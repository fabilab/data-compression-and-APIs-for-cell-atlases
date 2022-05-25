var dataForPlots = {};

function HeatmapCelltype(result_wrapper, html_element_id) {
        // const start = performance.now();
        // console.log("start plotting the heatmap");
        // check the button id active
        let useLog = dataForPlots['useLog'];

        if (result_wrapper === "") {
            result_wrapper = dataForPlots['result_wrapper'];
        } else {
            dataForPlots['result_wrapper'] = result_wrapper;
        }

        if (html_element_id === "") {
            html_element_id = "h5_data_plot";
        }

        let result = result_wrapper['result'];
        let celltypes;
        let celltypeOrder = dataForPlots['celltypeOrder'];

        if (!celltypeOrder) {
            celltypes = Object.keys(result[Object.keys(result)[0]]);
        } else {
            celltypes = result_wrapper['hierarchicalCelltypeOrder'];
        }

        if (!result) {
            alert("Error:Input gene name is invalid, please make sure you type in the corrent gene names")
        } else {
            // x-axis: cell types
            let x_axis = celltypes;
            // y-axis: genes
            let y_axis = Object.keys(result);
            
            let ngenes = y_axis.length;
            let ncelltypes = x_axis.length;

            let heatmap_width = 1300;
            let heatmap_height = 270 + 39 * ngenes;

            let data_content = [];
            for (var i = 0; i < ngenes; i++) {
                let each_gene_data = [];
                for (var j = 0; j < ncelltypes; j++) {
                    exp = result[y_axis[i]][x_axis[j]];
                    if (useLog) {
                        exp = Math.log10(exp + 0.5);
                    }
                    each_gene_data.push(exp);
                }
                data_content.push(each_gene_data);
            }

            var data = {
                type: 'heatmap',
                hoverongaps: false
            };
            var layout = {
                autosize: true, 
                title: 'Heatmap of gene expression level in selected cell types',
                xaxis: {
                    title: '<b>Cell types<b>',
                    automargin: true,
                    tickangle: 45
                },
                yaxis: {
                    title: '<b>Genes<b>',
                    automargin: true
                },
                with: heatmap_width,
                height: heatmap_height,
            };
            
            if ($('#'+html_element_id).text() === "") {
                data['z'] = data_content;
                data['x'] = x_axis;
                data['y'] = y_axis;
                Plotly.newPlot(document.getElementById(html_element_id), [data],layout);
            } else {
                data['z'] = [data_content];
                data['x'] = [x_axis];
                data['y'] = [y_axis];
                Plotly.update(document.getElementById(html_element_id), data);
            }
        };
        const duration = performance.now() - start;
        console.log("heatmap done:" + duration/1000 + " second");
    } 
