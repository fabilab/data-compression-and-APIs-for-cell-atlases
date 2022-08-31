var dataForPlotsDataset = {};

function generateHmDatasets(result_wrapper) {
    if (result_wrapper === "") {
        result_wrapper = dataForPlotsDataset['result_wrapper'];
    } else {
        dataForPlotsDataset['result_wrapper'] = result_wrapper;
    }
    
    let num = 1;
    for (index in Object.keys(result_wrapper['result'])) {
        let dataset = Object.keys(result_wrapper['result'])[index];
        let div_id = 'dataset_' + num;
        HeatmapDataset(result_wrapper, div_id, dataset);
        num++;
    }
}

function HeatmapDataset(result_wrapper, html_element_id,dataset_name) {
    let useLog = dataForPlotsDataset['useLog'];
    
    const result = result_wrapper['result'][dataset_name];
    
    let celltypes;
    let celltypeOrder = dataForPlotsDataset['celltypeOrder'];

    if (!celltypeOrder) {
        celltypes = Object.keys(result[Object.keys(result)[0]]);
    } else {
        celltypes = result_wrapper['hierarchicalCelltypeOrder'][dataset_name];
    }
    
    if (!result) {
        alert("Error:Input gene name is invalid, please make sure you type in the corrent gene names")
    } else {
        // x-axis: celltypes
        let x_axis = celltypes;
        // y-axis: timepoints
        let y_axis = Object.keys(result);
        let data_content = [];
        for (var i = 0; i < Object.keys(result).length; i++) {
            let timepoint = Object.keys(result)[i] // get the timepoint name as a string
            let all_gene_expression = [];        // find it from the dictionary as a key
            for (var j = 0; j < celltypes.length; j++) {
                let exp = result[timepoint][celltypes[j]];
                if (exp === -1) {
                    exp = null;
                }
                if (useLog && exp !== null) {
                    exp = Math.log10(exp + 0.5);
                }
                all_gene_expression.push(exp);
            }
            data_content.push(all_gene_expression);
        }

        let ncelltypes = x_axis.length;
        let nTimepoints = y_axis.length;
        let heatmap_width = 1300;
        let heatmap_height = 270 + 20 * nTimepoints;
        var data = {
            type: 'heatmap',
            hoverongaps: false,
            colorscale: 'Reds',
        };
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
        var tools = {
            modeBarButtonsToAdd: [{
                name: 'Download plot as an SVG',
                icon: Plotly.Icons.camera,
                click: function(gd) {
                  Plotly.downloadImage(gd, {format: 'svg'})
                }
              }],
              editable:true,
        };
        
        if ($('#'+html_element_id).text() === "") {
            data['z'] = data_content;
            data['x'] = x_axis;
            data['y'] = y_axis;
            Plotly.newPlot(document.getElementById(html_element_id), [data],layout,tools);
        } else {
            data['z'] = [data_content];
            data['x'] = [x_axis];
            data['y'] = [y_axis];
            Plotly.update(document.getElementById(html_element_id), data,tools);
        }
    };
} 

function AjaxDatasets() {
  
  var gene_name = $('#singleGene').val();
  // sent gene names to the API
  $.ajax({
    type:'GET',
    url:'http://127.0.0.1:5000/data_datasets',
    data: "gene=" + gene_name,
    success: function(result) {
        $("#displayPlotDataset").empty();
        generateHmDatasets(result);
    },
    error: function (e) {
        Swal.fire('Invalid input','please make sure you type in the correct gene name','error');

    }
    });
  }

