var dataForPlotsDataset = {};

function generateDpDatasets(result_wrapper) {
    if (result_wrapper === "") {
        result_wrapper = dataForPlotsDataset['result_wrapper'];
    } else {
        dataForPlotsDataset['result_wrapper'] = result_wrapper;
    }
    
    let num = 1;
    for (index in Object.keys(result_wrapper['exp_avg'])) {
        let dataset = Object.keys(result_wrapper['exp_avg'])[index];
        let div_id = 'dp_dataset_' + num;
        dotPlotDataset(result_wrapper, div_id, dataset);
        num++;
    }
}

function dotPlotDataset(result_wrapper, html_element_id,dataset_name) {
    let useLog = dataForPlotsDataset['useLog'];
    
    const result_avg = result_wrapper['exp_avg'][dataset_name];
    const result_pro = result_wrapper['exp_pro'][dataset_name];
    
    let celltypeOrder = dataForPlotsDataset['celltypeOrder'];

    let celltypes
    if (!celltypeOrder) {
        celltypes = Object.keys(result_avg[Object.keys(result_avg)[0]]);
    } else {
        celltypes = result_wrapper['hierarchicalCelltypeOrder'][dataset_name];
    }
    
    // x-axis: celltypes
    let all_x = [];
    // y-axis: timepoints
    let all_y = [];

    //  expression level
    let color = [];
    // proportion exression
    let size = [];
    let desired_maximum_marker_size = 6.2;
    for (var i = 0; i < Object.keys(result_avg).length; i++) {
        let timepoint = Object.keys(result_avg)[i]
        
        for (var j = 0; j < celltypes.length; j++) {
            let exp = result_avg[timepoint][celltypes[j]];
            all_y.push(timepoint);
            all_x.push(celltypes[j])
            if (exp === -1) {
                exp = null;
            }
            if (useLog && exp !== null) {
                exp = Math.log10(exp + 0.5);
            }
            color.push(exp);
            size.push(result_pro[timepoint][celltypes[j]]*100);
        }
    }

    let nTimepoints = all_y.length;
    var data = {
        mode:'markers',
        marker: {
            color: color,
            size: size,
            sizeref: 2 * Math.max(...size) / (desired_maximum_marker_size**2),
            colorscale: 'YlGnBu',
            reversescale:true,
            colorbar: {},
        }
    };
    var layout = {
        title: 'belongs to dataset: '+ dataset_name,
        autosize: true,
        showlegend:false,
        xaxis: {
            title: '<b>Celltypes<b>',
            automargin: true,
            tickangle: 45,
        },
        yaxis: {
            title: '<b>Timepoints<b>',
            automargin: true,
        },
        width: 1030,
        height: 270+20*(Object.keys(result_avg))
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
            responsive: true,
    };
    
    if ($('#'+html_element_id).text() === "") {
        data['x'] = all_x;
        data['y'] = all_y;
        Plotly.newPlot(document.getElementById(html_element_id), [data],layout,tools);
    } else {
        data['x'] = [all_x];
        data['y'] = [all_y];
        Plotly.update(document.getElementById(html_element_id), data,layout, tools);
    }
} 