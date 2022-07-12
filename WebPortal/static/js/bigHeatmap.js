dataForPlotsUnified = {}

function plotHeatmapUnified(result,html_element_id,gene_name) {
    // flags passed by events.js
    let useLog = dataForPlotsUnified['useLog'];
    let celltypeOrder = dataForPlotsUnified['celltypeOrder'];
    
    if (result === "") {
        result = dataForPlotsUnified['result'];
    } else {
        dataForPlotsUnified['result'] = result;
    }

    if (html_element_id === "") {
        html_element_id = "displayPlotUnified";
    }

    let celltypes;
    if (!celltypeOrder) {
        celltypes = result['cell_type'];
    } else {
        celltypes = result['hierarchicalCelltypeOrder'];
    }

    const expression = result['expression'];
    
    // x-axis: celltypes
    let x_axis = celltypes;
    //y-axis: timepoint
    let y_axis = [];
    //z: expression
    let data_content = [];
    //dataset
    let y_axis_dataset = [];
    // timepoint that with duplicate:
    let y_axis_time = [];
    // plot only the timepoint on Y-axis
    for (var i = 0; i < result['dataset_timepoint'].length; i++) {
        // dataset timepoint combination: e.g TMS_24m
        let dataset_timepoint = result['dataset_timepoint'][i];
        let time  = dataset_timepoint.split("_")[1];
        let dataset = dataset_timepoint.split("_")[0];
        if (y_axis.includes(time)) {
            y_axis.push(" ");
        } else {
            y_axis.push(time);
        }
        y_axis_time.push(time);
        y_axis_dataset.push(dataset);
        
        let expression_in_all_celltypes = [];
        for (var k = 0; k < celltypes.length; k++) {
            expression_in_all_celltypes.push(expression[dataset_timepoint][celltypes[k]]);
        }

        if (useLog) {
            for (var j = 0; j < expression_in_all_celltypes.length; j++) {
                if (expression_in_all_celltypes[j] !== -1) {
                    expression_in_all_celltypes[j] = Math.log10(expression_in_all_celltypes[j] + 0.5);
                }
            }
        }
        
        data_content.push(expression_in_all_celltypes);
    }
    // each with this format
    // Celltype: {ct}, Expression: {exp}, Dataset: {ds}, Timepoint: {tp}, 
    let hover_text = [];
    for (var i = 0; i < result['dataset_timepoint'].length; i++) {
        let temp = [];
        for (var j = 0; j < result['cell_type'].length; j++) {
            let dt = result['dataset_timepoint'][i];
            let ct = result['cell_type'][j];
            
            let exp = expression[dt][ct];
            // "ACZ_P21"
            let ds = dt.split("_")[0];
            let tp = dt.split("_")[1];
            temp.push('Celltype: '+ct+', Expression: '+exp+', Dataset: '+ds+', Timepoint: '+tp);
        }
        hover_text.push(temp);
    }

    let ncelltypes = x_axis.length;
    let ngenes = y_axis.length;
    let heatmap_width = 1300;
    let heatmap_height = 270 + 41 * ngenes;
    var data = {
        type: 'heatmap',
        text: hover_text,
        hoverinfo: 'text',
        colorscale: 'Reds',
    };
    console.log(result['gene']);
    var layout = {
        title: 'Expression of ' + result['gene'] + ' gene over time',
        xaxis: {
            title: '<b>Celltypes<b>',
            automargin: true,
            tickangle: 45,
        },
        yaxis: {
            title: '<b>Dataset Timepoints<b>',
            automargin: true,
        },
        with: heatmap_width,
        height: heatmap_height,
    };
    if ($('#'+html_element_id).text() === "") {
        data['z'] = data_content;
        data['x'] = x_axis;
        data['y'] = y_axis;
        Plotly.newPlot(document.getElementById(html_element_id), [data], layout);
    } else {
        data['z'] = [data_content];
        data['x'] = [x_axis];
        data['y'] = [y_axis];

        Plotly.update(document.getElementById(html_element_id), data, layout);
    }
};


function AssembleAjaxRequestUnified() {
    console.log("swicth to unified");
    cpm_is_active = $("#cpmTab").hasClass('is-active');
    orginal_is_active = $("#originalOrderTab").hasClass('is-active')
    var plot_type = 'original';
    var data_type = 'original';

    if (!cpm_is_active) {
        data_type = 'log10';
    } 

    if (!orginal_is_active) {
        plot_type = "hieracical";
    }
    // action here when clicking the search button
    var gene_name = $('#singleGene').val();

    $.ajax({
            type:'GET',
            url:'http://127.0.0.1:5000/data_unified',
            data: "gene=" + gene_name,
            success: function(result) {
                plotHeatmapUnified(result,"displayPlotUnified");
            },
            error: function (e) {
                alert('Request data Failed !')
            }
        });
        $("#originalTab").addClass('is-active');
}