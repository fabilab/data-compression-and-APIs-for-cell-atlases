dataForPlotsUnified = {}

function generateHmUnifiedGene(result,html_element_id) {
    // flags passed by events.js
    let useLog = dataForPlotsUnified['useLog'];
    let celltypeOrder = dataForPlotsUnified['celltypeOrder'];
    
    if (result === "") {
        result = dataForPlotsUnified['result'];
    } else {
        dataForPlotsUnified['result'] = result;
    }

    if (html_element_id === "") {
        html_element_id = "hm_unified_gene";
    }

    let celltypes;
    if (!celltypeOrder) {
        celltypes = result['cell_type'];
    } else {
        celltypes = result['hierarchicalCelltypeOrder'];
    }

    const expression = result['exp_avg'];
    
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
    
    for (var i = 0; i < result['dataset_timepoint'].length; i++) {
        // dataset timepoint combination: e.g TMS_24m
        let dataset_timepoint = result['dataset_timepoint'][i];
        let time  = dataset_timepoint.split("_")[1];
        let dataset = dataset_timepoint.split("_")[0];
        if (y_axis.includes(time)) {
            y_axis.push(time+'\'\'');
        } else {
            y_axis.push(time);
        }
        y_axis_time.push(time);
        y_axis_dataset.push(dataset);
        
        let expression_in_all_celltypes = [];
        for (var k = 0; k < celltypes.length; k++) {
            if (expression[dataset_timepoint][celltypes[k]] === -1) {
                expression_in_all_celltypes.push(null);
            } else {
                expression_in_all_celltypes.push(expression[dataset_timepoint][celltypes[k]]);
            }
        }

        if (useLog) {
            for (var j = 0; j < expression_in_all_celltypes.length; j++) {
                if (expression_in_all_celltypes[j] !== null) {
                    expression_in_all_celltypes[j] = Math.log10(expression_in_all_celltypes[j] + 0.5);
                }
            }
        }
        
        data_content.push(expression_in_all_celltypes);
    }
    // Generate hover text for each spot
    // Celltype: {ct}, Expression: {exp}, Dataset: {ds}, Timepoint: {tp}, 
    let hover_text = [];
    for (var i = 0; i < result['dataset_timepoint'].length; i++) {
        let temp = [];
        for (var j = 0; j < result['cell_type'].length; j++) {
            let dt = result['dataset_timepoint'][i];
            let ct = result['cell_type'][j];
            let exp;
            if(useLog && expression[dt][ct] !== null) {
                exp = Math.log10(expression[dt][ct]+0.5);
            } else{
                exp = expression[dt][ct];
            }
            // "ACZ_P21"
            let ds = dt.split("_")[0];
            let tp = dt.split("_")[1];
            temp.push('Celltype: '+ct+'<br>Dataset: '+ds+'<br>Timepoint: '+tp+'<br>Expression: '+exp);
        }
        hover_text.push(temp);
    }

    let ncelltypes = x_axis.length;
    let nTimepoints = y_axis.length;
    let heatmap_width = 1300;
    let heatmap_height = 270 + 30 * nTimepoints;
    var data = {
        type: 'heatmap',
        hoverinfo: 'text',
        colorscale: 'Reds',
    };
    var layout = {
        title: 'Expression profile of <b>' + result['gene'] + '</b> gene in all cell types over development time',
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
        hoverongaps: false,
    };
    var tools = {
        modeBarButtonsToAdd: [{
            name: 'Download plot as an SVG',
            icon: Plotly.Icons.camera,
            click: function(gd) {
              Plotly.downloadImage(gd, {format: 'svg'})
            }
          }]
    };

    if ($('#'+html_element_id).text() === "") {
        data['z'] = data_content;
        data['x'] = x_axis;
        data['y'] = y_axis;
        data['text'] = hover_text;
        Plotly.newPlot(document.getElementById(html_element_id), [data], layout,tools);
    } else {
        data['z'] = [data_content];
        data['x'] = [x_axis];
        data['y'] = [y_axis];
        data['text'] = [hover_text];
        Plotly.update(document.getElementById(html_element_id), data, layout,tools);
    }
};


function AjaxUnifiedGene() {
    // action here when clicking the search button
    var gene_name = $('#singleGene').val();

    $.ajax({
            type:'GET',
            url:'http://127.0.0.1:5000/data_unified',
            data: "gene=" + gene_name,
            success: function(result) {
                $("#hm_unified_gene").empty();
                generateHmUnifiedGene(result,"hm_unified_gene");
                $("#dp_unified_gene").empty();
                generateDpUnifiedGene(result, "dp_unified_gene");
            },
            error: function (e) {
                Swal.fire('Invalid input','please make sure you type in the correct gene name','error');
            }
        });
        $("#originalTab").addClass('is-active');
}