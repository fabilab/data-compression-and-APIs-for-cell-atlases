plotDataUnifiedByCell = {}

function HeatmapUnifiedByCell(result,html_element_id) {
    // flags passed by events.js
    let useLog = plotDataUnifiedByCell['useLog'];
    let geneOrder = plotDataUnifiedByCell['geneOrder'];
    
    if (result === "") {
        result = plotDataUnifiedByCell['result'];
    } else {
        plotDataUnifiedByCell['result'] = result;
    }

    if (html_element_id === "") {
        html_element_id = "displayPlotUnifiedByCell";
    }
    let celltype = result['cell_type'];
    
    let genes;
    if (!geneOrder) {
        genes = result['genes'];
    } else {
        genes = result['clustered_gene_order'];
    }
    const expression = result['exp_avg'];
    
    // x-axis: celltypes
    let x_axis = genes;
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
        let dataset = dataset_timepoint.split("_")[0];
        let time  = dataset_timepoint.split("_")[1];
        if (y_axis.includes(time)) {
            y_axis.push(time+'\'\'');
        } else {
            y_axis.push(time);
        }
        y_axis_time.push(time);
        y_axis_dataset.push(dataset);
        
        let avgExpAll = [];
        for (var k = 0; k < genes.length; k++) {
            // console.log(dataset_timepoint, genes[k]);
            if (expression[dataset_timepoint][genes[k]] === -1) {
                avgExpAll.push(null);
            } else {
                avgExpAll.push(expression[dataset_timepoint][genes[k]]);
            }
        }
        if (useLog) {
            for (var j = 0; j < avgExpAll.length; j++) {
                if (avgExpAll[j] !== null) {
                    avgExpAll[j] = Math.log10(avgExpAll[j] + 0.5);
                }
            }
        }
        
        data_content.push(avgExpAll);
    }
    // Generate hover text for each spot
    // Celltype: {ct}, Expression: {exp}, Dataset: {ds}, Timepoint: {tp}, 
    let hover_text = [];
    for (var i = 0; i < result['dataset_timepoint'].length; i++) {
        let temp = [];
        for (var j = 0; j < result['cell_type'].length; j++) {
            let dt = result['dataset_timepoint'][i];
            let gn = result['genes'][j];
            let exp;
            if(useLog && expression[dt][gn] !== null) {
                exp = Math.log10(expression[dt][gn]+0.5);
            } else{
                exp = expression[dt][gn];
            }
            // "ACZ_P21"
            let ds = dt.split("_")[0];
            let tp = dt.split("_")[1];
            temp.push('Gene: '+gn+'<br>Dataset: '+ds+'<br>Timepoint: '+tp+'<br>Expression: '+exp);
        }
        hover_text.push(temp);
    }

    let nTimepoints = y_axis.length;
    let nGenes = x_axis.length;
    let heatmap_width = 100+55*nGenes;
    let heatmap_height = 200 + 30 * nTimepoints;
    console.log(nGenes);
    var data = {
        type: 'heatmap',
        hoverinfo: 'text',
        colorscale: 'Reds',
    };
    var layout = {
        title: result['cell_type'],
        xaxis: {
            title: '<b>Genes<b>',
            automargin: true,
            tickangle: 45,
            scaleanchor: 'y',
            scaleratio: 1,
            type: 'category'
        },
        yaxis: {
            title: '<b>Timepoints<b>',
            automargin: true,
            scaleanchor: 'x',
            scaleratio: 1,
            type: 'category'
        },
        width: heatmap_width,
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
    console.log(data_content);
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

function pagesetupUnified() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/all_cell_types',
        success: function(result) {
            var celltype_categories = Object.keys(result);
            for(var i=0;i<celltype_categories.length;i++) {
                var category = celltype_categories[i];
                $("#celltypeSelection").append(`<option class='has-text-weight-bold mt-2 has-text-success' disabled>${category}</option><br>`);
                var celltypes = result[category];
                for (var j=0;j<celltypes.length;j++) {
                    // $("#celltypeSelection").append(`<label class='radio'><input type='radio' name='celltype_selection' value='${celltypes[j]}'/> ${celltypes[j]}</label><br>`);
                    $("#celltypeSelection").append(`<option value='${celltypes[j]}'>${celltypes[j]}</option>`);
                }
            }
        },
        error: function (e) {
            alert('Request data fail (no cell types available)')
        }
    });
}

$(document).ready(pagesetupUnified)

function AjaxExploreUnifiedByCell() {
    // action here when clicking the search button
    let celltype = $("#celltypeSelection option:selected").val();
    let genes = $('#listGenes').val();

    $.ajax({
            type:'GET',
            url:'http://127.0.0.1:5000/data_unified_by_cell',
            data: "celltype=" + celltype + "&genes=" + genes,
            success: function(result) {
                $("#displayPlotUnifiedByCell").empty();
                HeatmapUnifiedByCell(result,"displayPlotUnifiedByCell");
                // $("#dotPlotUnifiedByCell").empty();
                // DotplotProportionExpUnifed(result, "dotPlotUnified");
            },
            error: function (e) {
                Swal.fire('Invalid input','please make sure you type in the correct gene name','error');
            }
        });
        $("#originalTab").addClass('is-active');
}
