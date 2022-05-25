function plotHeatmapUnified(result) {
    const expression = result['expression']
    // x-axis: celltypes
    let x_axis = result['cell_type']
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
        let dataset_timepoint = result['dataset_timepoint'][i];   // get the dataset timepoint as a string
        let time  = dataset_timepoint.split("_")[1];
        let dataset = dataset_timepoint.split("_")[0];
        if (y_axis.includes(time)) {
            y_axis.push(" ");
        } else {
            y_axis.push(time);
        }
        y_axis_time.push(time);
        y_axis_dataset.push(dataset);
        expression_in_all_celltypes = expression[dataset_timepoint];  // find it from the dictionary as a key
        data_content.push(Object.values(expression_in_all_celltypes));
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
    var data = [
        {
            z: data_content,
            x: x_axis,
            y: y_axis,
            text: hover_text,
            type: 'heatmap',
            hoverinfo: 'text'
        }
    ];
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
        
    Plotly.newPlot(document.getElementById('bigHeatMap'), data,layout); 
}


function AssembleAjaxRequestUnified() {
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
    var gene_name = $('#searchGeneName').val();

    $.ajax({
            type:'GET',
            url:'http://127.0.0.1:5000/data_unified',
            data: "gene=" + gene_name + "&plottype=" + plot_type + "&datatype=" + data_type,
            success: function(result) {
                plotHeatmapUnified(result);
            },
            error: function (e) {
                alert('Request data Failed !')
            }
        });
        $("#originalTab").addClass('is-active');
}
$("#searchOnClick" ).click(AssembleAjaxRequestUnified)
$(document).ready(AssembleAjaxRequestUnified)