function plotHeatmapUnified(result) {
    const expression = result['expression']
    // x-axis: celltypes
    let x_axis = Object.keys(expression[Object.keys(expression)[0]]);
    // y-axis: dataset timepoint
    let y_axis = result['dataset_timepoint'];
    let data_content = [];
    
    let y_axis_time = [];
    let y_axis_dataset = [];
    // plot only the timepoint on Y-axis
    for (var i = 0; i < result['dataset_timepoint'].length; i++) {
        let dataset_timepoint = result['dataset_timepoint'][i];   // get the dataset timepoint as a string
        let time  = dataset_timepoint.split("_")[1];
        let dataset = dataset_timepoint.split("_")[0];
        if (y_axis_time.includes(time)) {
            y_axis_time.push(" ");
        } else {
            y_axis_time.push(time);
        }
        y_axis_dataset.push(dataset);
        expression_in_all_celltypes = expression[dataset_timepoint];  // find it from the dictionary as a key
        data_content.push(Object.values(expression_in_all_celltypes));
    }
    console.log(y_axis_dataset);
    let ncelltypes = x_axis.length;
    let ngenes = y_axis.length;
    let heatmap_width = 1300;
    let heatmap_height = 270 + 41 * ngenes;
    var data = [
        {
            z: data_content,
            x: x_axis,
            y: y_axis_time,
            text: y_axis_dataset,
            type: 'heatmap',
            hovertemplate: '<i>Expression</i>: %{z}' + '<br>' +
                            '<i>Cell type</i>: %{x}' + '<br>' +
                            '<i>Dataset</i>: %{text}' + '<br>' +
                            '<i>Timepoint</i>: %{y}',
        }
        ];
    var layout = {
        title: 'expression of ' + result['gene'] + ' gene over time',
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
        // $("#originalTab").addClass('is-active');
}
// $("#searchOnClick" ).click(AssembleAjaxRequestTimepoint)
$(document).ready(AssembleAjaxRequestUnified)