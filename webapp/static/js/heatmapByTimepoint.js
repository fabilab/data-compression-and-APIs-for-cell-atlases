function HeatMapTimepoint(result, html_element_id,dataset_name) {
    if (!result) {
        alert("Error:Input gene name is invalid, please make sure you type in the corrent gene names")
    } else {
        // x-axis: 5 genes of interest
        let x_axis = Object.keys(result[Object.keys(result)[0]]);
        // y-axis:41 cell types
        let y_axis = Object.keys(result);
        var ntimepoints = y_axis.length;
        var graph_width = 1300;
        var graph_height = 470 + 26 * ntimepoints;

        let data_content = [];
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
            title: 'Dataset: '+ dataset_name,
            width: graph_width,
            height: graph_height,
            xaxis: {
                //title: '<b>Celltypes<b>',
                automargin: true,
                tickangle: 60,
            },
            yaxis: {
                title: '<b>Timepoints<b>',
                automargin: true,
                scaleanchor: 'x',
                scaleratio: 1,
            },
        };
            
        Plotly.newPlot(document.getElementById(html_element_id), data,layout); 
    };
} 


function AssembleAjaxRequestTimepoint() {

  // When doing the search gene name action, we want it to be change immediatly without switching back to the original heatmap,
  //  for example, if we are looking at a log10 plot,and we do the search action, the tab stays at the log10 
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

  const gene_array = gene_name.split(",")
    // sent gene names to the API
  $.ajax({
    type:'GET',
    url:'/data_timepoint',
    data: "gene=" + gene_name + "&plottype=" + plot_type + "&datatype=" + data_type,
    success: function(result) {
        let div_id;
        for (const dataset in result) {
            div_id = 'dataset_' + dataset;
            HeatMapTimepoint(result[dataset], div_id, dataset);
        }
    },
    error: function (e) {
        alert('Request data Failed(in assemble timepoint)')
    }
    });
}

$("#searchOnClick" ).click(AssembleAjaxRequestTimepoint)

$(document).ready(function() {
    AssembleAjaxRequestTimepoint();
})
