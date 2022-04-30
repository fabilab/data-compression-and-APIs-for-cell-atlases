// Plot heatmap by celltype as a callback for the AJAX request
function HeatmapHyperoxia(result, html_element_id, title, dataScale) {
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
            if (dataScale === "log2FC") {
                data[0]['colorscale'] = 'RdBu';
            }

            var layout = {
                autosize: true,
                width: graph_width,
                height: graph_height,
                title: title,
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


// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
function AssembleAjaxRequest() {

    // Get the list of genes to plot from the search box
    var gene_names = $('#searchGeneName').val();

    // When doing the search gene name action, we want it to be change immediatly without switching back to the original heatmap,
    //  for example, if we are looking at a log10 plot,and we do the search action, the tab stays at the log10 
    let plot_type, data_type;
    
    if ($("#log2FCTab").hasClass('is-active')) {
        data_type = 'log2FC';
    } else if ($("#cpmTab").hasClass('is-active')) {
        data_type = "cpm";
    } else {
        data_type = "log10";
    }
    
    if ($("#originalOrderTab").hasClass('is-active')) {
      plot_type = "original";
    } else {
      plot_type = "hierachical";
    }

    // sent gene names to the API
    $.ajax({
        type:'GET',
        url:'/data_hyperoxia',
        data: "gene_names="+gene_names+"&plot_type="+plot_type+"&data_type="+data_type,
        dataType:'json',
        success: function(result) {
            // Clear mobile DOM elements
            $("#h5_data_plot").html("");

            // result is a JSON list
            // FIXME
            for(let i = 0; i < 1; i++) {
                let item = result[i];
                const dataset = item['dataset'];
                const timepoint = item['timepoint'];

                // Add DOM element
                const newId = "h5_data_plot_"+dataset+"_"+timepoint;
                $("#h5_data_plot").append("<div id="+newId+"></div>");

                // Plot inside DOM element
                HeatmapHyperoxia(
                      item['data'], 
                      newId,
                      dataset+", "+timepoint,
                      item['data_scale'],
                    );
            }
        },
        error: function (e) {
            console.log(e);
            alert('Request data Failed');
        }
    });

};

// Both on click and load, plot the heatmap
$("#searchOnClick").click(function() {
  // action here when clicking the search button
  AssembleAjaxRequest();
});

$(document).ready(function() {
  AssembleAjaxRequest();
});
