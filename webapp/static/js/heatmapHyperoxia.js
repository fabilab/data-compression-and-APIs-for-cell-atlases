// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var heatmapData = [];

function HeatmapHyperoxia(result, html_element_id, title, order) {
        if (!result) {
            alert("Error: Nothing to plot or gene names invalid")
        } else {
            var x_axis;
            if (order == "original") {
                x_axis = result['celltypes'];
            } else {
                x_axis = result['celltypes_hierarchical'];
            }
            var ngenes =  result['genes'].length;
            var graph_width = 1300;
            var graph_height = 370 + 26 * ngenes;

            let data_content = [];
            let i = 0;
            for (let gene in result['data']) {
                data_content.push([]);
                for (let j = 0; j < x_axis.length; j++) {
                    const ct = x_axis[j];
                    let gene_exp = result['data'][gene][ct];
                    data_content[i].push(gene_exp);
                }
                i+= 1;
            }
            var data = [
                {
                    z: data_content,
                    x: x_axis,
                    y: result['genes'],
                    type: 'heatmap',
                    hoverongaps: false,
                    colorscale: 'Reds',
                }
                ];
            if (result['data_scale'] === "log2FC") {
                data[0]['colorscale'] = 'RdBu';
                data[0]['zmid'] = 0;
            }

            var layout = {
                autosize: true,
                width: graph_width,
                height: graph_height,
                title: title,
                xaxis: {
                    //title: 'Cell types',
                    automargin: true,
                    tickangle: 60,
                    scaleanchor: 'y',
                    scaleratio: 1,
                    type: 'category',
                },
                yaxis: {
                    //title: 'Genes',
                    automargin: true,
                    autorange: "reversed",
                    type: 'category',
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
        url:'/data/hyperoxia',
        data: "gene_names="+gene_names+"&plot_type="+plot_type+"&data_type="+data_type,
        dataType:'json',
        success: function(result) {
            // Clear mobile DOM elements
            $("#h5_data_plot").html("");

            heatmapData = [];

            // result is a JSON list
            for(let i = 0; i < result.length; i++) {
                let item = result[i];
                const dataset = item['dataset'];
                const timepoint = item['timepoint'];
                const title = dataset+", "+timepoint;

                // Add DOM element
                const newId = "h5_data_plot_"+dataset+"_"+timepoint;
                $("#h5_data_plot").append("<div id="+newId+"></div>");

                // Add the data to persistency
                heatmapData.push({'item': item, 'div': newId, 'title': title});

                // Plot inside DOM element
                HeatmapHyperoxia(
                      item, 
                      newId,
                      title,
                      plot_type,
                    );
            }
        },
        error: function (e) {
            console.log(e);
            alert('Request data Failed');
        }
    });

};


function updateOrder(order) {
    // NOTE: heatmapData is the global persistent object
    for(let i = 0; i < heatmapData.length; i++) {
        HeatmapHyperoxia(
              heatmapData[i]['item'], 
              heatmapData[i]['div'],
              heatmapData[i]['title'],
              order,
            );
    }
}

// Both on click and load, plot the heatmap
$("#searchOnClick").click(function() {
  // action here when clicking the search button
  AssembleAjaxRequest();
});

$(document).ready(function() {
  AssembleAjaxRequest();
});


// normalization
$("#log2FCOnClick" ).click(function() {
    $("#log2FCTab").addClass('is-active');
    $("#cpmTab").removeClass('is-active');
    $("#logTab").removeClass('is-active');
    AssembleAjaxRequest()
});

$("#log10OnClick" ).click(function() {
    $("#log2FCTab").removeClass('is-active');
    $("#cpmTab").removeClass('is-active');
    $("#logTab").addClass('is-active');
    AssembleAjaxRequest()
});

$("#CPMOnClick" ).click(function() {
    $("#log2FCTab").removeClass('is-active');
    $("#cpmTab").addClass('is-active');
    $("#logTab").removeClass('is-active');
    AssembleAjaxRequest()
});


// order of cell types
$("#hClusterOnClick" ).click(function() {
    $("#hierachicalTab").addClass('is-active');
    $("#originalOrderTab").removeClass('is-active');
    updateOrder("hierarchical");
});


$("#originalOnClick" ).click(function() {
    $("#originalOrderTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    updateOrder("original");
});

