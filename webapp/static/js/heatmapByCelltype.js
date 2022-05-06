// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var heatmapData = {};

function HeatmapByCelltype(result, html_element_id, dataScale, celltypeOrder) {
        if (!result) {
            alert("Error: no data to plot");
        } else {

            let x_axis;
            if (celltypeOrder == "original") {
                x_axis = result['celltypes'];
                y_axis = result['genes'];
            } else {
                x_axis = result['celltypes_hierarchical'];
                y_axis = result['genes_hierarchical'];
            }

            var y_axis = result['genes'];
            var ngenes =  y_axis.length;
            var graph_width = 1300;
            var graph_height = 270 + 26 * ngenes;

            let data_content = [];
            for (let i = 0; i < y_axis.length; i++) {
                const gene = y_axis[i];
                data_content.push([]);
                for (let j = 0; j < x_axis.length; j++) {
                    const ct = x_axis[j];
                    let gene_exp = result['data'][gene][ct]; 
                    if (dataScale == "log10") {
                        gene_exp = Math.log10(gene_exp + 0.5);
                    }
                    data_content[i].push(gene_exp);
                }
            }
            var data = [
                {
                    type: 'heatmap',
                    hoverongaps: false,
                    colorscale: 'Reds',
                }
                ];

            // Make new plot if none is present
            if ($('#'+html_element_id).html() === "") {

                data[0]['z'] = data_content;
                data[0]['x'] = x_axis;
                data[0]['y'] = y_axis;

                var layout = {
                    autosize: true,
                    width: graph_width,
                    height: graph_height,
                    title: 'Heatmap of gene expression level in selected cell types',
                    xaxis: {
                        //title: 'Cell types',
                        automargin: true,
                        tickangle: 70,
                        scaleanchor: 'y',
                        scaleratio: 1,
                    },
                    yaxis: {
                        //title: 'Genes',
                        automargin: true,
                        autorange: "reversed",
                    },
                };
                    
                Plotly.newPlot(
                    document.getElementById(html_element_id),
                    data,
                    layout,
                ); 

            // Update existing plot if present
            } else {
                data[0]['z'] = [data_content];
                data[0]['x'] = [x_axis];
                data[0]['y'] = [y_axis];
                Plotly.update(
                    document.getElementById(html_element_id),
                    data[0],
                    {yaxis: {autorange: "reversed"}},
                    [0],
                ); 

                // FIXME: how to animate color changes in a heatmap?
                //data[0]['z'] = data_content;
                //data[0]['x'] = x_axis;
                //data[0]['y'] = y_axis;
                //Plotly.animate(
                //    document.getElementById(html_element_id),
                //    {
                //        data: data,
                //    },
                //    {
                //        transition: {
                //            duration: 1500,
                //            easing: 'cubic-in-out',
                //        },
                //        frame: {
                //            duration: 1500,
	        //        },
                //    },
                //); 
            }
        };
    } 


// SuggestGenes: create a div with a "suggest" button
function onClickSuggestions() {
    var gene_names = $('#searchGeneName').val();
    $.ajax({
        type:'GET',
        url:'/gene_friends',
        data: "gene_names=" + gene_names,
        success: function(result) {
            $('#searchGeneName').val(result);
            AssembleAjaxRequest();
        },
        error: function (e) {
          alert('Error: Could not find gene friends for '+gene_names+'.')
        }
    });
}

function SuggestGenes() {

    // Fill div
    const html_element_id = "gene_suggestions";
    $("#"+html_element_id).empty();
    $("#"+html_element_id).text('Suggest similar genes');

    // Add link
    $("#"+html_element_id).click(onClickSuggestions);
}

// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
function AssembleAjaxRequest() {
  if(! $('#scatter_plot').is('empty')) {
    $('#scatter_plot').empty();
  }

  // Get the list of genes to plot from the search box
  var gene_names = $('#searchGeneName').val();

  const gene_array = gene_names.split(",")
  if (gene_array.length == 2) {
    $.ajax({
      type:'GET',
      url:'/2_genes',
      data: "gene_names=" + gene_names,
      success: ScatterPlot,
      error: function (e) {
        alert('Request data Failed')
      }
    });
  }
    // sent gene names to the API
  $.ajax({
    type:'GET',
    url:'/data/by_celltype',
    data: "gene_names=" + gene_names,
    success: function(result) {
        // Store global variable
        heatmapData = {
            'result': result,
            'div': 'h5_data_plot',
        };

        // Create heatmap
        updatePlot();

        // Create gene suggestions DOM element
        SuggestGenes();

        //FIXME: ScatterPlot should happen here
    },
    error: function (e) {
      alert('Error:Input gene name is invalid, please make sure you type in the corrent gene names.')
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


// NOTE: this is why react was invented...
function updatePlot() {
    let dataScale = "original";
    if (!$("#cpmTab").hasClass('is-active')) {
        dataScale = "log10";
    }
    let celltypeOrder = "original";
    if (!$("#originalOrderTab").hasClass('is-active')) {
        celltypeOrder = "hierarchical";
    }

    // NOTE: heatmapData is the global persistent object
    HeatmapByCelltype(
        heatmapData['result'], 
        heatmapData['div'],
        dataScale,
        celltypeOrder,
    );
}


// Normalise the data with log10 and generate a new plot (when user click the button)
$("#log10OnClick" ).click(function() {
    // if User has input their gene of interest, generate the heatmap with that value
    // otherwise use the default data
    $("#logTab").addClass('is-active');
    $("#cpmTab").removeClass('is-active');
    updatePlot();
    
});

$("#CPMOnClick" ).click(function() {
    $("#logTab").removeClass('is-active');
    $("#cpmTab").addClass('is-active');
    updatePlot();
});

// Second set of buttons
$("#hClusterOnClick" ).click(function() {
    // if User has input their gene of interest, generate the heatmap with that value
    // otherwise use the default data
    $("#hierachicalTab").addClass('is-active');
    $("#originalOrderTab").removeClass('is-active');
    updatePlot();
});


$("#originalOnClick" ).click(function() {
    $("#originalOrderTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    updatePlot();
});

