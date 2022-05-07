// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var heatmapData = {};

function HeatmapByCelltype(result, html_element_id, dataScale, celltypeOrder) {
    if (!result) {
        alert("Error: no data to plot");
        return;
    }

    let x_axis;
    if (celltypeOrder == "original") {
        x_axis = result['celltypes'];
        y_axis = result['genes'];
    } else {
        x_axis = result['celltypes_hierarchical'];
        y_axis = result['genes_hierarchical'];
    }

    var ngenes =  y_axis.length;
    var graph_width = 1300;
    var graph_height = 270 + 26 * ngenes;

    let yticktext = [];
    for (let i = 0; i < y_axis.length; i++) {
        const gene = y_axis[i];
        const geneId = result['gene_ids'][gene];
        if (geneId === "") {
            yticktext.push(gene);
        } else {
            let geneUrl = gene;
            if (geneId.startsWith('MGI')) {
                geneUrl = 'http://www.informatics.jax.org/marker/'+geneId;
            } else {
                geneUrl = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='+geneId;
            }
            const tickText = '<a href="'+geneUrl+'">'+gene+'</a>'
            yticktext.push(tickText);
        }
    }

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
    var data = {
            type: 'heatmap',
            hoverongaps: false,
            colorscale: 'Reds',
        };

    // Make new plot if none is present
    if ($('#'+html_element_id).html() === "") {

        data['z'] = data_content;
        data['x'] = x_axis;
        data['y'] = y_axis;

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
                type: 'category',
            },
            yaxis: {
                //title: 'Genes',
                automargin: true,
                autorange: "reversed",
                type: 'category',
                tickvals: y_axis,
                ticktext: yticktext,
            },
        };
            
        Plotly.newPlot(
            document.getElementById(html_element_id),
            [data],
            layout,
        ); 

    // Update existing plot if present
    } else {
        data['z'] = [data_content];
        data['x'] = [x_axis];
        data['y'] = [y_axis];
        Plotly.update(
            document.getElementById(html_element_id),
            data,
            {
                height: graph_height,
                yaxis: {
                    autorange: "reversed",
                },
            },
            [0],
        ); 
    }
} 

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

function AssembleAjaxRequest() {
    // Get the list of genes to plot from the search box
    var gene_names = $('#searchGeneName').val();
  
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

            // Update search box: corrected gene names, excluding missing genes
            $('#searchGeneName').val(result['genes']);
  
            // Create heatmap
            updatePlot();
        },
        error: function (e) {
            console.log(e);
            alert('Error: Could not find some gene names.')
        },
    });
};

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


////////////////////
// EVENTS
////////////////////
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

// Both on click and load, plot the heatmap
$("#searchOnClick").click(AssembleAjaxRequest);
$(document).ready(AssembleAjaxRequest);
$("#geneSuggestions").click(onClickSuggestions);
