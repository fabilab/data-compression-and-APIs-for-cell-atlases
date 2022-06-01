var plotData = {};

function plotHeatmapUnified(result, scaleData, celltypeOrder) {
    let geneName = result['gene'];
    let geneId = result['gene_id'];
    let x_axis;
    if (celltypeOrder === "original") {
        x_axis = result['celltypes'];
    } else {
        x_axis = result['celltypes_hierarchical'];
    }
    let y_axis = result['row_labels'];
    let nx = x_axis.length;
    let ny = y_axis.length;

    let title;
    if (geneId === "") {
        title = geneName + ' expression over time';
    } else {
        if (geneId.startsWith('MGI')) {
            geneUrl = 'http://www.informatics.jax.org/marker/'+geneId;
        } else {
            geneUrl = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='+geneId;
        }
        title = '<a href="'+geneUrl+'">'+geneName+'</a> expression over time';
    }

    let x = [],
        y = [],
        tooltips = [],
        markersize = [],
        markeropacity = [],
        markercolor = [];
    let ms, opacity;
    for (let i = 0; i < y_axis.length; i++) {
        const label = y_axis[i];
        for (let j = 0; j < x_axis.length; j++) {
            const celltype = x_axis[j];
            let nc = result['ncells'][label][celltype]
            let ge = result['gene_expression'][label][celltype]
            if (scaleData == "log10") {
                ge = Math.log10(ge + 0.5);
            }
            const labelArray = label.split("_");
            const tooltip = "Expression: "+ge+", Dataset: "+labelArray[0]+", Time point: "+labelArray[1];
            if (nc < 5) {
                ms = 8;
                opacity = 0.5;
            } else if (nc < 40) {
                ms = 13;
                opacity = 0.7;
            } else {
                ms = 20;
                opacity = 0.9;
            }
            x.push(celltype)
            y.push(label)
            markercolor.push(ge);
            markeropacity.push(opacity);
            markersize.push(ms);
            tooltips.push(tooltip);
        }
    }

    let data = {
        mode: 'markers',
        marker: {
            symbol: 'square',
            colorscale: 'Reds',
            colorbar: {},
        },
        'hoverinfo': 'text',
    };

    // Make new plot if none is present
    if ($('#heatmapUnified').html() === "") {
        data['x'] = x;
        data['y'] = y;
        data['text'] = tooltips;
        data['marker']['color'] = markercolor;
        data['marker']['size'] = markersize;
        data['marker']['opacity'] = opacity;

        let layout = {
            automargin: true,
            autosize: true,
            width: 1300,
            height: 1000,
            title: {
                text: title,
                x: 0.5,
                xanchor: 'center',
            },
            xaxis: {
                tickangle: 70,
                automargin: true,
                linewidth: 0,
                type: 'category',
            },
            yaxis: {
                //automargin: true,
                autorange: 'reversed',
                type: 'category',
                tickvals: result['yticks'],
                ticktext: result['yticktext'],
            },
        };

        Plotly.newPlot(
            document.getElementById('heatmapUnified'),
            [data],
            layout,
        ); 

    // Update existing plot if present
    } else {
        data['x'] = [x];
        data['y'] = [y];
        data['text'] = [tooltips];
        data['marker']['color'] = markercolor;
        data['marker']['size'] = markersize;
        data['marker']['opacity'] = opacity;
        Plotly.update(
            document.getElementById('heatmapUnified'),
            data,
            {
                title: {
                    text: title,
                },
                yaxis: {
                autorange: "reversed",
                type: 'category',
                tickvals: result['yticks'],
                ticktext: result['yticktext'],
            }},
            [0],
        );
    }
}


function AssembleAjaxRequest() {
    var gene = $('#searchGeneName').val();
    let requestData = {
        gene: gene,
        species: species,
    }
    $.ajax({
        type:'GET',
        url:'/data/development',
        data: $.param(requestData),
        dataType:'json',
        success: function(result) {
            plotData = result;
            
            updateSimilarGenes();

            updatePlot();
        },
        error: function (e) {
          alert('Request data Failed')
        }
    });
}


function updateSimilarGenes() {
    let similarGenes = plotData['similarGenes'];
    $('#geneSuggestions').text('Similar genes:');
    for (let i = 0; i < similarGenes.length; i++) {
        const gene = similarGenes[i];
        $('#geneSuggestions').append(
            '<span class="geneSuggestion suggestButton" id="suggest'+gene+'">'+gene+'</span>'
        )
    }
    // Rebind the callback since the old elements are gone
    $(".geneSuggestion").click(onClickGeneSuggestions);
}

// Check another species, same gene
function onClickSpeciesSuggestions() {
    var geneName = $('#searchGeneName').val();
    const newSpecies = this.id.slice("suggest".length);
    let requestData = {
        newSpecies: newSpecies,
        gene: geneName,
        species: species,
    }
    $.ajax({
        type:'GET',
        url:'/data/development',
        data: $.param(requestData),
        success: function(result) {
            // Store global variable
            plotData = result;

            $("#suggest"+newSpecies).text(species.slice(0, 1).toUpperCase()+species.slice(1)).prop('id', "suggest"+species);
            species = newSpecies;

            updateSimilarGenes();

            // Update search box: corrected gene names, excluding missing genes
            setSearchBox(result['gene']);

            updatePlot();
        },
        error: function (e) {
          alert('Error: Could not find orthologs for '+geneName+'.')
        }
    });
}


// Check out a similar gene
function onClickGeneSuggestions() {
    var gene = $('#searchGeneName').val();
    var newGene = $(this).text();

    // swap button and search box
    // NOTE: this is not recursive, but easier to go back
    // API-wise, recursive approach would be more consistent
    $(this).text(gene);
    $('#searchGeneName').val(newGene);

    // Get new data and plot
    AssembleAjaxRequest();
}


function setSearchBox(text, gseaText = "") {
    $('#searchGeneName').val(text);
}



function updatePlot() {
    let scaleData, celltypeOrder;
    
    if ($("#cpmTab").hasClass('is-active')) {
      scaleData = "original";
    } else {
      scaleData = "log10";
    } 
    
    if ($("#originalOrderTab").hasClass('is-active')) {
        celltypeOrder = "original";
    } else {
        celltypeOrder = "hierarchical";
    }

    plotHeatmapUnified(plotData, scaleData, celltypeOrder);
}

$("#searchOnClick").click(AssembleAjaxRequest);
$(document).ready(AssembleAjaxRequest);
$(".geneSuggestion").click(onClickGeneSuggestions);
$(".speciesSuggestion").click(onClickSpeciesSuggestions);

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

$("body").keyup(function(event) {
    if (event.keyCode === 13) {
        $("#searchOnClick").click();
    }
});
