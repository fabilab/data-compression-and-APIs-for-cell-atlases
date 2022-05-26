// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var heatmapData = {};

function HeatmapByCelltype(result, html_element_id, dataScale, celltypeOrder) {
    if (!result) {
        alert("Error: no data to plot");
        return;
    }

    let x_axis, y_axis;
    if (celltypeOrder == "original") {
        x_axis = result['celltypes'];
        y_axis = result['genes'];
    } else {
        x_axis = [];
        for (let i = 0; i < result['celltypes_hierarchical'].length; i++) {
            const ct = result['celltypes'][result['celltypes_hierarchical'][i]];
            x_axis.push(ct);
        }
        y_axis = [];
        for (let i = 0; i < result['genes_hierarchical'].length; i++) {
            const gene = result['genes'][result['genes_hierarchical'][i]];
            y_axis.push(gene);
        }
    }

    var ngenes =  y_axis.length;
    var graph_width = 1300;
    var graph_height = 270 + 26 * ngenes;

    // Add hyperlinks to gene names
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

    // Fill heatmap data
    let data_content = [];
    if (celltypeOrder == "original") {
        for (let i = 0; i < y_axis.length; i++) {
            data_content.push([]);
            for (let j = 0; j < x_axis.length; j++) {
                let geneExp = result['data'][i][j];
                if (dataScale == "log10") {
                    geneExp = Math.log10(geneExp + 0.5);
                }
                data_content[i].push(geneExp);
            }
        }
    } else {
        for (let i = 0; i < y_axis.length; i++) {
            const ii = result['genes_hierarchical'][i];
            data_content.push([]);
            for (let j = 0; j < x_axis.length; j++) {
                const jj = result['celltypes_hierarchical'][j];
                let geneExp = result['data'][ii][jj];
                if (dataScale == "log10") {
                    geneExp = Math.log10(geneExp + 0.5);
                }
                data_content[i].push(geneExp);
            }
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
                    tickvals: y_axis,
                    ticktext: yticktext,
                },
            },
            [0],
        );
    }

    // Add tooltips to gene names
    $(".ytick > text").mousemove(function(evt) {
        let gene = $(this).text();
        let goTerms = result['GO_terms'][gene];
        if (goTerms === undefined) {
            return;
        }
        let text = "<b>GO terms:</b></br>";
        for (let i = 0; i < goTerms.length; i++) {
            text += '<div><a class="goHyperlink">'+goTerms[i]+"</a></div>";
        }
        showTooltip(evt, text);
    });
    $(".ytick > text").mouseout(function(evt) { hideTooltip(); });
}

function showTooltip(evt, text) {
    let tooltip = document.getElementById("tooltip");
    tooltip.innerHTML = text;
    tooltip.style.background = "white";
    tooltip.style.border = "1px solid black";
    tooltip.style.borderRadius = "5px";
    tooltip.style.padding = "5px";
    tooltip.style.display = "block";
    tooltip.style.left = evt.pageX - 2 + 'px';
    tooltip.style.top = evt.pageY - 2 + 'px';

    $(".goHyperlink").click(onClickGOTermSuggestions);
}

function hideTooltip() {
    var tooltip = document.getElementById("tooltip");
    setTimeout(function() {
      tooltip.style.display = "none";
    }, 4500);
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

function AssembleAjaxRequest( genestring = "" ) {
    // Get the list of genes to plot from the search box
    // If too long for the search box, it should be supplied as a specialGenestring
    let geneNames;
    if (genestring !== "") {
        geneNames = genestring;
    } else {
        geneNames = $('#searchGeneName').val();
    }

    let requestData = {
        gene_names: geneNames,
        species: species,
    }
    // HTML GET method length is capped at 4,000, but search box might be shorter
    let htmlVerb = (geneNames.length > 500) ? 'POST' : 'GET';

    // sent gene names to the API
    // FIXME: this fails at random times with large payloads?
    $.ajax({
        type: htmlVerb,
        url:'/data/by_celltype',
        data: $.param(requestData),
        success: function(result) {
            // Store global variable
            heatmapData = {
                'result': result,
                'div': 'h5_data_plot',
            };

            // Update search box: corrected gene names, excluding missing genes
            setSearchBox(result['genes']);

            // Create heatmap
            updatePlot();
        },
        error: function (e) {
            console.log(e);
            alert('Error: Could not find some gene names.')
        },
    });
};

// Check another species, same genes
function onClickSpeciesSuggestions() {
    var geneNames = $('#searchGeneName').val();
    const newSpecies = this.id.slice("suggest".length);
    let requestData = {
        newSpecies: newSpecies,
        gene_names: geneNames,
        species: species,
    }
    $.ajax({
        type:'GET',
        url:'/data/by_celltype',
        data: $.param(requestData),
        success: function(result) {
            // Store global variable
            heatmapData = {
                'result': result,
                'div': 'h5_data_plot',
            };
            $("#suggest"+newSpecies).text(species.slice(0, 1).toUpperCase()+species.slice(1)).prop('id', "suggest"+species);
            species = result['species'];

            // Update search box: corrected gene names, excluding missing genes
            setSearchBox(result['genes']);

            // Create heatmap
            updatePlot();
        },
        error: function (e) {
          alert('Error: Could not find orthologs for '+geneNames+'.')
        }
    });
}

// SuggestGenes: create a div with a "suggest" button
function onClickGeneSuggestions() {
    var geneNames = $('#searchGeneName').val();
    let requestData = {
        gene_names: geneNames,
        species: species,
    }
    $.ajax({
        type:'GET',
        url:'/data/gene_friends',
        data: $.param(requestData),
        success: function(result) {
            // Update search box: corrected gene names, excluding missing genes
            setSearchBox(result);

            // Request data
            AssembleAjaxRequest();
        },
        error: function (e) {
          alert('Error: Could not find gene friends for '+geneNames+'.')
        }
    });
}


// Request genes in a clicked GO term and update plot
function onClickGOTermSuggestions () {
    let goTerm = $(this).text();
    let requestData = {
        goTerm: goTerm,
        species: species,
    }

    hideTooltip();

    $.ajax({
        type:'GET',
        url:'/data/genes_in_go_term',
        data: $.param(requestData),
        success: function(result) {

            // Update search box: corrected gene names, excluding missing genes
            // NOTE: this is a "short text" HTML element, so it can crash if we have
            // many genes...
            if (result.length > 200) {
                setSearchBox(
                    result.slice(0, 10)+"...",
                    result.split(',').slice(0, 10).join(','),
                );
                // Request data
                AssembleAjaxRequest(result);
            } else {
                setSearchBox(result);
                AssembleAjaxRequest();
            }
        },
        error: function (e) {
          alert('Error: Could not find genes within GO term '+goTerm+'.')
        }
    });
}

function setSearchBox(text, gseaText = "") {
    $('#searchGeneName').val(text);
    // Sync with GSEA box
    if (gseaText == "") {
        gseaText = text;
    }
    $('#suggestGO > a').attr(
        'href', '/barplot_gsea?species='+species+'&genes='+gseaText);
    $('#suggestKEGG > a').attr(
        'href', '/barplot_gsea?species='+species+'&gene_set=KEGG&genes='+gseaText);
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
$("#searchOnClick").click(function() { AssembleAjaxRequest() });
$("body").keyup(function(event) {
    if (event.keyCode === 13) {
        $("#searchOnClick").click();
    }
});
$(document).ready(function() {
    $('#pathwaySuggestion > a').click(function() {
        $("body").addClass("loading");
    });
    AssembleAjaxRequest();
});
$("#geneSuggestions").click(onClickGeneSuggestions);
$(".speciesSuggestion").click(onClickSpeciesSuggestions);
