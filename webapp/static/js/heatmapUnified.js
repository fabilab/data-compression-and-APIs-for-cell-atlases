var plotData = {};

function plotHeatmapUnified(result, scaleData, celltypeOrder) {
    let gene_name = result['gene'];
    let x_axis;
    if (celltypeOrder === "original") {
        x_axis = result['celltypes'];
    } else {
        x_axis = result['celltypes_hierarchical'];
    }
    let y_axis = result['row_labels'];
    let nx = x_axis.length;
    let ny = y_axis.length;

    let x = [],
        y = [],
        tooltip = [],
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
            tooltip.push("Expression: "+ge, "Label: "+label);
        }
    }

    let data = {
        x: x,
        y: y,
        text: tooltip,
        mode: 'markers',
        marker: {
            size: markersize,
            opacity: opacity,
            symbol: 'square',
            colorscale: 'Reds',
            color: markercolor,
        },
    };

    let layout = {
        automargin: true,
        autosize: true,
        width: 1300,
        height: 1000,
        title: {
            text: gene_name + ' expression over time',
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
        document.getElementById('heatmap_unified'),
        [data],
        layout,
    ); 
}


function AssembleAjaxRequest() {

    var gene_name = $('#searchGeneName').val();

    $.ajax({
        type:'GET',
        url:'/data_heatmap_unified',
        data: "gene=" + gene_name,
        dataType:'json',
        success: function(result) {
            plotData = result;
            updatePlot();
        },
        error: function (e) {
          alert('Request data Failed')
        }
    });
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

