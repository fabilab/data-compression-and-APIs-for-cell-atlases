// Normalise the data with log10 and generate a new plot (when user click the button)
$("#log10OnClick" ).click(function() {
    // if User has input their gene of interest, generate the heatmap with that value
    // otherwise use the default data
    $("#logTab").addClass('is-active');
    $("#cpmTab").removeClass('is-active');
    dataForPlots['useLog'] = true;
    AssembleAjaxRequestTimepoint()
    
});

$("#CPMOnClick" ).click(function() {
    $("#logTab").removeClass('is-active');
    $("#cpmTab").addClass('is-active');
    dataForPlots['useLog'] = false;
    AssembleAjaxRequestTimepoint()
});

// Second set of buttons
$("#hClusterOnClick" ).click(function() {
    // if User has input their gene of interest, generate the heatmap with that value
    // otherwise use the default data
    $("#hierachicalTab").addClass('is-active');
    $("#originalOrderTab").removeClass('is-active');
    dataForPlots['celltypeOrder'] = true;
    AssembleAjaxRequestTimepoint()
});


$("#originalOnClick" ).click(function() {
    $("#originalOrderTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    dataForPlots['celltypeOrder'] = false;
    AssembleAjaxRequestTimepoint()
});
