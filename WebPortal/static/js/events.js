// Normalise the data with log10 and generate a new plot (when user click the button)
$("#log10OnClick" ).click(function() {
    $("#logTab").addClass('is-active');
    $("#cpmTab").removeClass('is-active');
    dataForPlots['useLog'] = true;
    HeatmapCelltype("", "");
    
});

$("#CPMOnClick" ).click(function() {
    $("#logTab").removeClass('is-active');
    $("#cpmTab").addClass('is-active');
    dataForPlots['useLog'] = false;
    HeatmapCelltype("", "");
});

// Second set of buttons
$("#hClusterOnClick" ).click(function() {
    $("#hierachicalTab").addClass('is-active');
    $("#originalOrderTab").removeClass('is-active');
    dataForPlots['celltypeOrder'] = true;
    HeatmapCelltype("", "");
});


$("#originalOnClick" ).click(function() {
    $("#originalOrderTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    dataForPlots['celltypeOrder'] = false;
    HeatmapCelltype("", "");
});
