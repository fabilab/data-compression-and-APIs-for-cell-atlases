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

$("#tabGeneral").click(function() {
    $("#tabGeneral").addClass('is-active');
    $("#tabDataset").removeClass('is-active');
    $("#tabTimepoints").removeClass('is-active');
    $("#tabMarker").removeClass('is-active');
});

$("#tabDataset").click(function() {
    $("#tabGeneral").removeClass('is-active');
    $("#tabDataset").addClass('is-active');
    $("#tabTimepoints").removeClass('is-active');
    $("#tabMarker").removeClass('is-active');
});

$("#tabTimepoints").click(function() {
    $("#tabGeneral").removeClass('is-active');
    $("#tabDataset").removeClass('is-active');
    $("#tabTimepoints").addClass('is-active');
    $("#tabMarker").removeClass('is-active');
});

$("#tabMarker").click(function() {
    $("#tabGeneral").removeClass('is-active');
    $("#tabDataset").removeClass('is-active');
    $("#tabTimepoints").removeClass('is-active');
    $("#tabMarker").addClass('is-active');
});

$("#selectDataType").change(function() {
    let type = $("#selectDataType option:selected").val();
    dataForPlots['useLog'] = (type === 'log');
    HeatmapCelltype("", "");
});

$("#selectDataOrder").change(function() {
    let order = $("#selectDataOrder option:selected").val();
    dataForPlots['celltypeOrder'] = (order === 'clustered');
    HeatmapCelltype("", "");
});
