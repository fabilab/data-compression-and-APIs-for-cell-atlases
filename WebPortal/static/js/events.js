$("#tabGeneral").click(function() {
    $("#tabGeneral").addClass('is-active');
    $("#tabDataset").removeClass('is-active');
    $("#tabTimepoints").removeClass('is-active');
    $("#tabMarker").removeClass('is-active');
    // remove the single search button
    $("#single_gene").hide();
    $("#multiple_genes").show();
    $("#displayPlot_dataset").hide();
    $("#displayPlot").show();
    $("#scatter_plot").show();
    $("#displayPlotUnified").hide();
});

$("#tabDataset").click(function() {
    $("#tabGeneral").removeClass('is-active');
    $("#tabDataset").addClass('is-active');
    $("#tabTimepoints").removeClass('is-active');
    $("#tabMarker").removeClass('is-active');
    // remove the multi-search button
    $("#multiple_genes").hide();
    $("#single_gene").show();
    $("#displayPlot_dataset").show();
    $("#displayPlot").hide();
    $("#scatter_plot").hide();
    $("#displayPlotUnified").hide();
    $("#searchOnClick_single").click(AssembleAjaxRequestTimepoint);
});

$("#tabTimepoints").click(function() {
    $("#tabGeneral").removeClass('is-active');
    $("#tabDataset").removeClass('is-active');
    $("#tabTimepoints").addClass('is-active');
    $("#tabMarker").removeClass('is-active');
    $("#multiple_genes").hide()
    $("#single_gene").show()
    $("#displayPlotUnified").show();
    $("#displayPlot_dataset").hide();
    $("#displayPlot").hide();
    $("#scatter_plot").hide();
    $("#searchOnClick_single" ).click(AssembleAjaxRequestUnified);
});

$("#tabMarker").click(function() {
    $("#tabGeneral").removeClass('is-active');
    $("#tabDataset").removeClass('is-active');
    $("#tabTimepoints").removeClass('is-active');
    $("#tabMarker").addClass('is-active');
});

$("#selectDataType").change(function() {
    let type = $("#selectDataType option:selected").val();
    dataForPlotsCellType['useLog'] = (type === 'log');
    dataForPlotsDataset['useLog'] = (type === 'log');
    dataForPlotsUnified['useLog'] = (type === 'log');

    if ($('#displayPlot').text() !== "") {
        HeatmapCelltype("", "");
    } 

    if ($('#dataset_1').text() !== "") {
        plotAll("");
    } 

    if ($('#displayPlotUnified').text() !== "") {
        plotHeatmapUnified("","");
    }

    
});

$("#selectDataOrder").change(function() {
    let order = $("#selectDataOrder option:selected").val();
    dataForPlotsCellType['celltypeOrder'] = (order === 'clustered');
    dataForPlotsDataset['celltypeOrder'] = (order === 'clustered');
    dataForPlotsUnified['celltypeOrder'] = (order === 'clustered');

    if ($('#displayPlot').text() !== "") {
        HeatmapCelltype("", "");
    } 

    if ($('#dataset_1').text() !== "") {
        plotAll("");
    } 

    if ($('#displayPlotUnified').text() !== "") {
        plotHeatmapUnified("","");
    }
});
