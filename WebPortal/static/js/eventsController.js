$("#tabGeneral").click(function() {
    $("#tabGeneral").addClass('is-active');
    $("#tabDataset").removeClass('is-active');
    $("#tabTimepoints").removeClass('is-active');
    $("#tabMarker").removeClass('is-active');
    // remove the single search button
    $("#single_gene").addClass('is-hidden');
    $("#multiple_genes").removeClass('is-hidden');
    $("#displayPlot_dataset").addClass('is-hidden');
    $("#displayPlot").removeClass('is-hidden');
    $("#scatter_plot").removeClass('is-hidden');
    $("#displayPlotUnified").addClass('is-hidden');
    $("#timepoint_info").addClass('is-hidden');
    $("#marker_genes").addClass('is-hidden');
    $("#displayPlotMarkers").addClass('is-hidden');
});

$("#tabDataset").click(function() {
    $("#tabGeneral").removeClass('is-active');
    $("#tabDataset").addClass('is-active');
    $("#tabTimepoints").removeClass('is-active');
    $("#tabMarker").removeClass('is-active');
    // remove the multi-search button
    $("#multiple_genes").addClass('is-hidden');
    $("#single_gene").removeClass('is-hidden');
    $("#displayPlot_dataset").removeClass('is-hidden');
    $("#displayPlot").addClass('is-hidden');
    $("#scatter_plot").addClass('is-hidden');
    $("#displayPlotUnified").addClass('is-hidden');
    $("#searchOnClick_single").click(AjaxExploreDatasets);
    $("#timepoint_info").removeClass('is-hidden');
    $("#marker_genes").addClass('is-hidden');
    $("#displayPlotMarkers").addClass('is-hidden');
});

$("#tabTimepoints").click(function() {
    $("#tabGeneral").removeClass('is-active');
    $("#tabDataset").removeClass('is-active');
    $("#tabTimepoints").addClass('is-active');
    $("#tabMarker").removeClass('is-active');
    $("#multiple_genes").addClass('is-hidden')
    $("#single_gene").removeClass('is-hidden')
    $("#displayPlotUnified").removeClass('is-hidden');
    $("#displayPlot_dataset").addClass('is-hidden');
    $("#displayPlot").addClass('is-hidden');
    $("#scatter_plot").addClass('is-hidden');
    $("#searchOnClick_single").click(AjaxExploreUnified);
    $("#timepoint_info").removeClass('is-hidden');
    $("#marker_genes").addClass('is-hidden');
    $("#displayPlotMarkers").addClass('is-hidden');
});

$("#tabMarker").click(function() {
    $("#tabGeneral").removeClass('is-active');
    $("#tabDataset").removeClass('is-active');
    $("#tabTimepoints").removeClass('is-active');
    $("#tabMarker").addClass('is-active');
    $("#timepoint_info").addClass('is-hidden');
    $("#marker_genes").removeClass('is-hidden');
    $("#multiple_genes").addClass('is-hidden')
    $("#single_gene").addClass('is-hidden');
    $("#applyOnClick").click(AjaxExploreMarkers);
    $("#displayPlotMarkers").removeClass('is-hidden');
    $("#displayPlot_dataset").addClass('is-hidden');
    $("#scatter_plot").addClass('is-hidden');
    $("#displayPlotUnified").addClass('is-hidden');
    $("#displayPlot").addClass('is-hidden');
});

$("#selectDataType").change(function() {
    let type = $("#selectDataType option:selected").val();
    dataAverageExp['useLog'] = (type === 'log');
    dataForPlotsDataset['useLog'] = (type === 'log');
    dataForPlotsUnified['useLog'] = (type === 'log');

    if ($('#displayPlot').text() !== "") {
        // HeatmapCelltype("", "");
        DotplotProportionExp("", "");
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
