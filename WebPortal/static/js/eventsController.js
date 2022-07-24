function ShowGeneral() {
    $("#tabGeneral").addClass('is-active');
    $("#multiple_genes").removeClass('is-hidden');
    $("#scatterPlot").removeClass('is-hidden');
    $("#plotType").removeClass('is-hidden');
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        $("#dotPlot").removeClass("is-hidden");
    }
    else {
        $("#displayPlot").removeClass("is-hidden");
    }
}

function HideGeneral() {
    $("#tabGeneral").removeClass('is-active');
    $("#multiple_genes").addClass('is-hidden');
    $("#scatterPlot").addClass('is-hidden');
    $("#plotType").addClass('is-hidden');
    $("#dotPlot").addClass("is-hidden");
    $("#displayPlot").addClass("is-hidden");
}

function ShowDataset() {
    $("#tabDataset").addClass('is-active');
    $("#single_gene").removeClass('is-hidden');
    $("#displayPlotDataset").removeClass('is-hidden');
    $("#timepoint_info").removeClass('is-hidden');
    $("#searchOnClick_single").click(AjaxExploreDatasets);
}

function HideDataset() {
    $("#tabDataset").removeClass('is-active');
    $("#single_gene").addClass('is-hidden');
    $("#displayPlotDataset").addClass('is-hidden');
    $("#timepoint_info").addClass('is-hidden');
}

function ShowUnified() {
    $("#tabTimepoints").addClass('is-active');
    $("#single_gene").removeClass('is-hidden');
    $("#displayPlotUnified").removeClass('is-hidden');
    $("#timepoint_info").removeClass('is-hidden');
    $("#searchOnClick_single").click(AjaxExploreUnified);
}

function HideUnified() {
    $("#tabTimepoints").removeClass('is-active');
    $("#single_gene").addClass('is-hidden');
    $("#displayPlotUnified").addClass('is-hidden');
    $("#timepoint_info").addClass('is-hidden');
}


function ShowMarker() {
    $("#tabMarker").addClass('is-active');
    $("#marker_genes").removeClass('is-hidden');
    $("#displayPlotMarkers").removeClass('is-hidden');
}

function HideMarker() {
    $("#tabMarker").removeClass('is-active');
    $("#marker_genes").addClass('is-hidden');
    $("#displayPlotMarkers").addClass('is-hidden');
}

$("#tabGeneral").click(function() {
    HideDataset();
    HideUnified();
    HideMarker();
    ShowGeneral();
});

$("#tabDataset").click(function() {
    HideGeneral();
    HideUnified();
    HideMarker();
    ShowDataset();
});

$("#tabTimepoints").click(function() {
    HideGeneral();
    HideDataset();
    HideMarker();
    ShowUnified();
});

$("#tabMarker").click(function() {
    HideGeneral();
    HideDataset();
    HideUnified();
    ShowMarker();
});

$("#selectPlotType").change(function() {
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        $("#dotPlot").removeClass("is-hidden");
        $("#displayPlot").addClass("is-hidden");

    }
    else {
        $("#dotPlot").addClass("is-hidden");
        $("#displayPlot").removeClass("is-hidden");
    }
});

$("#selectDataType").change(function() {
    let type = $("#selectDataType option:selected").val();
    dataAverageExp['useLog'] = (type === 'log');
    dataProportionExp['useLog'] = (type === 'log');
    dataForPlotsDataset['useLog'] = (type === 'log');
    dataForPlotsUnified['useLog'] = (type === 'log');

    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        if ($('#dotPlot').text() !== "") {
            DotplotProportionExp("", "");
            HeatmapAverageExp("", "");
        }
    } else{
        if ($('#displayPlot').text() !== "") {
            DotplotProportionExp("", "");
            HeatmapAverageExp("", "");
        }
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
    dataAverageExp['celltypeOrder'] = (order === 'clustered');
    dataProportionExp['celltypeOrder'] = (order === 'clustered');
    dataForPlotsDataset['celltypeOrder'] = (order === 'clustered');
    dataForPlotsUnified['celltypeOrder'] = (order === 'clustered');

    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        if ($('#dotPlot').text() !== "") {
            DotplotProportionExp("", "");
            HeatmapAverageExp("", "");
        }
    } else{
        if ($('#displayPlot').text() !== "") {
            DotplotProportionExp("", "");
            HeatmapAverageExp("", "");
        } 
    }

    if ($('#dataset_1').text() !== "") {
        plotAll("");
    } 

    if ($('#displayPlotUnified').text() !== "") {
        plotHeatmapUnified("","");
    }
});


$("#applyOnClick").click(AjaxExploreMarkers);
