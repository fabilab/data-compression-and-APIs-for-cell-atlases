function ShowGeneral() {
    $("#tabGeneral").addClass('is-active');
    $("#multiple_genes").removeClass('is-hidden');
    $("#scatterPlot").removeClass('is-hidden');
    $("#num_genes_div").removeClass('is-hidden');
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
    $("#dotPlot").addClass("is-hidden");
    $("#displayPlot").addClass("is-hidden");
    $("#num_genes_div").addClass('is-hidden');
}

function ShowDataset() {
    $("#tabDataset").addClass('is-active');
    $("#single_gene").removeClass('is-hidden');
    $("#timepoint_info").removeClass('is-hidden');
    $("#searchOnClick_single").click(AjaxExploreDatasets);
    $("#num_genes_div").addClass('is-hidden');
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        $("#dotPlotDataset").removeClass("is-hidden");
    }
    else {
        $("#displayPlotDataset").removeClass("is-hidden");
    }
}

function HideDataset() {
    $("#tabDataset").removeClass('is-active');
    $("#single_gene").addClass('is-hidden');
    $("#displayPlotDataset").addClass('is-hidden');
    $("#dotPlotDataset").addClass("is-hidden");
    $("#timepoint_info").addClass('is-hidden');
}

function ShowUnified() {
    $("#tabUnified").addClass('is-active');
    $("#single_gene").removeClass('is-hidden');
    $("#timepoint_info").removeClass('is-hidden');
    $("#searchOnClick_single").click(AjaxExploreUnified);
    $("#num_genes_div").addClass('is-hidden');
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        $("#dotPlotUnified").removeClass("is-hidden");
    }
    else {
        $("#displayPlotUnified").removeClass("is-hidden");
    }
}

function HideUnified() {
    $("#tabUnified").removeClass('is-active');
    $("#single_gene").addClass('is-hidden');
    $("#displayPlotUnified").addClass('is-hidden');
    $("#dotPlotUnified").addClass("is-hidden");
    $("#timepoint_info").addClass('is-hidden');

}

function ShowMarker() {
    $("#tabMarker").addClass('is-active');
    $("#marker_genes").removeClass('is-hidden');
    $("#sum_markers").removeClass("is-hidden");
    $("#num_genes_div").addClass('is-hidden');
    $("#dataOrder").addClass('is-hidden');
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        $("#dotPlotMarker").removeClass("is-hidden");
    }
    else {
        $("#displayPlotMarkers").removeClass("is-hidden");
    }
}

function HideMarker() {
    $("#tabMarker").removeClass('is-active');
    $("#marker_genes").addClass('is-hidden');
    $("#displayPlotMarkers").addClass('is-hidden');
    $("#sum_markers").addClass("is-hidden");
    $("#dataOrder").removeClass('is-hidden');
    $("#dotPlotMarker").addClass('is-hidden');
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

$("#tabUnified").click(function() {
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
        if($("#tabGeneral").hasClass('is-active')) {
            $("#dotPlot").removeClass("is-hidden");
            $("#displayPlot").addClass("is-hidden");
        } else if($("#tabDataset").hasClass('is-active')) {
            $("#displayPlotDataset").addClass("is-hidden");
            $("#dotPlotDataset").removeClass("is-hidden");
        } else if($("#tabUnified").hasClass('is-active')) {
            $("#displayPlotUnified").addClass("is-hidden");
            $("#dotPlotUnified").removeClass("is-hidden");
        } else if($("#tabMarker").hasClass('is-active')) {
            $("#displayPlotMarkers").addClass("is-hidden");
            $("#dotPlotMarker").removeClass("is-hidden");
        }
    }
    else {
        if($("#tabGeneral").hasClass('is-active')) {
            $("#dotPlot").addClass("is-hidden");
            $("#displayPlot").removeClass("is-hidden");
        } else if($("#tabDataset").hasClass('is-active')) {
            $("#displayPlotDataset").removeClass("is-hidden");
            $("#dotPlotDataset").addClass("is-hidden");
        } else if($("#tabUnified").hasClass('is-active')) {
            $("#displayPlotUnified").removeClass("is-hidden");
            $("#dotPlotUnified").addClass("is-hidden");
        } else if($("#tabMarker").hasClass('is-active')) {
            $("#displayPlotMarkers").removeClass("is-hidden");
            $("#dotPlotMarker").addClass("is-hidden");
        }
    }
});

$("#selectDataType").change(function() {
    let type = $("#selectDataType option:selected").val();
    dataAverageExp['useLog'] = (type === 'log');
    dataProportionExp['useLog'] = (type === 'log');
    dataForPlotsDataset['useLog'] = (type === 'log');
    dataForPlotsUnified['useLog'] = (type === 'log');
    dataMarker['useLog'] = (type === 'log');

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
        DotplotProportionExpUnifed("", "");
        plotHeatmapUnified("","");
    }

    if ($('#displayPlotMarkers').text() !== "") {
        HeatmapMarkerGenes("","","","", "","");
        DotplotProportionExpMarker("","","","","","", "");
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
        } else if ($('#displayPlotUnified').text() !== "") {
            plotHeatmapUnified("","");
            DotplotProportionExpUnifed("","");
        }
    } else{
        if ($('#displayPlot').text() !== "") {
            DotplotProportionExp("", "");
            HeatmapAverageExp("", "");
            plotHeatmapUnified("","");
            DotplotProportionExpUnifed("","");
        } else if ($('#displayPlotUnified').text() !== "") {
            plotHeatmapUnified("","");
            DotplotProportionExpUnifed("","");
        }
    }

    if ($('#dataset_1').text() !== "") {
        plotAll("");
    } 
});

$("#selectTopMarkers").change(function() {
    let showNum = $("#selectTopMarkers option:selected").val();
    if (showNum === 'top30') {
        HeatmapMarkerGenes("","","","", "", 30);
        DotplotProportionExpMarker("","","","","","", 30);
    } else if (showNum === 'top20') {
        HeatmapMarkerGenes("","","","", "", 20);
        DotplotProportionExpMarker("","","","","","", 20);

    } else if (showNum === 'top10') {
        HeatmapMarkerGenes("","","","", "", 10);
        DotplotProportionExpMarker("","","","","","", 10);
    } else {
        HeatmapMarkerGenes("","","","", "","");
        DotplotProportionExpMarker("","","","","","", "");
    }
})

$("#applyOnClick").click(AjaxExploreMarkers);
