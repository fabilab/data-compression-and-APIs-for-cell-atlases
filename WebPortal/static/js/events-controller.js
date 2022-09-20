function ShowGeneral() {
    $("#tabGeneral").addClass('is-active');
    $("#multiple_genes").removeClass('is-hidden');
    $("#scatter_plot").removeClass('is-hidden');
    $("#num_genes_div").removeClass('is-hidden');
    $("#timepoint_info").addClass('is-hidden');
    $("#searchOnClick_list").click(AjaxCompressed)
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        $("#dp_compressed").removeClass("is-hidden");
    }
    else {
        $("#hm_compressed").removeClass("is-hidden");
    }
}

function HideGeneral() {
    $("#tabGeneral").removeClass('is-active');
    $("#multiple_genes").addClass('is-hidden');
    $("#scatter_plot").addClass('is-hidden');
    $("#dp_compressed").addClass("is-hidden");
    $("#hm_compressed").addClass("is-hidden");
    $("#num_genes_div").addClass('is-hidden');
    $("#timepoint_info").removeClass('is-hidden');
}

function ShowDataset() {
    $("#tabDataset").addClass('is-active');
    $("#single_gene").removeClass('is-hidden');
    $("#searchOnClick_single").click(AjaxDatasets);
    $("#num_genes_div").addClass('is-hidden');
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        $("#dp_dataset").removeClass("is-hidden");
    }
    else {
        $("#hm_dataset").removeClass("is-hidden");
    }
}

function HideDataset() {
    $("#tabDataset").removeClass('is-active');
    $("#single_gene").addClass('is-hidden');
    $("#hm_dataset").addClass('is-hidden');
    $("#dp_dataset").addClass("is-hidden");
}

function ShowUnified() {
    $("#tabUnified").addClass('is-active');
    $("#single_gene").removeClass('is-hidden');
    $("#blank_info").removeClass('is-hidden');
    $("#searchOnClick_single").click(AjaxUnifiedGene);
    $("#num_genes_div").addClass('is-hidden');
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        $("#dp_unified_gene").removeClass("is-hidden");
    }
    else {
        $("#hm_unified_gene").removeClass("is-hidden");
    }
}

function HideUnified() {
    $("#tabUnified").removeClass('is-active');
    $("#single_gene").addClass('is-hidden');
    $("#hm_unified_gene").addClass('is-hidden');
    $("#dp_unified_gene").addClass("is-hidden");
    $("#blank_info").addClass('is-hidden');
}


function ShowUnifiedByCell() {
    $("#tabUnifiedByCell").addClass('is-active');
    $("#multiple_genes").removeClass('is-hidden');
    $("#selectCelltypeUnified").removeClass('is-hidden');
    $("#searchOnClick_list").click(AjaxUnifiedCell);
    $("#blank_info_2").removeClass('is-hidden');
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        $("#dp_unified_cell").removeClass("is-hidden");
    }
    else {
        $("#hm_unified_cell").removeClass("is-hidden");
    }
}

function HideUnifiedByCell() {
    $("#tabUnifiedByCell").removeClass('is-active');
    $("#multiple_genes").addClass('is-hidden');
    $("#selectCelltypeUnified").addClass('is-hidden');
    $("#hm_unified_cell").addClass('is-hidden');
    $("#blank_info_2").addClass('is-hidden');
    $("#dp_unified_cell").addClass('is-hidden');

}

function ShowMarker() {
    $("#tabMarker").addClass('is-active');
    $("#select_celltype").removeClass('is-hidden');
    $("#sum_markers").removeClass("is-hidden");
    $("#markers_info").removeClass("is-hidden");
    $("#num_genes_div").addClass('is-hidden');
    $("#timepoint_info").addClass('is-hidden');
    $("#dataOrder").addClass('is-hidden');
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        $("#dp_markers").removeClass("is-hidden");
    }
    else {
        $("#hm_markers").removeClass("is-hidden");
    }
}

function HideMarker() {
    $("#tabMarker").removeClass('is-active');
    $("#select_celltype").addClass('is-hidden');
    $("#hm_markers").addClass('is-hidden');
    $("#sum_markers").addClass("is-hidden");
    $("#markers_info").addClass("is-hidden");
    $("#dataOrder").removeClass('is-hidden');
    $("#dp_markers").addClass('is-hidden');
    $("#timepoint_info").removeClass('is-hidden');
}

$("#tabGeneral").click(function() {
    HideDataset();
    HideUnified();
    HideMarker();
    HideUnifiedByCell();
    ShowGeneral();
});

$("#tabDataset").click(function() {
    HideGeneral();
    HideUnified();
    HideMarker();
    HideUnifiedByCell();
    ShowDataset();
    
});

$("#tabUnified").click(function() {
    HideGeneral();
    HideDataset();
    HideMarker();
    HideUnifiedByCell();
    ShowUnified();
});

$("#tabUnifiedByCell").click(function() {
    HideGeneral();
    HideDataset();
    HideMarker();
    HideUnified();
    ShowUnifiedByCell();
});

$("#tabMarker").click(function() {
    HideGeneral();
    HideDataset();
    HideUnified();
    HideUnifiedByCell();
    ShowMarker();
});

$("#selectPlotType").change(function() {
    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        if($("#tabGeneral").hasClass('is-active')) {
            $("#hm_compressed").addClass("is-hidden");
            $("#dp_compressed").removeClass("is-hidden");
        } else if($("#tabDataset").hasClass('is-active')) {
            $("#hm_dataset").addClass("is-hidden");
            $("#dp_dataset").removeClass("is-hidden");
        } else if($("#tabUnified").hasClass('is-active')) {
            $("#hm_unified_gene").addClass("is-hidden");
            $("#dp_unified_gene").removeClass("is-hidden");
        } else if($("#tabUnifiedByCell").hasClass('is-active')) {
            $("#hm_unified_cell").addClass("is-hidden");
            $("#dp_unified_cell").removeClass("is-hidden");
        } else if($("#tabMarker").hasClass('is-active')) {
            $("#hm_markers").addClass("is-hidden");
            $("#dp_markers").removeClass("is-hidden");
        }
    }
    else {
        if($("#tabGeneral").hasClass('is-active')) {
            $("#hm_compressed").removeClass("is-hidden");
            $("#dp_compressed").addClass("is-hidden");
        } else if($("#tabDataset").hasClass('is-active')) {
            $("#hm_dataset").removeClass("is-hidden");
            $("#dp_dataset").addClass("is-hidden");
        } else if($("#tabUnified").hasClass('is-active')) {
            $("#hm_unified_gene").removeClass("is-hidden");
            $("#dp_unified_gene").addClass("is-hidden");
        } else if($("#tabUnifiedByCell").hasClass('is-active')) {
            $("#hm_unified_cell").removeClass("is-hidden");
            $("#dp_unified_cell").addClass("is-hidden");
        } else if($("#tabMarker").hasClass('is-active')) {
            $("#hm_markers").removeClass("is-hidden");
            $("#dp_markers").addClass("is-hidden");
        }
    }
});

$("#selectDataType").change(function() {
    let type = $("#selectDataType option:selected").val();
    dataAverageExp['useLog'] = (type === 'log');
    dataProportionExp['useLog'] = (type === 'log');
    dataForPlotsDataset['useLog'] = (type === 'log');
    dataForPlotsUnified['useLog'] = (type === 'log');
    dataForPlotsUnifiedCell['useLog'] = (type === 'log');
    dataMarker['useLog'] = (type === 'log');
    plotDataUnifiedByCell['useLog'] = (type === 'log');

    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        if ($('#dp_compressed').text() !== "") {
            generateDpCompressed("", "");
            generateHmCompressed("", "");
        }
    } else{
        if ($('#hm_compressed').text() !== "") {
            generateDpCompressed("", "");
            generateHmCompressed("", "");
        }
    }

    if ($('#dataset_1').text() !== "") {
        generateHmDatasets("", "");
        generateDpDatasets("", "");
    } 
    if ($('#dp_dataset_1').text() !== "") {
        generateHmDatasets("", "");
        generateDpDatasets("", "");
    } 

    if ($('#hm_unified_gene').text() !== "") {
        generateDpUnifiedGene("", "");
        generateHmUnifiedGene("","");
    }

    if ($('#hm_unified_cell').text() !== "") {
        generateDpUnifiedCell("", "");
        generateHmUnifiedCell("","");
    }

    if ($('#hm_markers').text() !== "") {
        generateHmMarkers("","","","", "","");
        generateDpMarkers("","","","","","", "");
    }
});

$("#selectDataOrder").change(function() {
    let order = $("#selectDataOrder option:selected").val();
    dataAverageExp['celltypeOrder'] = (order === 'clustered');
    dataProportionExp['celltypeOrder'] = (order === 'clustered');
    dataForPlotsDataset['celltypeOrder'] = (order === 'clustered');
    dataForPlotsUnified['celltypeOrder'] = (order === 'clustered');
    plotDataUnifiedByCell['geneOrder'] = (order === 'clustered')
    dataForPlotsUnifiedCell['geneOrder'] = (order === 'clustered');

    let plot = $("#selectPlotType option:selected").val();
    if (plot === "dot") {
        if ($('#dp_compressed').text() !== "") {
            generateDpCompressed("", "");
            generateHmCompressed("", "");
        } else if ($('#dataset_1').text() !== "") {
            generateHmDatasets("");
            generateDpDatasets("");
        } else if ($('#dp_dataset_1').text() !== "") {
            generateHmDatasets("");
            generateDpDatasets("");
        } else if ($('#hm_unified_gene').text() !== "") {
            generateHmUnifiedGene("","");
            generateDpUnifiedGene("","");
        } else if ($('#hm_unified_cell').text() !== "") {
            generateDpUnifiedCell("", "");
            generateHmUnifiedCell("","");
        }
    } else{
        if ($('#hm_compressed').text() !== "") {
            generateDpCompressed("", "");
            generateHmCompressed("", "");
            generateHmUnifiedGene("","");
            generateDpUnifiedGene("","");
        } else if ($('#dataset_1').text() !== "") {
            generateHmDatasets("");
            generateDpDatasets("");
        } else if ($('#hm_unified_gene').text() !== "") {
            generateHmUnifiedGene("","");
            generateDpUnifiedGene("","");
        } else if ($('#hm_unified_cell').text() !== "") {
            generateDpUnifiedCell("", "");
            generateHmUnifiedCell("","");
        }
    }
});

$("#selectTopMarkers").change(function() {
    let showNum = $("#selectTopMarkers option:selected").val();
    if (showNum === 'top30') {
        generateHmMarkers("","","","", "", 30);
        generateDpMarkers("","","","","","", 30);
    } else if (showNum === 'top20') {
        generateHmMarkers("","","","", "", 20);
        generateDpMarkers("","","","","","", 20);

    } else if (showNum === 'top10') {
        generateHmMarkers("","","","", "", 10);
        generateDpMarkers("","","","","","", 10);
    } else {
        generateHmMarkers("","","","", "","");
        generateDpMarkers("","","","","","", "");
    }
})

$("#applyOnClick").click(AjaxMarkers);
