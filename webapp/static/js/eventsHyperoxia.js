// normalization
$("#log2FCOnClick" ).click(function() {
    $("#log2FCTab").addClass('is-active');
    $("#cpmTab").removeClass('is-active');
    $("#logTab").removeClass('is-active');
    AssembleAjaxRequest()
});

$("#log10OnClick" ).click(function() {
    $("#log2FCTab").removeClass('is-active');
    $("#cpmTab").removeClass('is-active');
    $("#logTab").addClass('is-active');
    AssembleAjaxRequest()
});

$("#CPMOnClick" ).click(function() {
    $("#log2FCTab").removeClass('is-active');
    $("#cpmTab").addClass('is-active');
    $("#logTab").removeClass('is-active');
    AssembleAjaxRequest()
});


// order of cell types
$("#hClusterOnClick" ).click(function() {
    $("#hierachicalTab").addClass('is-active');
    $("#originalOrderTab").removeClass('is-active');
    AssembleAjaxRequest()
});


$("#originalOnClick" ).click(function() {
    $("#originalOrderTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    AssembleAjaxRequest()
});
