// Normalise the data with log10 and generate a new plot (when user click the button)
$("#log10OnClick" ).click(function() {
    // if User has input their gene of interest, generate the heatmap with that value
    // otherwise use the default data
    $("#logTab").addClass('is-active');
    $("#cpmTab").removeClass('is-active');
    AssembleAjaxRequest()
    
});

$("#CPMOnClick" ).click(function() {
    $("#logTab").removeClass('is-active');
    $("#cpmTab").addClass('is-active');
    AssembleAjaxRequest()
});

// Second set of buttons
$("#hClusterOnClick" ).click(function() {
    // if User has input their gene of interest, generate the heatmap with that value
    // otherwise use the default data
    $("#hierachicalTab").addClass('is-active');
    $("#originalOrderTab").removeClass('is-active');
    AssembleAjaxRequest()
});


$("#originalOnClick" ).click(function() {
    $("#originalOrderTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    AssembleAjaxRequest()
});
