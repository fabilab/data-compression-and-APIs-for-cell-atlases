// Normalise the data with log10 and generate a new plot (when user click the button)
$("#log10OnClick" ).click(function() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/dataLog',
        data: "gene_names=Car4,Vwf,Col1a1,Ptprc,Ms4a1",
        dataType:'json',
        success: HeatMap,
        error: function (e) {
        alert('Request data Failed')
        }
    });

    $("#logTab").addClass('is-active');
    $("#originalTab").removeClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    
});