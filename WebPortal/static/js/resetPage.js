$( "#resetOnClick" ).click(function() {
    if(! $('#scatter_plot').is('empty')) {
        $('#scatter_plot').empty();
      }
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data',
        data: "gene_names=Car4,Vwf,Col1a1,Ptprc,Ms4a1&plot_type=original&data_type=original",
        dataType:'json',
        success: HeatMap,
        error: function (e) {
        alert('Request data Failed')
        }
    });
    $("#searchGeneName").val("");
    $("#originalTab").addClass('is-active');
    $("#logTab").removeClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
})