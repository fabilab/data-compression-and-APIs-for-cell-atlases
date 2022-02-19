$("#hClusterOnClick" ).click(function() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/dataHierachical',
        data: "gene_names=Car4,Vwf,Col1a1,Ptprc,Ms4a1",
        dataType:'json',
        success: HeatMap,
        error: function (e) {
        alert('Request data Failed')
        }
    });
    $("#logTab").removeClass('is-active');
    $("#originalTab").removeClass('is-active');
    $("#hierachicalTab").addClass('is-active');
});