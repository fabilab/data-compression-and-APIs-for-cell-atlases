$("#originalOnClick" ).click(function() {
    var gene_name = $('#searchGeneName').val();
    if (gene_name === '') {
        gene_name = "Car4,Vwf,Col1a1,Ptprc,Ms4a1";
    }
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data',
        data: "gene_names=" + gene_name + "&plot_type=original&data_type=original",
        dataType:'json',
        success: HeatMap,
        error: function (e) {
        alert('Request data Failed')
        }
    });
    $("#logTab").removeClass('is-active');
    $("#originalTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
})