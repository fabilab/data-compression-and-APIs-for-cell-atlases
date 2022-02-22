// Normalise the data with log10 and generate a new plot (when user click the button)
$("#log10OnClick" ).click(function() {
    // if User has input their gene of interest, generate the heatmap with that value
    // otherwise use the default data
    var gene_name = $('#searchGeneName').val();
    if (gene_name === '') {
        gene_name = "Car4,Vwf,Col1a1,Ptprc,Ms4a1";
    }
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/dataLog',
        data: "gene_names=" + gene_name,
        dataType:'json',
        success: HeatMap,
        error: function (e) {
        alert('Request data Failed')
        }
    });

    $("#logTab").addClass('is-active');
    $("#originalTab").removeClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    $("#dendrogramTab").removeClass('is-active');
    
});