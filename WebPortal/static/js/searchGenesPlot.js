// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
$( "#searchOnClick" ).click(function() {
  if(! $('#scatter_plot').is('empty')) {
    $('#scatter_plot').empty();
  }
  // action here when clicking the search button
  var gene_name = $('#searchGeneName').val();
  const gene_array = gene_name.split(",")
  if (gene_array.length == 2) {
    $.ajax({
      type:'GET',
      url:'http://127.0.0.1:5000/2_genes',
      data: "gene_names=" + gene_name,
      success: ScatterPlot,
      error: function (e) {
        alert('Request data Failed')
      }
    });
  }
    // sent gene names to the API
  $.ajax({
    type:'GET',
    url:'http://127.0.0.1:5000/dataOrigin',
    data: "gene_names=" + gene_name,
    success: HeatMap,
    error: function (e) {
      alert('Request data Failed')
    }
    });
});