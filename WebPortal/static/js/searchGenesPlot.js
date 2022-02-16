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
      success: function(result){
        // compare the expression level of 2 specific genes across all cell types
        let expr1 = result["gene1_expr"];
        let expr2 = result["gene2_expr"];
        let cell_types =result["cell_types"];
        var trace1 = {
            x: expr1,
            y: expr2,
            mode: 'markers+text',
            type: 'scatter',
            text: cell_types,
            textposition: 'bottom center',
            textfont: {
              family:  'Raleway, sans-serif'
            },
            marker: { size: 12 }
        };
            
        var data = [ trace1 ];
        
        var layout = {
            xaxis: {
                title: {
                  text: result['gene1_name']+' expression',
                },
                range: [ 0, 5 ]
            },
            yaxis: {
                title: {
                  text: result['gene2_name']+' expression',
              },
              range: [0, 5]
            },
            title:'Correlation between ' + result['gene1_name'] + ' and ' + result['gene2_name']
        };
        
        Plotly.newPlot(document.getElementById('scatter_plot'), data, layout);
      }
    });
  }
    // sent gene names to the API
  $.ajax({
    type:'GET',
    url:'http://127.0.0.1:5000/data',
    data: "gene_names=" + gene_name,
    success: HeatMap,
    error: function (e) {
      alert('Request data Failed')
    }
    });
});