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
    success: function(result){
      console.log(result)
      // x-axis: 5 genes of interest
      let x_axis = Object.keys(result[Object.keys(result)[0]]);
      // y-axis:41 cell types
      let y_axis = Object.keys(result);
      let data_content = [];
      for (var i = 0; i < Object.keys(result).length; i++) {
          cell_type = Object.keys(result)[i] // get the cell_type name as a string
          all_gene_expression = result[cell_type]         // find it from the dictionary as a key
          data_content.push(Object.values(all_gene_expression))
      }
      var data = [
          {
              z: data_content,
              x: x_axis,
              y: y_axis,
              type: 'heatmap',
              hoverongaps: false
          }
        ];
      var layout = {
          title: 'Heatmap of gene expression level in selected cell types',
          xaixs: {
              title: {
                  text: 'Genes of interest',
                  standoff: 20
              }},
          yaxis: {
              title: {
                  text: 'Cell types',
                  standoff: 40
              }},
      };
        
        Plotly.newPlot(document.getElementById('h5_data_plot'), data,layout);  
    },
    error: function (e) {
      alert('Request data Failed')
    }
    });
});