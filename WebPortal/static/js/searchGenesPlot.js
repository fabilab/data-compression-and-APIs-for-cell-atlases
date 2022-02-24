// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
function AssembleAjaxRequest() {
  if(! $('#scatter_plot').is('empty')) {
    $('#scatter_plot').empty();
  }

  // When doing the search gene name action, we want it to be change immediatly without switching back to the original heatmap,
  //  for example, if we are looking at a log10 plot,and we do the search action, the tab stays at the log10 
  cpm_is_active = $("#cpmTab").hasClass('is-active');
  orginal_is_active = $("#originalOrderTab").hasClass('is-active')
  if (cpm_is_active && orginal_is_active) {
    var path_name = "dataOrigin";
  } else if (cpm_is_active) {
    var path_name = "dataHierachicalOriginal"
  } else if (orginal_is_active) {
    var path_name = "dataLog"
  } else {
    var path_name = "dataHierachical"
  }
  
  // action here when clicking the search button
  var gene_name = $('#searchGeneName').val();
  if (gene_name === '') {
    gene_name = "Car4,Vwf,Col1a1,Ptprc,Ms4a1";
  }

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
    url:'http://127.0.0.1:5000/' + path_name,
    data: "gene_names=" + gene_name,
    success: HeatMap,
    error: function (e) {
      alert('Error:Input gene name is invalid, please make sure you type in the corrent gene names.')
    }
    });
  }
$("#searchOnClick" ).click(AssembleAjaxRequest)
