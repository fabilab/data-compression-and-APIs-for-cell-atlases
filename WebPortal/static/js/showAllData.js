$(document).ready(function() {
  $.ajax({
    type:'GET',
    url:'http://127.0.0.1:5000/data',
    dataType:'json',
    success: function (result) {
      // remove the first column (depends on the given data set)
      delete result['Unnamed: 0'];
      let num_rows = Object.keys(result[Object.keys(result)[0]]).length;
      // display header to the table
      for(var i = 0; i < Object.keys(result).length; i++) {
        var header_field = $('<th>' + Object.keys(result)[i] + '</th>');
        // console.log(header_field);
        $('#table_header').append(header_field);
      }

      // display values to the table:
      // each row is a record for a patient
      let plot_data = [];
      for(var j = 0; j <num_rows; j++) {
        let patient_expr_values = [];
        var patient_row = $('<tr></tr>');
        for (var i = 0; i < Object.keys(result).length; i++) {
          let field = result[Object.keys(result)[i]][j];
          patient_row.append($('<td>' + field + '</td>'));
          // only store/push the "numerical" data 
          if (i >= 3) {
            patient_expr_values.push(field);
          }
        }
        $('#data_table').append(patient_row);
        plot_data.push(patient_expr_values);
      }

      // source code: https://plotly.com/javascript/heatmaps/
      // slice() to get only gene names
      let x_axis_label = Object.keys(result).slice(3);
      let y_axis_label = Object.values(result["Patient_ID"]);
      var information = [{
        z: plot_data,
        x: x_axis_label,
        y: y_axis_label,
        type: 'heatmap',
        hoverongraphs:false
      }];
      var layout = {
        title: 'Heatmap of gene expression level in different patients',
      };
      Plotly.newPlot(document.getElementById('data_plot'), information, layout);
    },
    error: function (e) {
      alert('Request Failed')
    }
  })
});