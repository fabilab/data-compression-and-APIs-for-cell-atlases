$(document).ready(function() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data',
        dataType:'json',
        success: function (result) {
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
                title: 'Heatmap of gene expression level in different cell types',
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
})
