$(document).ready(function() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data2',
        dataType:'json',
        success: function (result) {
            let x_axis_label = Object.keys(result);
            let y_axis_label = Object.keys(result[Object.keys(result)[0]]);
            let z_data_all = [];
            for (var i = 0; i < Object.keys(result).length; i++) {
                // get the current gene name, ie: LDHB
                let gene = Object.keys(result)[i]
                // access avergae expression level of current gene in all 8 cell types
                let gene_expr_all_cell_type = result[gene]
                let z_data_each = []
                for (var j = 0; j < Object.keys(gene_expr_all_cell_type).length; j++) {
                    let current_cell_type = Object.keys(gene_expr_all_cell_type)[j];
                    z_data_each.push(gene_expr_all_cell_type[current_cell_type]);
                }
                z_data_all.push(z_data_each);
            }
            var data = [
                {
                    z: z_data_all,
                    x: x_axis_label,
                    y: y_axis_label,
                    type: 'heatmap',
                    hoverongaps: false
                }
              ];
            var layout = {
                title: 'Heatmap of Average gene expression level in different cell types',
            };
              
              Plotly.newPlot(document.getElementById('data_2_plot'), data,layout);  
        },
        error: function (e) {
        alert('Request data2 Failed')
        }
    });
})
