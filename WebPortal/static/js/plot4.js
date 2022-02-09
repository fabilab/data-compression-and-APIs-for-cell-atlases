$(document).ready(function() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data3',
        dataType:'json',
        success: function (result) {
            let all_genes = Object.keys(result);
            gene_data = []
            for (var i = 0; i < Object.keys(result).length; i++) {
                // get the current gene name, ie: LDHB
                let gene = Object.keys(result)[i]
                gene_data.push(result[gene])
            }
            var data = [ {
                type: 'vioin',
                x: all_genes,
                y: gene_data,
                box: {
                    visible: true
                },
                meanline: {
                visible: true
                }

            }]
            var layout = {
                title: 'Vilon plot of expression level of all genes across B cell',
            }
            Plotly.newPlot(document.getElementById('data_4_plot'), data, layout);
        },
        error: function (e) {
        alert('Request data3 Failed')
        }
    });
})
