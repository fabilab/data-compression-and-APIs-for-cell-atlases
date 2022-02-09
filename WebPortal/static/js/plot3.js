$(document).ready(function() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data2',
        dataType:'json',
        success: function (result) {
            // compare the expression level of 2 specific genes across all cell types
            let expr1 = result["LDHB"];
            let expr2 = result["HLA-DRA"];
            let cell_types = Object.keys(result[Object.keys(result)[0]]);
            var trace1 = {
                x: Object.values(expr1),
                y: Object.values(expr2),
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
                        text: 'LDHB expression',
                        font: {
                            family: 'Courier New, monospace',
                            size: 18,
                            color: '#7f7f7f'
                        }
                    },
                    range: [ 0, 3 ]
                },
                yaxis: {
                    title: {
                        text: 'HLA-DRA expression',
                        font: {
                            family: 'Courier New, monospace',
                            size: 18,
                            color: '#7f7f7f'
                        }
                    },
                    range: [0, 6]
                },
                title:'Average expression level of LDHB and HLA-DRA in different cell types'
            };
            
            Plotly.newPlot(document.getElementById('data_3_plot'), data, layout);
        },
        error: function (e) {
        alert('Request data2 Failed')
        }
    });
})
