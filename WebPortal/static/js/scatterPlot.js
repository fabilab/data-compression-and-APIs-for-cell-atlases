function ScatterPlot(result) {
    // compare the expression level of 2 specific genes across all cell types
    let expr1 = result["gene1_expr"];
    let expr2 = result["gene2_expr"];
    // let expr1_log = expr1.map(Math.log10);
    // let expr2_log = expr2.map(Math.log10);
    let cell_types =result["cell_types"];
    var trace1 = {
        x: expr1,
        y: expr2,
        mode: 'markers',
        type: 'scatter',
        text: cell_types,
        // textposition: 'bottom center',
        textfont: {
          family:  'Raleway, sans-serif'
        },
        marker: { size: 12 }
    };
        
    var data = [ trace1 ];
    var x_max = Math.max(expr1) + 1;
    var x_min = Math.min(expr1) - 1;
    var y_max = Math.max(expr2) + 1;
    var y_min = Math.min(expr2) - 1;
    
    var layout = {
        xaxis: {
            title: {
              text: result['gene1_name']+' expression',
            },
            range: [ x_min, x_max ]
        },
        yaxis: {
            title: {
              text: result['gene2_name']+' expression',
          },
          range: [y_min, y_max]
        },
        title:'Comparison of gene expression level between ' + result['gene1_name'] + ' and ' + result['gene2_name'] + ' (Log10 Applied)',
        height: 600
    };
    
    Plotly.newPlot(document.getElementById('scatterPlot'), data, layout);
}