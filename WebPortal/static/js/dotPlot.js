var dataProportionExp = {};
function DotplotProportionExp(result_wrapper, html_element_id) {
    let useLog = dataProportionExp['useLog'];
    if (result_wrapper === "") {
      result_wrapper = dataProportionExp['result_wrapper'];
    } else {
      dataProportionExp['result_wrapper'] = result_wrapper;
    }
  
    if (html_element_id === "") {
      html_element_id = "dotPlot";
    }
  
    let result_avg = result_wrapper['result_average'];
    let result_proportion = result_wrapper['result_proportion'];
    let celltypes;
    let genes = Object.keys(result_avg);
    let celltypeOrder = dataProportionExp['celltypeOrder'];
  
    if (!celltypeOrder) {
      celltypes = Object.keys(result_avg[Object.keys(result_avg)[0]]);
    } else {
      celltypes = result_wrapper['hierarchicalCelltypeOrder'];
    }
  
    let ngenes = genes.length;
    let ncelltypes = celltypes.length;
  
    let all_x = [];
    let all_y = [];
    let all_color = [];
    let all_size = [];
    let all_hovertext = [];
    var desired_maximum_marker_size = 6.3;


    for (var i = 0; i < ngenes; i++) {
      let this_gene = genes[i];
      
      for (var j = 0; j < ncelltypes; j++) {
        
        let exp = result_avg[this_gene][celltypes[j]];
        if (useLog) {
            exp = Math.log10(exp + 0.5);
        }
        // marker_color.push(exp);
        all_x.push(celltypes[j])
        all_hovertext.push(`Gene: ${this_gene}<br>Cell Type: ${celltypes[j]}<br>Avg Exp: ${exp.toPrecision(3)}<br>Proportion: ${result_proportion[this_gene][celltypes[j]].toPrecision(3)*100}%`);
        all_y.push(`<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${genes[i]}">${genes[i]}</a>`);
        all_color.push(exp);
        all_size.push(result_proportion[this_gene][celltypes[j]]*100);
      }
    }

    let data = {
      x: all_x,
      y: all_y,
      mode: 'markers',
      marker: {
        color: all_color,
        size: all_size,
        sizeref: 2 * Math.max(...all_size) / (desired_maximum_marker_size**2),
        colorscale: 'YlGnBu',
        reversescale:true,
        colorbar: {},
      },
      text: all_hovertext,
      hoverinfo: 'text',
    };
    
    var layout = {
        title:"Dot Plot showing expression level and population faction of selected genes across 41 different cell types",
        autosize: true,
        showlegend: false,
        xaxis: {
            title: '<b>Cell types<b>',
            automargin: true,
            tickangle: 45
        },
        yaxis: {
            title: '<b>Genes<b>',
            automargin: true
        },
        width: 1050,
        height: 400 + 25 * ngenes,
    };
    
    Plotly.newPlot(document.getElementById(html_element_id), [data],layout);
  
  }