var dataForPlot = {};

function DotplotProportionExpMarker(result_original,result_scaled,exp_proportion,html_element_id,selected_cell, gene_order, showNumMarkers) {

    if (html_element_id === "") {
        html_element_id = "dotPlotMarker";
    }

    if (result_original === '') {
        result_original = dataForPlot['result_original']; 
    } else {
        dataForPlot['result_original'] = result_original; 
    }

    if (result_scaled === '') {
        result_scaled = dataForPlot['result_scaled']; 
    } else {
        dataForPlot['result_scaled'] = result_scaled; 
    }

    if (exp_proportion === '') {
        exp_proportion = dataForPlot['exp_proportion']; 
    } else {
        dataForPlot['exp_proportion'] = exp_proportion; 
    }

    if (selected_cell === '') {
        selected_cell = dataForPlot['selected_cell']; 
    } else {
        dataForPlot['selected_cell'] = selected_cell; 
    }

    if (gene_order === '') {
        gene_order = dataForPlot['gene_order']; 
    } else {
        dataForPlot['gene_order'] = gene_order; 
    }

    //  show the top num markers selected by the user, show all as default
    if (showNumMarkers === "") {
        showNumMarkers = gene_order.length
    }

    let exp_scaled = dataForPlot['result_scaled'];
    let result_proportion = dataForPlot['exp_proportion'];
    let celltypes = Object.keys(result_scaled);
    let genes = gene_order;
  
    let ngenes = genes.length;
    let ncelltypes = celltypes.length;
  
    let all_x = [];
    let all_y = [];
    let all_color = [];
    let all_size = [];
    let all_hovertext = [];
    var desired_maximum_marker_size = 6.5;

    for (var i = showNumMarkers-1; i >= 0; i--) {
      let this_gene = genes[i];
      for (var j = 0; j < ncelltypes; j++) {
        
        let exp = exp_scaled[celltypes[j]][this_gene];

        // marker_color.push(exp);
        all_x.push(celltypes[j])
        all_hovertext.push(`Gene: ${this_gene}<br>Cell Type: ${celltypes[j]}<br>Avg Exp: ${exp.toPrecision(3)}<br>Proportion: ${result_proportion[celltypes[j]][this_gene].toPrecision(3)*100}%`);
        all_y.push(`<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${genes[i]}">${genes[i]}</a>`);
        all_color.push(exp);
        all_size.push(result_proportion[celltypes[j]][this_gene]*100);
      }
    }

    let data = {
      mode: 'markers',
      marker: {
        color: all_color,
        size: all_size,
        sizeref: 2 * Math.max(...all_size) / (desired_maximum_marker_size**2),
        colorscale: 'YlGnBu',
        reversescale:true,
        colorbar: {},
      },
      hoverinfo: 'text',
    };
    
    var layout = {
        title:"Dot Plot showing expression level and population faction of <b>"+selected_cell+"'s</b> marker genes across 41 different cell types",
        autosize: false,
        showlegend: false,
        xaxis: {
            title: '<b>Cell types<b>',
            tickangle: 45,
            automargin: true
        },
        yaxis: {
            title: '<b>Genes<b>', 
            automargin: true,
            range: [-1, showNumMarkers]
        },
        width: 1000,
        height: 25 * showNumMarkers
    };
    // Plotly.newPlot(document.getElementById(html_element_id), [data],layout);
    if ($('#'+html_element_id).text() === "") {
        data['x'] = all_x;
        data['y'] = all_y;
        data['text'] = all_hovertext;
        Plotly.newPlot(document.getElementById(html_element_id), [data],layout);   
    } else {
        data['x'] = [all_x];
        data['y'] = [all_y];
        data['text'] = [all_hovertext];
        Plotly.update(
            document.getElementById(html_element_id), 
            data, 
            {
                height: 400 + 25 * showNumMarkers,
                yaxis: {
                    title: '<b>Genes<b>', 
                    automargin: true,
                    range: [-1, showNumMarkers]
                }
            },
        layout);
    }
  
  }