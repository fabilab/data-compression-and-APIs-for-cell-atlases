var dataAverageExp = {};

function HeatmapAverageExp(result_wrapper, html_element_id) {

        let useLog = dataAverageExp['useLog'];

        if (result_wrapper === "") {
            result_wrapper = dataAverageExp['result_wrapper'];
        } else {
            dataAverageExp['result_wrapper'] = result_wrapper;
        }

        if (html_element_id === "") {
            html_element_id = "displayPlot";
        }

        let result = result_wrapper['result_average'];
        let celltypes;
        let genes = Object.keys(result);
        let celltypeOrder = dataAverageExp['celltypeOrder'];

        if (!celltypeOrder) {
            celltypes = Object.keys(result[Object.keys(result)[0]]);
        } else {
            celltypes = result_wrapper['hierarchicalCelltypeOrder'];
        }

        if (!result) {
            alert("Error:Input gene name is invalid, please make sure you type in the corrent gene names")
        } else {
            // x-axis: cell types
            let x_axis = celltypes;
            // y-axis: genes
            let y_axis = genes;
            let y_ticks = [];
            for (var i = 0; i < genes.length; i++) {
                y_ticks[i] = `<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${genes[i]}">${genes[i]}</a>`;
            }
            
            let ngenes = y_axis.length;
            let ncelltypes = x_axis.length;

            let heatmap_width = 1050;
            let heatmap_height = 270 + 25 * ngenes;

            let data_content = [];
            for (var i = 0; i < ngenes; i++) {
                let each_gene_data = [];
                for (var j = 0; j < ncelltypes; j++) {
                    exp = result[y_axis[i]][x_axis[j]];
                    if (useLog) {
                        exp = Math.log10(exp + 0.5);
                    }
                    each_gene_data.push(exp);
                }
                data_content.push(each_gene_data);
            }

            var data = {
                type: 'heatmap',
                hoverongaps: false,
                colorscale: 'Reds',
            };
            var layout = {
                autosize: true, 
                title: 'Expression profile of selected genes in all cell types',
                xaxis: {
                    title: '<b>Cell types<b>',
                    automargin: true,
                    tickangle: 45
                },
                yaxis: {
                    title: '<b>Genes<b>',
                    automargin: true
                },
                width: heatmap_width,
                height: heatmap_height,
            };
            
            if ($('#'+html_element_id).text() === "") {
                data['z'] = data_content;
                data['x'] = x_axis;
                data['y'] = y_ticks;
                Plotly.newPlot(document.getElementById(html_element_id), [data],layout);
            } else {
                data['z'] = [data_content];
                data['x'] = [x_axis];
                data['y'] = [y_ticks];
                Plotly.update(document.getElementById(html_element_id), data);
            }
        };
    } 

function DotplotProportionExp(result_wrapper, html_element_id) {
  let useLog = dataAverageExp['useLog'];
  if (result_wrapper === "") {
    result_wrapper = dataAverageExp['result_wrapper'];
  } else {
    dataAverageExp['result_wrapper'] = result_wrapper;
  }

  
  if (html_element_id === "") {
    html_element_id = "displayPlot";
  }

  let result_avg = result_wrapper['result_average'];
  let result_proportion = result_wrapper['result_proportion'];
  let celltypes;
  let genes = Object.keys(result_avg);
  let celltypeOrder = dataAverageExp['celltypeOrder'];

  if (!celltypeOrder) {
    celltypes = Object.keys(result_avg[Object.keys(result_avg)[0]]);
  } else {
    celltypes = result_wrapper['hierarchicalCelltypeOrder'];
  }

  // console.log(genes);
  // console.log(celltypes);

  let ngenes = genes.length;
  let ncelltypes = celltypes.length;

  let data = [];
  var desired_maximum_marker_size = 8;
  let max_expr_value = result_wrapper['max_expression'];
  console.log(max_expr_value);
  for (var i = 0; i < ngenes; i++) {
    let this_gene = genes[i];
    
    let marker_color = [];
    let marker_size = [];
    let hovertext = [];
    for (var j = 0; j < ncelltypes; j++) {
      marker_size.push(result_proportion[this_gene][celltypes[j]]*100);
      
      let exp = result_avg[this_gene][celltypes[j]] / max_expr_value;
      if (useLog) {
          exp = Math.log10(exp + 0.5);
      }
      marker_color.push(exp);
      hovertext.push(`Gene: ${this_gene}<br>Celltype: ${celltypes[j]}<br>AvgExp: ${exp.toPrecision(3)}<br>Proportion: ${result_proportion[this_gene][celltypes[j]].toPrecision(3)*100}%`);
    }

    let gene_trace = {
      x: celltypes,
      y: Array(ncelltypes).fill(this_gene),
      mode: 'markers',
      marker: {
        size: marker_size,
        sizeref: 2 * Math.max(...marker_size) / (desired_maximum_marker_size**2),
        color: marker_color,
        // colorscale: 'YlGnBu',
      },
      text: hovertext,
      hoverinfo: 'text',
    }
    if (i === 0) {
      gene_trace['marker']['colorbar'] = {thickness:15}
    }
    data.push(gene_trace);
  }
  
  var layout = {
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
      width: 1000,
      height: 400 + 25 * ngenes,
  };
  
  Plotly.newPlot(document.getElementById(html_element_id), data,layout);

}

function AjaxExploreGeneral() {

  if(! $('#scatter_plot').is('empty')) {
    $('#scatter_plot').empty();
  }

  // When doing the search gene name action, we want it to be change immediatly without switching back to the original heatmap,
  //  for example, if we are looking at a log10 plot,and we do the search action, the tab stays at the log10 
  cpm_is_active = $("#cpmTab").hasClass('is-active');
  orginal_is_active = $("#originalOrderTab").hasClass('is-active')
  var plot_type = 'original';
  var data_type = 'original';
  
  if (!cpm_is_active) {
    data_type = 'log10';
  } 
  
  if (!orginal_is_active) {
    plot_type = "hieracical";
  }
  
  // action here when clicking the search button
  var genes_string = $('#listGenes').val();
  const gene_array = genes_string.split(",")
  if (gene_array.length == 2) {
    $.ajax({
      type:'GET',
      url:'http://127.0.0.1:5000/data_scatter',
      data: "gene_names=" + genes_string,
      success: ScatterPlot,
      error: function (e) {
        alert('Request data Failed')
      }
    });
  }
    // sent gene names to the API
  $.ajax({
    type:'GET',
    url:'http://127.0.0.1:5000/data_general',
    data: "gene_names=" + genes_string + "&plot_type=" + plot_type + "&data_type=" + data_type,
    success: function(result) { 
      $("#displayPlot").empty();
      // HeatmapAverageExp(result, "displayPlot");
      DotplotProportionExp(result, "displayPlot");
    },
    error: function (e) {
      alert('Error:Input gene name is invalid, please make sure you type in the corrent gene names.')
    }
    });
  }
$("#searchOnClick_list" ).click(AjaxExploreGeneral)