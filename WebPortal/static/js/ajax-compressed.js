var base_url = window.location.href.replace('/dataExplore', '');

var dataAverageExp = {};

function generateHmCompressed(result_wrapper, html_element_id) {
        $("#hm_compressed").empty();
        let useLog = dataAverageExp['useLog'];

        if (result_wrapper === "") {
            result_wrapper = dataAverageExp['result_wrapper'];
        } else {
            dataAverageExp['result_wrapper'] = result_wrapper;
        }

        if (html_element_id === "") {
            html_element_id = "hm_compressed";
        }
        let result = result_wrapper['result_average'];
        let celltypes;
        let genes;
        let celltypeOrder = dataAverageExp['celltypeOrder'];

        if (!celltypeOrder) {
            celltypes = Object.keys(result[Object.keys(result)[0]]);
            genes = Object.keys(result);
        } else {
            celltypes = result_wrapper['hierarchicalCelltypeOrder'];
            genes = result_wrapper['hierarchicalGeneOrder'];
        }

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
        $("#num_genes").empty();
        $("#num_genes").append(ngenes);

        let heatmap_width = 1050;
        let heatmap_height = 280 + 25 * ngenes;

        let data_content = [];
        let hover_text = [];
        for (var i = 0; i < ngenes; i++) {
            let each_gene_data = [];
            let temp = [];
            for (var j = 0; j < ncelltypes; j++) {
                exp = result[y_axis[i]][x_axis[j]];
                let gn = genes[i]
                let ct = celltypes[j];
                if (useLog) {
                    exp = Math.log10(exp + 0.5);
                }
                each_gene_data.push(exp);
                temp.push('Celltype: '+ct+'<br>Gene: '+gn+'<br>Expression: '+exp);
            }
            data_content.push(each_gene_data);
            hover_text.push(temp);
        }

        var data = {
            type: 'heatmap',
            hoverongaps: false,
            colorscale: 'Reds',
            text:hover_text,
            hoverinfo: 'text',
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
        // source code from: https://codepen.io/etpinard/pen/zzzBXv?editors=0010
        var tools = {
            modeBarButtonsToAdd: [
                {
                    name: 'Download plot as an SVG',
                    icon: Plotly.Icons.camera,
                    click: function(gd) {
                    Plotly.downloadImage(gd, {format: 'svg'})
                    }
                },
                {
                    name: "downloadCsv",
                    title: "Download data as csv",
                    icon: Plotly.Icons.disk,
                    click: (gd) => {
                      let data = [
                        [gd.layout.xaxis.title.text, gd.layout.yaxis.title.text].join(
                          ","
                        ),
                      ];
                      gd.data[0].x.forEach((xvalue, i) =>
                        data.push([xvalue, gd.data[0].y[i]].join(","))
                      );
                      let blob = new Blob([data.join("\r\n")], { type: "text/csv" });
                      // import { saveAs } from "file-saver";
                      saveAs(
                        blob,
                        "export.csv"
                      );
                    },
                  }
            
            
            ],
              editable:true,
              responsive: true,
        }

        // var config = {responsive:true}
        
        if ($('#'+html_element_id).text() === "") {
            data['z'] = data_content;
            data['x'] = x_axis;
            data['y'] = y_ticks;
            Plotly.newPlot(document.getElementById(html_element_id), [data],layout,tools);
        } else {
            data['z'] = [data_content];
            data['x'] = [x_axis];
            data['y'] = [y_ticks];
            Plotly.update(document.getElementById(html_element_id),data,tools);
        }
    } 


function AjaxCompressed() {

  if(! $('#scatter_plot').is('empty')) {
    $('#scatter_plot').empty();
  }

  var genes_string = $('#listGenes').val();
  const gene_array = genes_string.split(",")
  if (gene_array.length == 2) {
    $.ajax({
      type:'GET',
      url:`${base_url}/data_scatter`,
      data: "gene_names=" + genes_string,
      success: ScatterPlot,
      error: function (e) {
        Swal.fire('Invalid input', `${e.responseText.substring(1, e.responseText.length - 2)}is invalid, please make sure you type in the corrent gene names.`, 'error')
      }
    });
  }
  $.ajax({
    type:'GET',
    url:`${base_url}/data_general`,
    data: "gene_names=" + genes_string,
    success: function(result) { 
      $("#hm_compressed").empty();
      generateHmCompressed(result, "hm_compressed");
      $("#dp_compressed").empty();
      generateDpCompressed(result, "dp_compressed");
    },
    error: function (e) {
        Swal.fire('Invalid input', `${e.responseText.substring(1, e.responseText.length - 2)}is invalid, please make sure you type in the correct gene names.`, 'error')
    }
    });
  }

function clearTextArea() {
    $("#listGenes").val('');
}
$("#searchOnClick_list").click(AjaxCompressed)
$("#clearOnClick_list").click(clearTextArea)