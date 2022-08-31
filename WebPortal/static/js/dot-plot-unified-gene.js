var dataForPlotsUnified = {};

function generateDpUnifiedGene(result,html_element_id) {
    // flags passed by events.js
    let useLog = dataForPlotsUnified['useLog'];
    let celltypeOrder = dataForPlotsUnified['celltypeOrder'];
    
    if (result === "") {
        result = dataForPlotsUnified['result'];
    } else {
        dataForPlotsUnified['result'] = result;
    }

    if (html_element_id === "") {
        html_element_id = "dp_unified_gene";
    }

    let celltypes;
    if (!celltypeOrder) {
        celltypes = result['cell_type'];
    } else {
        celltypes = result['hierarchicalCelltypeOrder'];
    }

    let result_avg = result['exp_avg'];
    let result_pro = result['exp_pro'];
    
    // x-axis: celltypes
    let all_x = [];
    //y-axis: timepoint
    let all_y = [];
    // expression level
    let all_color = [];
    // proportion expression
    let all_size = [];
    //dataset
    let y_axis_dataset = [];
    // timepoint that with duplicate:
    let y_axis_time = [];
    var desired_maximum_marker_size = 6.2;
    
    for (var i = 0; i < result['dataset_timepoint'].length; i++) {
        // dataset timepoint combination: e.g TMS_24m
        let dataset_timepoint = result['dataset_timepoint'][i];
        let time  = dataset_timepoint.split("_")[1];
        let dataset = dataset_timepoint.split("_")[0];
        y_axis_time.push(time);
        y_axis_dataset.push(dataset);
        
        for (var k = 0; k < celltypes.length; k++) {
            if (time === 'P7' && dataset === 'ACZ') {
                all_y.push(time+'\'\'');
            } else {
                all_y.push(time);
            }
            all_x.push(celltypes[k]);
            if (result_avg[dataset_timepoint][celltypes[k]] === -1) {
                all_color.push(null);
                all_size.push(null);
            } else {
                if (useLog) {
                    all_color.push(Math.log10(result_avg[dataset_timepoint][celltypes[k]] + 0.5));
                } else {
                    all_color.push(result_avg[dataset_timepoint][celltypes[k]]);
                }
                all_size.push(result_pro[dataset_timepoint][celltypes[k]]*100);
            }
        }
    }
    // Generate hover text for each spot
    // Celltype: {ct}, Expression: {exp}, Dataset: {ds}, Timepoint: {tp}, 
    let hover_text = [];
    for (var i = 0; i < result['dataset_timepoint'].length; i++) {
        let temp = [];
        for (var j = 0; j < result['cell_type'].length; j++) {
            let dt = result['dataset_timepoint'][i];
            let ct = result['cell_type'][j];
            let exp;
            let pro = result_pro[dt][ct]*100+'%'
            if(useLog && result_avg[dt][ct] !== null) {
                exp = Math.log10(result_avg[dt][ct]+0.5);
            } else{
                exp = result_avg[dt][ct];
            }
            // "ACZ_P21"
            let ds = dt.split("_")[0];
            let tp = dt.split("_")[1];
            hover_text.push('Celltype: '+ct+'<br>Dataset: '+ds+'<br>Timepoint: '+tp+'<br>Expression: '+exp+'<br>Proportion: '+ pro);
        }
        // hover_text.push(temp);
    }
  
    let nTimepoints = y_axis_time.length;
    let data = {
        mode:'markers',
        marker: {
            color: all_color,
            size: all_size,
            sizeref: 2 * Math.max(...all_size) / (desired_maximum_marker_size**2),
            colorscale: 'YlGnBu',
            reversescale:true,
            colorbar: {},
        },
        hoverinfo:'text',
    };
    
    var layout = {
        title: 'Proportion expression profile of <b>' + result['gene'] + '</b> gene in all cell types over development time',
        autosize: false,
        showlegend:false,
        xaxis: {
            title: '<b>Cell types<b>',
            automargin: true,
            tickangle: 45,
        },
        yaxis: {
            title: '<b>Dataset Timepoints<b>',
            automargin: true,
        },
        width: 1000,
        height: 400+25*nTimepoints,
    };
    var tools = {
        modeBarButtonsToAdd: [{
            name: 'Download plot as an SVG',
            icon: Plotly.Icons.camera,
            click: function(gd) {
              Plotly.downloadImage(gd, {format: 'svg'})
            }
          }],
          editable:true,
    };
    // Plotly.newPlot(document.getElementById(html_element_id), [data],layout);
    if ($('#'+html_element_id).text() === "") {
        data['x'] = all_x;
        data['y'] = all_y;
        data['text'] = hover_text;
        Plotly.newPlot(document.getElementById(html_element_id), [data],layout,tools);   
    } else {
        data['x'] = [all_x];
        data['y'] = [all_y];
        data['text'] = [hover_text];
        Plotly.update(
            document.getElementById(html_element_id), 
            data, 
            {
                yaxis: {
                    title: '<b>Timepoints<b>', 
                    automargin: true,
                }
            },
        layout,
        tools);
    }};