dataForPlotsUnifiedCell = {}

function generateDpUnifiedCell(result,htmlElementId) {
    // flags passed by events.js    
    let useLog = dataForPlotsUnifiedCell['useLog'];
    let geneOrder = dataForPlotsUnifiedCell['geneOrder'];

    if (result === "") {
        result = dataForPlotsUnifiedCell['result'];
    } else {
        dataForPlotsUnifiedCell['result'] = result;
    }

    if (htmlElementId === "") {
        htmlElementId = "dp_unified_cell";
    }

    let genes;
    if (!geneOrder) {
        genes = result['genes'];
    } else {
        genes = result['clustered_gene_order'];
    }
    const exp_avg = result['exp_avg'];
    const exp_pro = result['exp_pro'];

    // x-axis: genes
    let x_axis = [];
    //y-axis: timepoint
    let y_axis = [];
    // expression level
    let color = [];
    // proportion expression
    let size = [];
    //dataset
    let y_axis_dataset = [];
    // timepoint that with duplicate:
    let y_axis_time = [];
    var desired_maximum_marker_size = 6.2;

    for (var i = 0; i < result['dataset_timepoint'].length; i++) {
        // dataset timepoint combination: e.g TMS_24m
        let dataset_timepoint = result['dataset_timepoint'][i];
        let dataset = dataset_timepoint.split("_")[0];
        let time  = dataset_timepoint.split("_")[1];
        y_axis_time.push(time);
        y_axis_dataset.push(dataset);
        
        for (var k = 0; k < genes.length; k++) {
            if (time === 'P7' && dataset === 'ACZ') {
                y_axis.push(time+'\'\'');
            } else {
                y_axis.push(time);
            }
            x_axis.push(genes[k]);
            if (exp_avg[dataset_timepoint][genes[k]] === -1) {
                color.push(null);
                size.push(null);
            } else {
                if (useLog) {
                    color.push(Math.log10(exp_avg[dataset_timepoint][genes[k]] + 0.5));
                } else {
                    color.push(exp_avg[dataset_timepoint][genes[k]]);
                }
                size.push(exp_pro[dataset_timepoint][genes[k]]*100);
            }
        }
    }

    // Generate hover text for each spot
    // Celltype: {ct}, Expression: {exp}, Dataset: {ds}, Timepoint: {tp}, 
    let hover_text = [];
    for (var i = 0; i < result['dataset_timepoint'].length; i++) {
        for (var j = 0; j < genes.length; j++) {
            let dt = result['dataset_timepoint'][i];
            let gn = genes[j];
            let exp;
            let pro = exp_pro[dt][gn]*100+'%'
            if(useLog && exp_avg[dt][gn] !== null) {
                exp = Math.log10(exp_avg[dt][gn]+0.5);
            } else{
                exp = exp_avg[dt][gn];
            }
            // "ACZ_P21"
            let ds = dt.split("_")[0];
            let tp = dt.split("_")[1];
            hover_text.push('Gene: '+gn+'<br>Dataset: '+ds+'<br>Timepoint: '+tp+'<br>Expression: '+exp+'<br>Proportion: '+ pro);
        }
    }
    let nTimepoints = y_axis_time.length;
    let data = {
        mode:'markers',
        marker: {
            color: color,
            size: size,
            sizeref: 2 * Math.max(...size) / (desired_maximum_marker_size**2),
            colorscale: 'YlGnBu',
            reversescale:true,
            colorbar: {},
        },
        hoverinfo:'text',
    };
    
    var layout = {
        title: 'Proportion expression profile of selected  Genes in <b>' + result['cell_type'] + '</b> over development time',
        autosize: false,
        showlegend:false,
        xaxis: {
            title: '<b>Genes<b>',
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
          responsive: true,
    };
    if ($('#'+htmlElementId).text() === "") {
        data['x'] = x_axis;
        data['y'] = y_axis;
        data['text'] = hover_text;
        Plotly.newPlot(document.getElementById(htmlElementId), [data],layout,tools);   
    } else {
        data['x'] = [x_axis];
        data['y'] = [y_axis];
        data['text'] = [hover_text];
        Plotly.update(
            document.getElementById(htmlElementId), 
            data, 
            {
                yaxis: {
                    title: '<b>Timepoints<b>', 
                    automargin: true,
                }
            },
        layout,
        tools);
    }
};