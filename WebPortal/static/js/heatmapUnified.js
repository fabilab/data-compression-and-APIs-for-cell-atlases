function range(start, end) {
    var ans = [];
    for (let i = start; i < end; i++) {
        ans.push(i);
    }
    return ans;
}

function plotHeatmapUnified(gene_name, result) {
    var data = result['data'];
    data['type'] = 'scatter';
    data['marker']['symbol'] = 'square';

    var xticks = result['xticks'];
    var yticks = result['yticks'];
    var nx = xticks.length;
    var ny = yticks.length;

    var layout = {
        automargin: true,
        autosize: true,
        width: 1300,
        height: 1000,
        title: {
            text: gene_name + ' expression over time',
            x: 0.5,
            xanchor: 'center',
        },
        xaxis: {
            tickangle: 60,
            automargin: true,
            linewidth: 0,
            type: 'category',
        },
        yaxis: {
            //automargin: true,
            autorange: 'reversed',
            type: 'category',
            tickvals: result['yticks'],
            ticktext: result['yticktext'],
        },
    };

    Plotly.newPlot(
        document.getElementById('heatmap_unified'),
        [data],
        layout,
    ); 
}


function AssembleAjaxRequest() {

    // action here when clicking the search button
    var gene_name = $('#searchGeneName').val();
    if (gene_name === '') {
      gene_name = "Car4";
    }
    cpm_is_active = $("#cpmTab").hasClass('is-active');
    orginal_is_active = $("#originalOrderTab").hasClass('is-active')
    var use_log, use_hierarchical;
    
    if (cpm_is_active) {
      use_log = '0';
    } else {
      use_log = '1';
    } 
    
    if (orginal_is_active) {
      use_hierarchical = "0";
    } else {
      use_hierarchical = "1";
    }

    $.ajax({
        type:'GET',
        url:'/data_heatmap_unified',
        data: "gene=" + gene_name + "&use_log=" + use_log + "&use_hierarchical=" + use_hierarchical,
        dataType:'json',
        success: function(result) {
            plotHeatmapUnified(gene_name, result)
        },
        error: function (e) {
          alert('Request data Failed')
        }
    });
}

$("#searchOnClick").click(function() {
  // action here when clicking the search button
  AssembleAjaxRequest();
})

$(document).ready(function () {
    AssembleAjaxRequest();
});
