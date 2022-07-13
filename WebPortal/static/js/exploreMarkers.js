function HeatmapMarkerGenes(result,html_element_id,selected_cell) {
    // const start = performance.now();
    // console.log("start plotting the heatmap");
    // check the button id active

    if (html_element_id === "") {
        html_element_id = "displayPlotMarkers";
    }

    let genes = Object.keys(result[Object.keys(result)[0]])
    let celltypes = Object.keys(result);
    
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

        // let heatmap_width = 1300;
        let heatmap_height = 270 + 20 * ngenes;

        let data_content = [];
        for (var i = 0; i < ngenes; i++) {
            let each_gene_data = [];
            for (var j = 0; j < ncelltypes; j++) {
                exp = result[celltypes[j]][genes[i]];
                // exp = Math.log10(exp + 0.5);
                each_gene_data.push(exp);
            }
            data_content.push(each_gene_data);
        }

        var data = {
            type: 'heatmap',
            // hoverongaps: false,
            colorscale: 'Reds',
        };
        var layout = {
            autosize: true, 
            title: 'Expression profile of ' + selected_cell +'\'s marker genes in all cell types ',
            xaxis: {
                title: '<b>Cell types<b>',
                automargin: true,
                tickangle: 45
            },
            yaxis: {
                title: '<b>Genes<b>',
                automargin: true,
                autotick: false,
            },
            // with: heatmap_width,
            height: heatmap_height,
        };
        
        data['z'] = data_content;
        data['x'] = x_axis;
        data['y'] = y_ticks;
        Plotly.newPlot(document.getElementById(html_element_id), [data],layout);
    };
} 

function pagesetup() {
    $('#loadingtext').hide();
    $('#loadingbar').hide();
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/all_cell_types',
        success: function(result) {
            var celltype_categories = Object.keys(result);

            for(var i=0;i<celltype_categories.length;i++) {
                var category = celltype_categories[i];
                $("#celltype_selection").append(`<label class='radio has-text-weight-bold mt-2 has-text-success' disabled>${category}</label><br>`);
                var celltypes = result[category];
                for (var j=0;j<celltypes.length;j++) {
                    $("#celltype_selection").append(`<label class='radio'><input type='radio' name='celltype_selection' value='${celltypes[j]}'/>${celltypes[j]}</label><br>`);
                }
            }
        },
        error: function (e) {
            alert('Request data fail (no cell types available)')
        }
    });
}

function AssembleAjaxRequestMarker() {

    var selected_cell = $('input[name="celltype_selection"]:checked').val();

    if (selected_cell) {
    // Generate the plot
        $.ajax({
            type:'GET',
            url:'http://127.0.0.1:5000/markers_page',
            data: "celltype=" + selected_cell.replace('+','%2b'),
            success: function(result) {
                HeatmapMarkerGenes(result,'',selected_cell);
            },
            error: function (e) {
                alert('Request data fail (no cell types available)')
            }
        });
    }
}
$(document).ready(pagesetup)

function clearCheckbox() {
    $('input[name="celltype_selection"]:checked').each(function(){
        $(this).prop('checked', false);
    });
}

$("#clearOnClick").click(clearCheckbox)