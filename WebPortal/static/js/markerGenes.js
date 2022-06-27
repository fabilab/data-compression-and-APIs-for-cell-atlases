function AssembleAjaxRequestMarker() {

    // Load a list of cell types check box
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/all_cell_types',
        success: function(result) {
            for(var i=0;i<result.length;i++) {
                var celltype = result[i];
                // var $input = $("<input>", {id:`${celltype}_check`,type:'checkbox'});
                var $label = $("<div>",{class:'block mgb-small',style:"margin-bottom: 0.0rem"});
                $label.append(`<input id='${celltype}_check' type='checkbox'> ${celltype}`);

                $("#celltypeFilter").append($label);
            }
        },
        error: function (e) {
            alert('Request data fail (no cell types available)')
        }
    });

    // update the plot when a new cell type from the checkbox is selected 
    var selected_id = ''
    $('input[type=checkbox]').each(function () {
        if (this.checked === true) {
            selected_id  = this.id;
        }
    });
    var selected_cell = selected_id.split('_')[0];

    $('#loadingtext').hide();
    $('#loadingbar').hide();
    if (selected_id !== '') {
        $('#markerWelcome').hide();
        $('#loadingtext').show();
        $('#loadingbar').show();
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
        $('#loadingtext').hide();
        $('#loadingbar').hide();
    }
}
$("#applyOnClick").click(AssembleAjaxRequestMarker)
$(document).ready(AssembleAjaxRequestMarker)

function clearCheckbox() {
    $('input[type=checkbox]').prop('checked', false);
}

$("#clearOnClick").click(clearCheckbox)