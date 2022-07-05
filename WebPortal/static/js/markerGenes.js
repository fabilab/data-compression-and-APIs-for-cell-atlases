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
                $("#celltype_selection").append(`<label class='radio has-text-weight-bold' disabled>===== ${category} =====</label><br>`);
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

    $('#loadingtext').hide();
    $('#loadingbar').hide();
    if (selected_cell) {
        $('#markerWelcome').hide();
        $('#loadingtext').show();
        $('#loadingbar').show();
    // Generate the plot
        $.ajax({
            type:'GET',
            url:'http://127.0.0.1:5000/markers_page',
            data: "celltype=" + selected_cell.replace('+','%2b'),
            success: function(result) {
                $('#loadingtext').hide();
                $('#loadingbar').hide();
                HeatmapMarkerGenes(result,'',selected_cell);
            },
            error: function (e) {
                alert('Request data fail (no cell types available)')
            }
        });
    }
}
$("#applyOnClick").click(AssembleAjaxRequestMarker)
$(document).ready(pagesetup)

function clearCheckbox() {
    $('input[name="celltype_selection"]:checked').each(function(){
        $(this).prop('checked', false);
    });
}

$("#clearOnClick").click(clearCheckbox)