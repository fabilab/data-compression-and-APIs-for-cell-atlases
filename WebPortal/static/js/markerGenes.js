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
    var selected_id = 'B cell'
    $('input[type=checkbox]').each(function () {
        if (this.checked === true) {
            selected_id  = this.id;
        }
    });
    var selected_cell = selected_id.split('_')[0];


    // Generate the plot
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/markers_page',
        data: "celltype=" + selected_cell,
        success: function(result) {
            $('#loadingtext').show();
            $('#loadingbar').show();
            HeatmapMarkerGenes(result,'',selected_cell);
            $('#loadingtext').hide();
            $('#loadingbar').hide();
        },
        error: function (e) {
            alert('Request data fail (no cell types available)')
        }
    });
}
$("#applyOnClick" ).click(AssembleAjaxRequestMarker)
$(document).ready(AssembleAjaxRequestMarker)