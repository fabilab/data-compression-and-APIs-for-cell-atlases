function AssembleAjaxRequestMarker() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/markers_page',
        data: "celltype=" + 'Car4%2b capillaries',
        success: function(result) {
            console.log(result)
            // for(var i=0;i<result.length;i++) {
            //     var celltype = result[i];
            //     // var $input = $("<input>", {id:`${celltype}_check`,type:'checkbox'});
            //     var $label = $("<div>",{class:'block mgb-small',style:"margin-bottom: 0.0rem"});
            //     $label.append(`<input id='${celltype}_check' type='checkbox'> ${celltype}`);

            //     $("#celltypeFilter").append($label);
            // }

            HeatmapMarkerGenes(result,'');
        },
        error: function (e) {
            alert('Request data fail (no cell types available)')
        }
        });
}

$(document).ready(AssembleAjaxRequestMarker)