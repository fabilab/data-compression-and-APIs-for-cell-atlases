var gene_name = $('#searchGeneName').val();
if (gene_name === '') {
  gene_name = "Car4";
}

$(document).ready(function() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data_timepoint_dataset',
        data: "gene=" + gene_name,
        dataType:'json',
        success: function(result) {
            HeatMap(result, "bigHeatMap");
        },
        error: function (e) {
        alert('Request data Failed')
        }
    });
    $("#originalTab").addClass('is-active');
})