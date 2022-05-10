// action here when clicking the search button
var gene_name = $('#searchGeneName').val();

$(document).ready(function() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data_timepoint',
        data: "gene=" + gene_name,
        dataType:'json',
        success: function(result) {
            num = 1
            for (key in Object.keys(result)) {
                dataset = Object.keys(result)[key]
                let div_id = 'dataset_' + num
                HeatMapTimepoint(result[dataset], div_id,dataset);
                num++
            }
        },
        error: function (e) {
        alert('Request data Failed')
        }
    });
    $("#originalTab").addClass('is-active');
})