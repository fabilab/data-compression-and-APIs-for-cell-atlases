$(document).ready(function() {
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data',
        data: "gene_names=Car4,Vwf,Col1a1,Ptprc,Ms4a1&plot_type=original&data_type=original&data_set=celltype",
        dataType:'json',
        success: function(result) {
            HeatMap(result, "h5_data_plot_1");
        },
        error: function (e) {
        alert('Request data Failed')
        }
    });
    $("#originalTab").addClass('is-active');
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data',
        data: "gene_names=Car4,Vwf,Col1a1,Ptprc,Ms4a1&plot_type=original&data_type=original&data_set=celltype_dataset",
        dataType:'json',
        success: function(result) {
            HeatMap(result, "h5_data_plot_2");
        },
        error: function (e) {
        alert('Request data Failed')
        }
    });
    $.ajax({
        type:'GET',
        url:'http://127.0.0.1:5000/data',
        data: "gene_names=Car4,Vwf,Col1a1,Ptprc,Ms4a1&plot_type=original&data_type=original&data_set=celltype_dataset_timepoint",
        dataType:'json',
        success: function(result) {
            HeatMap(result, "h5_data_plot_3");
        },
        error: function (e) {
        alert('Request data Failed')
        }
    });
})