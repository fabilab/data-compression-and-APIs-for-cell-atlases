function AssembleAjaxRequestUnified() {
    cpm_is_active = $("#cpmTab").hasClass('is-active');
    orginal_is_active = $("#originalOrderTab").hasClass('is-active')
    var plot_type = 'original';
    var data_type = 'original';

    if (!cpm_is_active) {
        data_type = 'log10';
    } 

    if (!orginal_is_active) {
        plot_type = "hieracical";
    }
    // action here when clicking the search button
    var gene_name = $('#searchGeneName').val();
    $(document).ready(function() {
        $.ajax({
            type:'GET',
            url:'http://127.0.0.1:5000/unified',
            data: "gene=" + gene_name + "&plottype=" + plot_type + "&datatype=" + data_type,
            dataType:'json',
            success: function(result) {
                // HeatmapCelltype(result, "bigHeatMap");
                console.log(result)
            },
            error: function (e) {
                alert('Request data Failed !')
            }
        });
        $("#originalTab").addClass('is-active');
    });
}
// $("#searchOnClick" ).click(AssembleAjaxRequestTimepoint)
$(document).ready(AssembleAjaxRequestUnified)