$(document).ready(function() {
    $('#searchGeneName').keyup(function(){
        $('#matchNames').html('');  // Clear the previou list
        var searchField = $('#searchGeneName').val();
        var expression = new RegExp(searchField,'i');
        $.ajax({
            type:'GET',
            url:'http://127.0.0.1:5000/all_gene_names',
            success: function(result) {
                // for earch geneName in result
                $.each(result, function(index){
                    var geneName = result[index];
                    if(geneName.search(expression) != -1) {
                        console.log(geneName);
                        $('#matchNames').append('<option value="' + geneName + '">' + geneName + '</option>')
                    }
                });
            }
          });
    });
});