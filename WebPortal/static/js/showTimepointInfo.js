$.ajax({
    url: "/info.txt",
    dataType:"text",
    success: function(data){
    $("#data-container").text(data);
    }
});