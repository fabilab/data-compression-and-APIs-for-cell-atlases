function AssembleAjaxRequest() {
  var command = $('#textCommand').val();

  $.ajax({
      type: 'GET',
      url: '/submit_text',
      data: "text="+command,
      dataType: 'json',
      success: function(result) {
          console.log(result);
      },
      error: function(e) {
          console.log(e);
          alert('Text command not understood');
      }
  })

}


$("#askOnClick").click(function() {
  // action here when clicking the search button
  AssembleAjaxRequest();
});
