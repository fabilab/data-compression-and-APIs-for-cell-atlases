let chunks = [];
var rec;
var stream;

//add events to the button
$("#recordButton").mousedown(() => {
    $(window).off();
    $(window).mouseup(stopRecording);
    console.log("recordButton clicked");

    // Switch to stop button until we hear back from getUserMedia
    $("#recordButton").attr("src", "/static/images/stop_button.png")
                      .attr("alt", "stop button");

    /*
        Simple constraints object, for more advanced audio features see
        https://addpipe.com/blog/audio-constraints-getusermedia/
    */
    var constraints = { audio: true, video:false };
    navigator.mediaDevices.getUserMedia(constraints).then(localStream => {
        console.log("getUserMedia() success, stream created, ready for recording ...");

        // make stream and recorder global
        stream = localStream;
        rec = new MediaRecorder(stream, { mimeType: "audio/flac" });

        // start recording
        rec.start();
        
        rec.ondataavailable = (e) => {
            // Push the recorded media data to the chunks array
            chunks.push(e.data);
        };

        rec.onstop = () => {
            const blob = new Blob(chunks, {type: "audio/flac"});
            chunks = [];
            postAudio(blob);
        }

        console.log("Recording started");

    }).catch(function(err) {
        //enable the record button if getUserMedia() fails
          $("#recordButton").attr("src", "/static/images/rec_button.png")
                            .attr("alt", "rec button");
    });
});


function stopRecording() {
    console.log("recordButton released");

    // Remove the mouseup event handler
    $(window).off();

    // enable the record button back to signal the user we got the message
    $("#recordButton").attr("src", "/static/images/rec_button.png")
                      .attr("alt", "rec button");

    //tell the recorder to stop the recording
    //rec.stop();

    //stop microphone access
    // this triggers the onstop callback
    stream.getAudioTracks().forEach(track => track.stop());

    console.log("end of stopRecording");
}


function postAudio(blob) {

    // Do this with jQuery?
    var xhr=new XMLHttpRequest();
    xhr.onload=function(e) {
        if(this.readyState === 4) {
            //FIXME
            var response = e.target.responseText;
            console.log("Server returned: ", response);
            if (response != "") {
              window.location.href = response;
            } else {
                alert("Voice command not understood.");
            }
        }
    };

    // This uses a form with POST
    var fd = new FormData();
    fd.append("audio_data", blob);
    xhr.open("POST","/submit_audio", true);
    xhr.send(fd);

}
