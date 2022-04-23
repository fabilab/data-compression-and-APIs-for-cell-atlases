(function() {

// from https://github.com/arkochatterjee/flask-audio-recorder/blob/main/static/js/app.js
var gumStream;                      //stream from getUserMedia()
var rec;                            //Recorder.js object
var input;                          //MediaStreamAudioSourceNode we'll be recording

// shim for AudioContext when it's not avb. 
var AudioContext = window.AudioContext || window.webkitAudioContext;
var audioContext; //audio context to help us record

//add events to the button
$("#recordButton").mousedown(function() {
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

    /*
        We're using the standard promise based getUserMedia() 
        https://developer.mozilla.org/en-US/docs/Web/API/MediaDevices/getUserMedia
    */

    navigator.mediaDevices.getUserMedia({audio: true, video:false}).then(function(stream) {
        console.log("getUserMedia() success, stream created, initializing Recorder.js ...");

        /*
            create an audio context after getUserMedia is called
            sampleRate might change after getUserMedia is called, like it does on macOS when recording through AirPods
            the sampleRate defaults to the one set in your OS for your playback device
        */
        audioContext = new AudioContext();

        //update the format 
        //document.getElementById("formats").innerHTML="Format: 1 channel pcm @ "+audioContext.sampleRate/1000+"kHz"

        /*  assign to gumStream for later use  */
        gumStream = stream;

        /* use the stream */
        input = audioContext.createMediaStreamSource(stream);

        /* 
            Create the Recorder object and configure to record mono sound (1 channel)
            Recording 2 channels  will double the file size
        */
        rec = new Recorder(input,{numChannels:1});

        //start the recording process
        rec.record();

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
    rec.stop();

    //stop microphone access
    gumStream.getAudioTracks()[0].stop();

//    //create the wav blob and pass it on to createDownloadLink
//    rec.exportWAV(postAudio);
}


function postAudio(blob) {

    // Do this with jQuery?
    var xhr=new XMLHttpRequest();
    xhr.onload=function(e) {
        if(this.readyState === 4) {
            console.log("Server returned: ", e.target.responseText);
        }
    };

    // This uses a form with POST, I guess one could make an API call
    // or does that become insecure?
    var fd=new FormData();
    fd.append("audio_data", blob, filename);
    xhr.open("POST","/voice_control",true);
    xhr.send(fd);

}

}());
