function play_Callback(hObject, eventdata, handles)
fs = handles.fs*(1 + handles.FSQ);
sound(handles.x, fs);


% --- Executes on button press in formant_freq.
function formant_freq_Callback(hObject, eventdata, handles)
    h = spectrum.welch;
    hs = psd(h,handles.x,'fs',handles.fs);
    figure;
plot(hs);

function exit_Callback(hObject, eventdata, handles)
cl = questdlg('Do you want to EXIT?','EXIT',...
            'Yes','No','No');
switch cl
    case 'Yes'
        close();
        clear all;
        return;
    case 'No'
        quit cancel;
end 

function record_Callback(hObject, eventdata, handles)
    fs = 44100;
    y = wavrecord(88200,fs);
    [filename, pathname] = uiputfile('*.wav', 'Pick an M-file');
    cd (pathname);
    wavwrite(y,fs,filename);
    sound(y,fs);
    handles.x = y;
    handles.fs = fs;
    axes(handles.axes1);
    time = 0:1/fs:(length(handles.x)-1)/fs;
    plot(time,handles.x);
    title('Original Signal');
    axes(handles.axes2);
    specgram(handles.x, 1024, handles.fs);
    title('Spectrogram of Original Signal');
guidata(hObject, handles);

function load_file_Callback(hObject, eventdata, handles)
    clc;
    [FileName,PathName] = uigetfile({'*.wav'},'Load Wav File');
    [x,fs] = wavread([PathName '/' FileName]);
    handles.x = x;
    handles.fs = fs;
    axes(handles.axes1);
    time = 0:1/fs:(length(handles.x)-1)/fs;
    plot(time,handles.x);
    title('Original Signal');
    axes(handles.axes2);
    specgram(handles.x(:,1), 1024, handles.fs);
    title('Spectrogram of Original Signal');
guidata(hObject,handles);

function load_random_Callback(hObject, eventdata, handles)
    clc;
    fs = 8200;
    x = randn(5*fs,1);
	handles.x = x;
    handles.fs = fs;
    axes(handles.axes1);
    time = 0:1/fs:(length(handles.x)-1)/fs;
    plot(time,handles.x);
    title('Original Signal');
    axes(handles.axes2);
    specgram(handles.x, 1024, handles.fs);
    title('Spectrogram of Original Signal');
guidata(hObject, handles);

function slider2_Callback(hObject, eventdata, handles)
    handles.FSQ = (get(hObject,'Value'));
    set(handles.edit1, 'String', [sprintf('%.1f',handles.FSQ) ''] );
guidata(hObject, handles);
