% Copyright (c) 2012 Romain Hennequin (http://www.romain-hennequin.fr/).
% script to synthesize audio data from a MIDI file.
clear all

% directory of MIDI files
inputDirectory = 'sounds\separation\score informed\Midi';
directory = dir(inputDirectory);

for p = 3:length(directory)
    
    % name of the MIDI file
    MIDIfileName = directory(p).name;
    
    % name of the WAV file
    dots = strfind(MIDIfileName,'.');
    WAVfileName = MIDIfileName(1:dots(end)-1);
    
    % directory of output
    outputDirectory = 'sounds\separation\score informed\separated tracks';
    
    % sampling frequency of sthe synthesized sound.
    fs = 11025;
    
    % instrument of each track (for a list of the possible instrument see help of synth.m)
    
    
    load instrumentNames.mat;
    %instrTrack = {'none','Piano','Bass','Trumpet.vib','none','BbClar','none','none','none','none','none','none','none','none','none','none','none','none'};
    %instrTrack = {'none',instrumentNames{34},instrumentNames{33},instrumentNames{8},instrumentNames{13}};
    instrTrack = {'none','clarinet','flute','violin','trumpet'};
        
    
    % Separate Track Synthesis
    
    % synthesis of audio data
    y = midi2audioTrackByTrack(MIDIfileName,fs,'multi',instrTrack);
    
    
    ymix = sum(cell2mat(y));
    normalizationFactor = max(abs(ymix));
    
    activeTrack = cell(length(y),1);
    
    % export
    for k = 1:length(y)
        if length(y{k})>0
            mkdir(outputDirectory,WAVfileName);
            wavwrite(y{k}/normalizationFactor,fs,16,[outputDirectory '\' WAVfileName '\' WAVfileName '.track' int2str(k) '.' instrTrack{k} '.wav']);
            activeTrack{k} = 1;
        else
            disp(['Track ' int2str(k) ' is empty']);
            activeTrack{k} = 0;
        end
    end
    
    save([outputDirectory '\activeTrack.' WAVfileName '.mat'],'activeTrack');
    
end
