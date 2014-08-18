function [y synthOK] = synth(freq,dur,amp,Fs,type)
% [y synthOK] = synth(freq,dur,amp,Fs,type)
%
% Synthesize a single note
%
% Inputs:
%  freq - frequency in Hz
%  dur - duration in seconds
%  amp - Amplitude in range [0,1]
%  Fs -  sampling frequency in Hz
%  type - string to select synthesis type
%         current options: 'fm', 'sine', 'saw', 'Oboe', 'EbClar', 'BbClar',
%         'altosax.novib', 'altosax.vib', 'trumpet.novib',
%         'trumpet.vib', 'tuba', 'violin.arco', 'basstrombone', 'piano',
%         'bass', 'electricpiano', 'classicguitar', 'organfull',
%         'organfullvib', 'organcheap'

% Copyright (c) 2009 Ken Schutte
% more info at: http://www.kenschutte.com/midi

synthOK = 0;

if nargin<5
    error('Five arguments required for synth()');
end

load instrumentNames.mat;

N = floor(dur*Fs);

if N == 0
    warning('Note with zero duration.');
    y = [];
    return;
    
elseif N < 0
    warning('Note with negative duration. Skipping.');
    y = [];
    return;
end

n=0:N-1;
if (strcmp(type,'sine'))
    y = amp.*sin(2*pi*n*freq/Fs);
    synthOK = 1;
elseif (strcmp(type,'saw'))
    
    T = (1/freq)*Fs;     % period in fractional samples
    ramp = (0:(N-1))/T;
    y = ramp-fix(ramp);
    y = amp.*y;
    y = y - mean(y);
    synthOK = 1;
    
elseif (strcmp(type,'fm'))
    
    t = 0:(1/Fs):dur;
    envel = interp1([0 dur/6 dur/3 dur/5 dur], [0 1 .75 .6 0], 0:(1/Fs):dur);
    I_env = 5.*envel;
    y = envel.*sin(2.*pi.*freq.*t + I_env.*sin(2.*pi.*freq.*t));
    synthOK = 1;
    
elseif sum(ismember(instrumentNames,type))

    velocity = amp * 127;
    if velocity>90
        dynamic = 'ff';
    elseif velocity>41
        dynamic = 'mf';
    else
        dynamic = 'pp';
    end
    
    midinum = 69 + 12*log2(freq/440);
    t = 0:(1/Fs):dur;
    
    fileName = ['sounds\instruments\BDD perso\son isolés\' type '\' type '.' dynamic '.' int2str(midinum) '.wav'];
    if exist(fileName)
        [x,fsSample] = wavread(fileName);
        x = resample(x,Fs,fsSample);
        x = x./max(abs(x));
        synthOK = 1;
    else
        dynamic
        x = zeros(length(t),1);
    end
    
    if length(x)<length(t)
        x = [x;zeros(length(t)-length(x)+1,1)];
    end
    y = amp*x(1:length(t))';
    %y = x(1:length(t))';
    
else
    disp('Unknown synthesis type');
    t = 0:(1/Fs):dur;
    y = zeros(1,length(t));
end

    
    
%     (strcmp(type,'oboe'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['oboe.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'ebclar'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['EbClar.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'bbclar'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['BbClar.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'altosax.novib'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['AltoSax.NoVib.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'altosax.vib'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['AltoSax.Vib.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'altosax.novib'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['AltoSax.NoVib.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'trumpet.novib'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['Trumpet.novib.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'trumpet.vib'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['Trumpet.vib.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'tuba'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['Tuba.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'violin.arco'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['Violin.arco.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'basstrombone'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['BassTrombone.ff.' int2str(midinum) '.wav']);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'piano'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['piano.ff.' int2str(midinum) '.wav'],44100*5);
%         x = toMono(x);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
%     
%     
% elseif (strcmp(type,'bass'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['Bass.' int2str(midinum) '.wav']);
%         x = toMono(x);
%         x = resample(x,Fs,fsSample);
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
%        
% elseif (strcmp(type,'electricpiano'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['electricPiano.ff.' int2str(midinum) '.wav']);
%         x = toMono(x);
%         x = resample(x,Fs,fsSample)/4;
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'organfull'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['organFull.ff.' int2str(midinum) '.wav']);
%         x = toMono(x);
%         x = resample(x,Fs,fsSample)/12;
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'organfullvib'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['organFullVib.ff.' int2str(midinum) '.wav']);
%         x = toMono(x);
%         x = resample(x,Fs,fsSample)/6;
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'organcheap'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['organCheap.ff.' int2str(midinum) '.wav']);
%         x = toMono(x);
%         x = resample(x,Fs,fsSample)/12;
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% elseif (strcmp(type,'classicguitar'))
%     
%     midinum = 69 + 12*log2(freq/440);
%     t = 0:(1/Fs):dur;
%     try
%         [x,fsSample] = wavread(['classicGuitar.ff.' int2str(midinum) '.wav']);
%         x = toMono(x);
%         x = resample(x,Fs,fsSample)/4;
%     catch ME
%         if strcmp(ME.identifier,'wavread:InvalidFile')
%             x = zeros(length(t),1);
%         else
%             throw(ME);
%         end
%     end
%     if length(x)<length(t)
%         x = [x;zeros(length(t)-length(x)+1,1)];
%     end
%     y = amp*x(1:length(t))';
% 
% else
%     %disp('Unknown synthesis type');
%     t = 0:(1/Fs):dur;
%     y = zeros(1,length(t));
% end

% smooth edges w/ 10ms ramp
if (dur > .1)
    L = 2*fix(.05*Fs)+1;  % L odd
    ramp = bartlett(L)';  % odd length
    L = ceil(L/2);
    y(1:L) = y(1:L) .* ramp(1:L);
    y(end-L+1:end) = y(end-L+1:end) .* ramp(end-L+1:end);
end
