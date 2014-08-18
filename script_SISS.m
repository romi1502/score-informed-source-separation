% Test script for Score informed source separation as proposed in
% R. Hennequin, B. David, R. Badeau. "Score informed audio source separation 
% using a parametric model of non-negative spectrogram", ICASSP11, May 2011
%
% Copyright (c) 2010 Romain Hennequin (http://www.romain-hennequin.fr/).


close all
clear all

%% name of the files

midiFileName = 'roundMidnight.mid';
wavFileName = 'roundMidnightOriginalMix.wav';


%% parameters

% first sample of the soundw
firstSample = 1;

% last sample of the sounds
lastSample = 11025*8;

% decimation factor
decimateFactor = 1;

% size of fft in the spectrogram
Nfft = 1024;

% overlap factor (in [0,1[)
overlapFactor = 0.75;

% beta-divergence used
beta = 1;

% number of templates (atoms)
R = 80;

% pitch space (pitch difference between to successive atoms in semitones)
pitchSpace = 1;

% shift with reference pitch
atomShift = 1;

% reference (lower) MIDI note (atoms correspond to MIDI notes : noteMIDIref, noteMIDIref+1 ... noteMIDIref+R-1)
noteMIDIref = 21;

% number of iteration
Niter = 15;

% number of non harmonic templates
Rnh = 0;

% Constraints
mu.muCorrel12 = 0;
mu.muCorrel19 = 0;
mu.muCorrel24 = 0;
mu.muSmoothH = 0;
mu.muSmoothW = 0;
mu.muSparse = 0;


%% data preparation

% Separate track synthesis
Notes = midiInfo(readmidi(midiFileName), 0);

[x,fs] = wavread(wavFileName);

% transforming x in mono signal
x = toMono(x);

% selecting the sound to analyze :
x = x(firstSample:min(end,lastSample));

% computation of the spectrogram
sp = stft(x,Nfft,hanning(Nfft,'periodic'),Nfft*(1-overlapFactor));
V = abs(sp).^2;

% size of spectrogram
M = size(V,1);
N = size(V,2);

% reference pitch (normalized frequency)
f0ref = 440/(fs/decimateFactor)*2^((noteMIDIref-69)/12);

% central pitch of each band initialization
f0center = ones(1,R);
for r=1:R
    f0center(r) = f0ref*2^((r-atomShift)*pitchSpace/12);
end


%% initialization


tracks = unique(Notes(:,1));

% associate an onset template to the instrument (0->no attack, 1->attack (onset))
associatedTemplateAttack = zeros(size(tracks));
initVal.associatedTemplateAttack = associatedTemplateAttack;

% associate a "stationary noise" (in the sense that its spectral profile is stationary) template to the instrument (0->stationary noise, 1->no stationary noise)
associatedTemplateNoise = zeros(size(tracks));
initVal.associatedTemplateNoise = associatedTemplateNoise;

% time-varying fundamental frequency for each track  (0-> fixed frequency over time, 1->varying fundamental frequency)
initVal.timeVaryingF0 = ones(size(tracks));

% inharmonicity for each track  (0-> no inharmonicity for this instrument, 1->inharmonicity)
inharmonicityEnabled = zeros(size(tracks));
initVal.inharmonicityEnabled = inharmonicityEnabled;

% select tracks
nbTracks = length(tracks);
track = cell(nbTracks,1);
midiTrack = cell(nbTracks,1);
pianoRollTrack = cell(nbTracks,1);
t = cell(nbTracks,1);
nn = cell(nbTracks,1);

% partial amplitudes initialization
initVal.ak = cell(nbTracks,1);

% fundamental frequency initialization
initVal.f0 = cell(nbTracks,1);


% Non harmonic templates initialization
Rnh = sum(associatedTemplateNoise);
Ronset = sum(associatedTemplateAttack);
indexInstrumentOnset = 1;
indexInstrumentNoise = 1;
onsetDuration = 50;% onset duration in ms
Tonset = round(onsetDuration/(Nfft*(1 - overlapFactor)/11025*1000));
initVal.W = rand(M,Rnh);
initVal.H = rand(Rnh,N);
initVal.Wonset = rand(Tonset,M,Ronset);
initVal.Honset = zeros(Ronset,N);
initVal.h = cell(nbTracks,1);

for k = 1:nbTracks
    
    % Build Piano roll from MIDI notes
    track{k} = Notes(Notes(:,1) == tracks(k),:);
    [pianoRollTrack{k},t{k},nn{k}] = piano_roll(track{k},0,length(x)/fs*decimateFactor/(N+1));
    
    pianoRollTrack{k} = pianoRollTrack{k}(:,1:min(N,end));
    t{k} = t{k}(1:min(N,end));
    
    % Remove notes outside the defined range
    if sum(nn{k}>noteMIDIref+R-1)
        warning(['notes outside (above) defined range removed ' int2str(nn{k}(end))])
        pianoRollTrack{k} = pianoRollTrack{k}(1:find(nn{k}==noteMIDIref+R-1,1),:);
        nn{k} = nn{k}(1:find(nn{k}==noteMIDIref+R-1,1));
    end
    if sum(nn{k}<noteMIDIref)
        warning(['notes outside (below) defined range suppressed ' int2str(nn{k}(1))])
        pianoRollTrack{k} = pianoRollTrack{k}(find(nn{k}==noteMIDIref,1):end,:);
        nn{k} = nn{k}(find(nn{k}==noteMIDIref,1):end);
    end
    
    % activation of harmonic templates initialization
    initVal.h{k} = [zeros(nn{k}(1)-noteMIDIref,N); pianoRollTrack{k}; zeros(noteMIDIref+R-1-nn{k}(end),N)];
    initVal.nn{k} = nn{k} - noteMIDIref + 1;
    initVal.h{k} = (smooth(initVal.h{k}',5)'>0)*1.00;
    
    
    
    % associated non harmonic template
    if associatedTemplateAttack(k)
        % The template only takes non null value around the attack time
        % (onset template)
        initVal.Honset(indexInstrumentOnset,:) = smooth([(sum(max(diff(initVal.h{k},1,2),0)))~=0,0],4);
        indexInstrumentOnset = indexInstrumentOnset + 1;
    end
    if associatedTemplateNoise(k)
        initVal.H(indexInstrumentNoise,:) = sum(initVal.h{k});
        indexInstrumentNoise = indexInstrumentNoise + 1;
    end
    
    % pitch initialization
    if initVal.timeVaryingF0(k)
        initVal.f0{k} = ones(N,R);
        for r = initVal.nn{k}
            initVal.f0{k}(:,r) = f0center(r);
        end
    else
        for r = initVal.nn{k}
            initVal.f0{k}(r) = f0center(r);
        end
    end
    
    initVal.ak{k} =  (1./(1:M)').^2*ones(1,R);%ones(M,R);
    
    
    % inharmonicity initialization
    if inharmonicityEnabled(k);
        
        inharmoSup = 3*10^-4*ones(100,1);
        inharmoSup(1:18) = 3*10^-4*10.^((17:-1:0)/22);
        inharmoSup(30:100) = 3*10^-4*10.^((0:70)/22);
        inharmoSup = inharmoSup(noteMIDIref-21+1:noteMIDIref-21+R);
        
        inharmoInit = inharmoSup/3;
        inharmoInit = inharmoInit(noteMIDIref-21+1:noteMIDIref-21+R);
        
        initVal.Binharm{k} = inharmoInit;
    else
        initVal.Binharm{k} = zeros(R,1);
    end
    
end



%% Computation of the decomposition

[f0, ak, w, h, W, H, Wonset, Honset, Binharm, Lambda] = variablePitchNMFMultiInstrument(V,R,f0center,Niter,mu,Rnh,Ronset,beta,initVal);


%% Separation


% compute the harmonic associated to each track
LambdaTracks = cell(nbTracks,1);
Lambda = zeros(size(V));
for k=1:nbTracks
    LambdaTracks{k} = zeros(size(V));
    
    if initVal.timeVaryingF0(k)
        for t=1:N
            for r = 1:R;
                if initVal.h{k}(r,t)
                    LambdaTracks{k}(:,t) = LambdaTracks{k}(:,t) + h{k}(r,t)*w{k,r}(:,t);
                end
            end
        end
    else
        LambdaTracks{k} = LambdaTracks{k} + cell2mat(w(k,:))*h{k}(initVal.nn{k},:);
    end
    Lambda = Lambda + LambdaTracks{k};
end

LambdaOnset = zeros(M,N);
for k = 1:Ronset
    for t=0:Tonset-1
        LambdaOnset = LambdaOnset + reshape(Wonset(t+1,:,k),M,1)* shiftLR(Honset(k,:),-t);
    end
end


Lambda = Lambda + W*H + LambdaOnset;


% sound synthesis of harmonic part
xr = cell(nbTracks,1);
for k = 1:nbTracks
    xr{k} = istft(sp.*LambdaTracks{k}(:,:)./(Lambda + eps),Nfft,hanning(Nfft,'periodic')',Nfft*(1-overlapFactor));
end


% sound synthesis of "stationary" noise part
indexNoise = find(associatedTemplateNoise);
xNoise = cell(Rnh,1);
for k = 1:Rnh
    xNoise{indexNoise(k)} = istft(sp.*(W(:,k)*H(k,:))./(Lambda + eps),Nfft,hanning(Nfft,'periodic')',Nfft*(1-overlapFactor));
end

% sound synthesis of onset part
indexOnset = find(associatedTemplateAttack);
xOnset = cell(Ronset,1);
for k = 1:Ronset
    LambdaOnsetk{k} = zeros(M,N);
    for t=0:Tonset-1
        LambdaOnsetk{k} = LambdaOnsetk{k} + reshape(Wonset(t+1,:,k),M,1)* shiftLR(Honset(k,:),-t);
    end
    
    xOnset{indexOnset(k)} = istft(sp.*(LambdaOnsetk{k})./(Lambda + eps),Nfft,hanning(Nfft,'periodic')',Nfft*(1-overlapFactor));
end

xSep = cell(nbTracks,1);
for k = 1:nbTracks
    xSep{k} = xr{k};
    
    if associatedTemplateNoise(k)
        xSep{k} = xSep{k} + xNoise{k};
    end
    if associatedTemplateAttack(k)
        xSep{k} = xSep{k} + xOnset{k};
    end
end


