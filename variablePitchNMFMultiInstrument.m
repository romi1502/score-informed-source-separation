function [f0, ak, w, h, W, H, Wonset, Honset, Binharm, Lambda] = variablePitchNMFMultiInstrument(V,R,f0center,Niter,mu,Rnh,Ronset,beta,initVal)
% [f0, ak, w, h, W, H, Wonset, Honset, Binharm, Lambda] = variablePitchNMFMultiInstrument(V,R,f0center,Niter,mu,Rnh,Ronset,beta,initVal)
%       - V : power spectrogram to decompose
%       - R : number of templates (atoms)
%       - f0center : center frequency of each band
%       - Niter : number of iteration
%       - mu : various constraints (default: all = 0).
%              a struct with fields : - muSparse: sparsness constraints on the columns of h
%                                     - muCorrel12,muCorrel19 and muCorrel24: uncorrelation constraints between activations of octaves/twelves/double-octave
%                                     - muSmoothH: smoothness onstraints on the rows of h
%                                     - muSmoothW: smoothness onstraints on the successive amplitudes of harmonics
%       - Rnh : number of non harmonic templates (classic NMF templates), optionnal : default = 0
%       - beta : beta used for beta-divergence, optionnal : default = 1
%       - h0 : time/frequency mask  for each instruments (h0 is cell-array. h0{k} is the T/F mask associated to the kth instrument)
%       - initVal : initial values (struct with field f0 ak h W H Binharm),
%       optionnal, default = random initialization
%
% Copyright (c) 2010 Romain Hennequin (http://www.romain-hennequin.fr/).


% plot at each iteration?
plotResults = true;

% free harmonic amplitude, fixed template
freeHarmonicTemplate = true;

% same harmonic distribution for every atom
sameHarmonicDistribution = true;

% beta-divergence used
if nargin<7
    beta = 1;
end
if nargin<6
    Rnh = 0;
end
if nargin<5
    mu.muSparse = 0;
    mu.muCorrel12 = 0;
    mu.muCorrel19 = 0;
    mu.muCorrel24 = 0;
    mu.muSmoothH = 0;
    mu.muSmoothW = 0;
end

% Minimum value for h to avoid numerical instabilities
hThreshold = 10^-10;

% size of input spectrogram
M = size(V,1);
N = size(V,2);

%% Window computation

% window length
T = 2*(M-1);

% window type ('hanning' or 'hamming')
windowType = 'hamming';

switch lower(windowType)
    case {'hanning','hann'}
        a = 0.5;
        b = 0.5;
    case 'hamming'
        a = 0.54;
        b = 0.46;
    otherwise
        disp('Unknown window, replaced by Hann window')
        a = 0.5;
        b = 0.5;
end

% Quantized frequencies
division = T*4000;
f = (-1/4:1/division:1/4)';

% Quantized pulses
w = 2*pi*f;

% Squared analysis window
g = (2-2*cos(T*w)).*(T^2*w.^2*(b-a)+4*pi^2*a).^2./(T^2*w.^3 - 4*pi^2*w).^2;
% Singular value (limits)
g(abs(f)<10^-10) = a^2*T^2;
g(abs(f-1/T)<10^-10) = b^2*T^2/4;
g(abs(f+1/T)<10^-10) = b^2*T^2/4;

% Squared analysis window derivative
derivat = 1/(4*pi^2)*( (2-2*cos(2*pi*T*f)).*2.*(f.^2*T^2*(b-a)+a)./(f.^2.*(f.^2*T^2-1).^2)   .*...
    ( 2*f*T^2*(b-a) - (2*f*T^2.*(f.^2*T^2*(b-a)+a))./(f.^2*T^2 - 1) -  (f.^2*T^2*(b-a)+a)./f) + ...
    4*pi*T*sin(2*pi*T*f) .* (T^2*f.^2*(b-a) + a).^2./(f.^2.*(T^2*f.^2-1).^2) );

% Squared analysis window derivative divided by f
dg = derivat./f;
% Singular value (limits)
[mini imini] = min(abs(f));
parabol = polyfit([imini-5 imini+5 imini+4]',dg([imini-5 imini+5 imini+4]),2);
dg(imini-2:imini+2) = polyval(parabol,imini-2:imini+2);
[mini imini] = min(abs(f-1/T));
dg(imini) = dg(imini-1)/2+dg(imini+1)/2;
[mini imini] = min(abs(f+1/T));
dg(imini) = dg(imini-1)/2+dg(imini+1)/2;

len = length(g);



%% initialization of algorithm variables

% normalized frequency index
freq = (0:1/(M-1):1)'/2;


%% initialization of the parameters of the decomposition

if nargin<8 % no initial values : random initialization
    
    
    error('number of arguments lower than 8 not supported yet')
    
    % nb instrument
    nbInstr = 1;
    
    % pitch initialization
    f0{1} = ones(N,R);
    for r=1:R
        f0{1}(:,r) = f0center(r);
    end
    
    % Non harmonic templates initialization
    W = ones(M,Rnh);
    H = rand(Rnh,N);
    
    % partial amplitudes initialization
    ak{1} = ones(M,R);
        
    % activation of harmonic templates initialization
    h{1} = rand(R,N);
    h{1}(h{1}<hThreshold) = 10*hThreshold;
        
    if inharmonicityEnabled
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TO BE CHANGED fs NOT FIXED!!!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fs = 11025;
        refMIDInote = 12*log2(f0center(1)*fs/440) + 69;
        inharmoSup = 3*10^-4*ones(100,1);
        inharmoSup(30:88) = 3*10^-4*10.^((0:58)/22);
        
        inharmoSup = inharmoSup(refMIDInote-27+1:refMIDInote-27+R);
        
        inharmoInit = 0.5*10^-4*ones(88,1);
        inharmoInit(30:88) = 0.5*10^-4*10.^((0:58)/22);
        
        inharmoInit = inharmoInit(refMIDInote-27+1:refMIDInote-27+R);
        
        Binharm{1} = inharmoInit;
        
    else
        Binharm{1} = zeros(R,1);
    end
    
else % take initial values
    
    % pitch initialization
    f0 = initVal.f0;
    
    % Non harmonic templates initialization
    W = initVal.W;
    H = initVal.H;
    
    % Onset templates initialization
    Wonset = initVal.Wonset;
    Honset = initVal.Honset;
    Ronset = size(Wonset,3);
    Tonset = size(Wonset,1);
    One = ones(M,N);
    Ht = zeros(Tonset,Ronset,N);
    LambdaOnset = zeros(M,N);
    
    % computation of LambdaOnset
    for t=0:Tonset-1
        LambdaOnset = LambdaOnset + reshape(Wonset(t+1,:,:),M,Ronset)* shiftLR(Honset,-t);
    end
    
    % partial amplitudes initialization
    ak = initVal.ak;
    
    % activation of harmonic templates initialization
    h = initVal.h;
    
    nbInstr = length(h);
    
    nn = initVal.nn;
    
    % inharmonicity initialization
    inharmoSup = cell(nbInstr,1);
    inharmonicityEnabled = zeros(nbInstr,1);
    Binharm = cell(nbInstr,1);
    for instr = 1:nbInstr
        inharmoSup{instr} = 6*initVal.Binharm{instr};
        inharmonicityEnabled(instr) = initVal.inharmonicityEnabled(instr);
        if inharmonicityEnabled(instr)
            Binharm{instr} = initVal.Binharm{instr};
        else
            Binharm{instr} = zeros(R,1);
        end
    end
    
    timeVaryingF0 = initVal.timeVaryingF0;
end

% minimum value of harmonic templates
wThreshold = 10^-10;

% width of peak
Deltamr = min(ceil(f0center*T*2),T/8);

% harmonic templates and estimated spectrogram initialization
w = cell(nbInstr,R);

for instr = 1:nbInstr
    if timeVaryingF0(instr)
        for r = nn{instr}
            if sum(h{instr}(r,:))>hThreshold*sum(h{instr}(r,:)>0)
                w{instr,r} = zeros(M,N);
            end
        end
    else
        for r = nn{instr}
            w{instr,r} = zeros(M,1);
        end
    end
end


% normalization of ak
for instr = 1:nbInstr
    for r = nn{instr}
        sumak = sum(ak{instr}(1:floor(0.5/f0center(r)),r));
        ak{instr}(:,r) = ak{instr}(:,r)/sumak;
    end
end


% Estimate of V
Lambda = zeros(M,N);
for instr = 1:nbInstr
    if timeVaryingF0(instr) %time varying fundamental frequency
        
        for r = nn{instr}
            for t = 1:N
                if h{instr}(r,t)>0
                    kf0 = f0{instr}(t,r)*2*M;
                    for k = 1:floor(0.5/f0center(r))
                        deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                        fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                        w{instr,r}(deb:fin,t) = w{instr,r}(deb:fin,t) + ak{instr}(k,r)* g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(t,r)) + len/2));
                    end
                    w{instr,r}(:,t) = w{instr,r}(:,t) + wThreshold;
                    Lambda(:,t) = Lambda(:,t)+ w{instr,r}(:,t)*h{instr}(r,t);
                end
            end
        end
        
    else %non time varying fundamental frequency
        
        for r = nn{instr}
            if sum(h{instr}(r,:))>0
                kf0 = f0{instr}(r)*2*M;
                for k = 1:floor(0.5/f0center(r))
                    deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                    fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                    w{instr,r}(deb:fin) = w{instr,r}(deb:fin) + ak{instr}(k,r)* g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(r)) + len/2));
                end
                w{instr,r} = w{instr,r} + wThreshold;           
            end
        end
        Lambda = Lambda + cell2mat(w(instr,:))*h{instr}(nn{instr},:);
        
    end
end

Lambda = Lambda + W*H + LambdaOnset;


if R>1
    pitchSpace = 12*(log(f0center(2)) - log(f0center(1)))/log(2);
else
    pitchSpace = 1;
end



%% iterations

for iter = 1:Niter
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update of fundamental frequency (f0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for instr = 1:nbInstr
        
        if timeVaryingF0(instr) % time-varying fundamental frequency
            
            for r = nn{instr}
                for t = 1:N
                    if h{instr}(r,t)>hThreshold
                        
                        f0current = f0{instr}(t,r);
                        G = 0;
                        F = 0;
                        Lambdabeta2 = Lambda(:,t).^(beta-2);
                        LambdaFreq = Lambda(:,t).*freq;
                        VFreq = V(:,t).*freq;
                        kf0 = f0{instr}(t,r)*2*M;
                        
                        for k = 1:floor(0.5/f0center(r))
                            kxf0current = k*f0current;
                            deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                            fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                            
                            vloc = Lambdabeta2(deb:fin).*(-dg(round(division*(freq(deb:fin)-sqrt(1+k^2*Binharm{instr}(r))*kxf0current) + len/2)));
                            
                            F = F + sum((LambdaFreq(deb:fin)+(kxf0current)*V(deb:fin,t)).*vloc)*(ak{instr}(k,r)*k)*h{instr}(r,t);
                            G = G + sum((VFreq(deb:fin)+(kxf0current)*Lambda(deb:fin,t)).*vloc)*(ak{instr}(k,r)*k)*h{instr}(r,t);
                        end
                        
                        FSmooth = 4/f0center(r)^2*f0{instr}(t,r)^2;
                        if t>2&&t<N
                            GSmooth = 2/f0center(r)^2*f0{instr}(t,r)*(f0{instr}(t+1,r)+f0{instr}(t-1,r));
                        elseif t>2
                            GSmooth = 2/f0center(r)^2*f0{instr}(t,r)*f0{instr}(t-1,r);
                        else
                            GSmooth = 2/f0center(r)^2*f0{instr}(t,r)*f0{instr}(t+1,r);
                        end
                        muSmoothF0 = 0;
                        
                        f0{instr}(t,r) = f0{instr}(t,r)* (G + muSmoothF0*GSmooth)/ (F + muSmoothF0*FSmooth + eps);
                        
                        % ensure that fundamental frequencies remain in the right band
                        if (f0{instr}(t,r)/f0center(r))>2^(pitchSpace/(2*12))||(f0{instr}(t,r)/f0center(r))<2^(-pitchSpace/(2*12))
                            f0{instr}(t,r) = f0center(r)*2^((-pitchSpace+rand(1))/(2*12));                            
                        end
                        if isnan(f0{instr}(t,r))
                            f0{instr}(t,r) = f0center(r);
                            h{instr}(r,t) = 0;
                        end
                    end
                end
            end
            
        else % non time-varying fundamental frequency
            for r = nn{instr}
                G = 0;
                F = 0;
                f0current = f0{instr}(r);
                kf0 = f0{instr}(r)*2*M;
                
                for t = 1:N
                    if h{instr}(r,t)>hThreshold
                        
                        Lambdabeta2 = Lambda(:,t).^(beta-2);
                        LambdaFreq = Lambda(:,t).*freq;
                        VFreq = V(:,t).*freq;
                        
                        for k = 1:floor(0.5/f0center(r))
                            kxf0current = k*f0current;
                            deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                            fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                            
                            vloc = Lambdabeta2(deb:fin).*(-dg(round(division*(freq(deb:fin)-sqrt(1+k^2*Binharm{instr}(r))*kxf0current) + len/2)));
                     
                            F = F + sum((LambdaFreq(deb:fin)+(kxf0current)*V(deb:fin,t)).*vloc)*(ak{instr}(k,r)*k);
                            G = G + sum((VFreq(deb:fin)+(kxf0current)*Lambda(deb:fin,t)).*vloc)*(ak{instr}(k,r)*k);
                        end                                             
                    end
                end
                
                if F>0
                    f0{instr}(r) = f0{instr}(r)*G / (F+eps);
                end
                
                % ensure that fundamental frequencies remain in the right band
                if (f0{instr}(r)/f0center(r))>2^(pitchSpace/(2*12))||(f0{instr}(r)/f0center(r))<2^(-pitchSpace/(2*12))
                    f0{instr}(r) = f0center(r)*2^((-pitchSpace+rand(1))/(2*12));
                end
                if isnan(f0{instr}(r))
                    f0{instr}(r) = f0center(r);
                    h{instr}(r,:) = 0;
                    error('f0 is NaN');
                end
            end
        end
        
    end
    
    
    % recomputation of w (templates) and Lambda (estimated power spectrogram)
    Lambda(:) = 0;
    for instr = 1:nbInstr
        
        if timeVaryingF0(instr) % time-varying fundamental frequency
            
            for r = nn{instr}
                w{instr,r}(:) = 0;
                for t = 1:N
                    if h{instr}(r,t)>hThreshold
                        kf0 = f0{instr}(t,r)*2*M;
                        for k=1:floor(0.5/f0center(r))
                            deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                            fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                            w{instr,r}(deb:fin,t) = w{instr,r}(deb:fin,t) + ak{instr}(k,r)* g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(t,r)) + len/2));
                        end
                        w{instr,r}(:,t) = w{instr,r}(:,t) + wThreshold;
                        Lambda(:,t) = Lambda(:,t)+ w{instr,r}(:,t)*h{instr}(r,t);
                    end
                end
                if sum(isnan(w{instr,r}(:)))
                    disp('w is NaN, function returns')
                    return
                end
            end
            
        else % non time-varying fundamental frequency
            
            for r = nn{instr}
                w{instr,r}(:) = 0;
                if sum(h{instr}(r,:))>hThreshold*sum(h{instr}(r,:)>0)
                    kf0 = f0{instr}(r)*2*M;
                    for k=1:floor(0.5/f0center(r))
                        deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                        fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                        w{instr,r}(deb:fin) = w{instr,r}(deb:fin) + ak{instr}(k,r)* g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(r)) + len/2));
                    end
                    w{instr,r} = w{instr,r} + wThreshold;
                end
                
                if sum(isnan(w{instr,r}(:)))
                    disp('w is NaN, function returns')
                    return
                end
            end
            Lambda = Lambda + cell2mat(w(instr,:))*h{instr}(nn{instr},:);
            
        end
    end
    
    Lambda = Lambda + W*H + LambdaOnset;
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update of inharmonicity (Binharm)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for instr = 1:nbInstr
        
        if ~timeVaryingF0(instr) % time-varying fundamental frequency
            
            if inharmonicityEnabled(instr)
                
                % type of update for inharmonicity
                updateInharmoType = 'grid';
                
                switch updateInharmoType
                    case 'grid' % grid update
                        
                        if mod(iter,10)==0
                            for r = nn{instr}
                                wr = zeros(M,1);
                                % Estimate of Vr (spectrogram associated to the rth template)
                                Lambdar = zeros(M,N);
                                
                                if sum(h{instr}(r,:))>hThreshold*sum(h{instr}(r,:)>0)
                                    kf0 = f0{instr}(r)*2*M;
                                    for k=1:floor(0.5/f0center(r))
                                        deb = max(floor(k*kf0*sqrt(1+k^2*Binharm{instr}(r))-Deltamr(r)),1);
                                        fin = min(ceil(k*kf0*sqrt(1+k^2*Binharm{instr}(r))+Deltamr(r)),M);
                                        wr(deb:fin) = wr(deb:fin) + ak{instr}(k)* g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(r)) + len/2));
                                    end
                                    wr = wr + wThreshold;
                                    
                                end
                                
                                Lambdar = wr*h{instr}(r,:);
                                
                                if (Binharm{instr}(r)>0)&&(sum(h{instr}(r,:))>hThreshold*sum(h{instr}(r,:)>0))
                                    inharmTest = 0:min(inharmoSup{instr}(r),(3-iter/Niter*1.5)*Binharm{instr}(r))/100:min(inharmoSup{instr}(r),(3-iter/Niter*1.5)*Binharm{instr}(r));
                                    bdivinharm = zeros(size(inharmTest));
                                    
                                    p = 0;
                                    for inharmonicity = inharmTest
                                        p = p+1;
                                        bdivinharm(p) = betaDivInharm(inharmonicity,f0{instr},h{instr},ak{instr},V,Lambda - Lambdar,beta,g,Deltamr,f0center,freq,len,hThreshold,wThreshold,division,M,N,r);
                                    end
                                    figure(7)
                                    
                                    plot(inharmTest,bdivinharm);
                                    handleLine = line([Binharm{instr}(r) Binharm{instr}(r)],[min(bdivinharm)*0.999 1.001*max(bdivinharm)]);
                                    set(handleLine,'color','r','LineStyle','--');
                                    
                                    [mi,imi] = min(bdivinharm);
                                    Binharm{instr}(r) = inharmTest(imi);
                                    drawnow
                                end
                            end
                        end
                        
                    case 'golden'% golden search update
                        
                        options = optimset('TolFun',1e-10,'MaxIter',25);
                        
                        if mod(iter,5)==0
                            for r=1:R
                                wr = zeros(M,1);
                                % Estimate of Vr (spectrogram associated to the rth template)
                                Lambdar = zeros(M,N);
                                
                                if sum(h{instr}(r,:))>hThreshold*sum(h{instr}(r,:)>0)
                                    kf0 = f0{instr}(r)*2*M;
                                    for k=1:floor(0.5/f0center(r))
                                        deb = max(floor(k*kf0*sqrt(1+k^2*Binharm(r))-Deltamr(r)),1);
                                        fin = min(ceil(k*kf0*sqrt(1+k^2*Binharm(r))+Deltamr(r)),M);
                                        wr(deb:fin) = wr(deb:fin) + ak(k,r)* g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm(r))*f0{instr}(r)) + len/2));
                                    end
                                    wr(:) = wr(:) + wThreshold;
                                end
                                
                                Lambdar = wr*h{instr};
                                
                                [Binharm,bdiv1] = fminbnd(@(inharmonicity) betaDivInharm(inharmonicity,f0{instr},h{instr},ak{instr},V,Lambda - Lambdar,beta,g,Deltamr,f0center,freq,len,hThreshold,wThreshold,division,M,N,r),Binharm(r)/2,Binharm(r)*2,options);
                                
                                Binharm{instr}(r) = Binharm;
                                
                                if Binharm{instr}(r)>inharmoSup{instr}(r)*0.99
                                    Binharm{instr}(r) = 0;
                                    disp('inharmonicity above threshold')
                                end
                            end
                        end
                        
                    case 'mult'
                        
                        % multiplicative update
                        for r = nn{instr}
                            G = 0;
                            F = 0;
                            
                            if Binharm{instr}(r)~=0
                                for t=1:N
                                    
                                    f0current = f0{instr}(t,r);
                                    Lambdabeta2 = Lambda(:,t).^(beta-2);
                                    LambdaFreq = Lambda(:,t).*freq;
                                    VFreq = V(:,t).*freq;
                                    kf0 = f0{instr}(t,r)*2*M;
                                    
                                    for k=1:floor(0.5/f0center(r))
                                        kxf0current = k*f0current;
                                        deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                                        fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                                        
                                        vloc = Lambdabeta2(deb:fin).*(-dg(round(division*(freq(deb:fin)-kxf0current*sqrt(1+k^2*Binharm{instr}(r))) + len/2)));
                                        
                                        F = F + sum((LambdaFreq(deb:fin)+(kxf0current*sqrt(1+k^2*Binharm{instr}(r)))*V(deb:fin,t)).*vloc)*(ak{instr}(k,r)*f0{instr}(t,r)*k^3*h{instr}(r,t)/(sqrt(1+k^2*Binharm{instr}(r))));
                                        G = G + sum((VFreq(deb:fin)+(kxf0current*sqrt(1+k^2*Binharm{instr}(r)))*Lambda(deb:fin,t)).*vloc)*(ak{instr}(k,r)*f0{instr}(t,r)*k^3*h{instr}(r,t)/(sqrt(1+k^2*Binharm{instr}(r))));
                                    end
                                    
                                end
                                if isnan(F)||isnan(G)
                                    disp('NaN update')
                                    return;
                                end
                                Binharm{instr}(r) = Binharm{instr}(r)*(G/(F+eps))^100;
                                
                                % Limit inharmonicity
                                if Binharm{instr}(r)>inharmoSup{instr}(r)%0.005
                                    Binharm{instr}(r) = 0;
                                end
                                
                                
                            end
                        end
                end 
            end            
        end
    end
    
    % recomputation of w (templates) and Lambda (estimated power spectrogram)
    Lambda(:) = 0;
    for instr = 1:nbInstr
        
        if timeVaryingF0(instr) % time-varying fundamental frequency
            
            for r = nn{instr}
                w{instr,r}(:) = 0;
                for t=1:N
                    if h{instr}(r,t)>hThreshold
                        kf0 = f0{instr}(t,r)*2*M;
                        for k=1:floor(0.5/f0center(r))
                            deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                            fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                            w{instr,r}(deb:fin,t) = w{instr,r}(deb:fin,t) + ak{instr}(k,r)* g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(t,r)) + len/2));
                        end
                        w{instr,r}(:,t) = w{instr,r}(:,t) + wThreshold;
                        Lambda(:,t) = Lambda(:,t)+ w{instr,r}(:,t)*h{instr}(r,t);
                    end
                end
                
                if sum(isnan(w{instr,r}(:)))
                    disp('w is NaN, function returns')
                    return
                    
                end
            end
            
        else % non time-varying fundamental frequency
            
            for r = nn{instr}
                w{instr,r}(:) = 0;
                if sum(h{instr}(r,:))>hThreshold*sum(h{instr}(r,:)>0)
                    kf0 = f0{instr}(r)*2*M;
                    for k=1:floor(0.5/f0center(r))
                        deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                        fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                        w{instr,r}(deb:fin) = w{instr,r}(deb:fin) + ak{instr}(k,r)* g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(r)) + len/2));
                    end
                    w{instr,r} = w{instr,r} + wThreshold;
                end
                
                if sum(isnan(w{instr,r}(:)))
                    disp('w is NaN, function returns')
                    return
                    
                end
            end
            Lambda = Lambda + cell2mat(w(instr,:))*h{instr}(nn{instr},:);
            
        end
        
    end
    
    
    Lambda = Lambda + W*H + LambdaOnset;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update of activation (h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    htemp = h;
    muSparse = mu.muSparse;
    muCorrel12 = mu.muCorrel12;
    muCorrel19 = mu.muCorrel19;
    muCorrel24 = mu.muCorrel24;
    muSmoothH = mu.muSmoothH;
    
    
    for instr = 1:nbInstr
        
        if timeVaryingF0(instr) % time-varying fundamental frequency
            
            NormH = sqrt(sum(h{instr}.^2,2));
            dotH = h{instr}*h{instr}';
            for r=nn{instr}
                if sum(h{instr}(r,:))>hThreshold*sum(h{instr}(r,:)>0)
                    squaredNormh = sum(h{instr}(r,:).^2);
                    squaredDiffh = sum(diff(h{instr}(r,:)).^2);
                    for t=1:N
                        if h{instr}(r,t)>hThreshold
                            vloc = w{instr,r}(:,t).*(Lambda(:,t).^(beta-2));
                 
                            % Distance function
                            G = sum(vloc.*V(:,t));
                            F = sum(vloc.*Lambda(:,t));
                            
                           
                            htemp{instr}(r,t) = h{instr}(r,t).* G ./ (F + eps) ;
                        end
                    end
                end
            end
            if sum(isnan(h{instr}(:)))
                disp('h is Nan. Function returns')
                return;
            end
            
        else % non time-varying fundamental frequency
            
            NormH = sqrt(sum(h{instr}.^2,2));
            dotH = h{instr}*h{instr}';
            for r=nn{instr}
                squaredNormh = sum(h{instr}(r,:).^2);
                squaredDiffh = sum(diff(h{instr}(r,:)).^2);
                for t=1:N
                    if h{instr}(r,t)>hThreshold
                        vloc = w{instr,r}(:).*(Lambda(:,t).^(beta-2));
                                               
                        % Distance function
                        G = sum(vloc.*V(:,t));
                        F = sum(vloc.*Lambda(:,t));
                        
                        if F >0
                            htemp{instr}(r,t) = h{instr}(r,t).* G ./ (F + eps) ;
                        end
                    end
                end
            end
            
        end
    end
    
    h = htemp;
    
    if sum(isnan(h{instr}(:)))
        disp('h is Nan. Function returns')
        return;
    end
    
    
    
    % recomputation of Lambda
    Lambda(:) = 0;
    for instr = 1:nbInstr
        
        if timeVaryingF0(instr) % time-varying fundamental frequency
            
            for t = 1:N
                for r = nn{instr}
                    if h{instr}(r,t)>0
                        Lambda(:,t) = Lambda(:,t)+ w{instr,r}(:,t)*h{instr}(r,t);
                    end
                end
            end
            
        else % non time-varying fundamental frequency
            
            Lambda = Lambda + cell2mat(w(instr,:))*h{instr}(nn{instr},:);
            
        end
    end
    
    LambdaHarmo = Lambda;
    
    Lambda = Lambda + W*H + LambdaOnset;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update of ak
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if freeHarmonicTemplate
        
        % one harmonic distribution for each atom
        if ~sameHarmonicDistribution
            
            muSmoothW = mu.muSmoothW;
            
            LambdaBetam2 = Lambda.^(beta-2);
            
            for instr = 1:nbInstr
                
                if timeVaryingF0(instr) % time-varying fundamental frequency
                    
                    for k = 1:floor(0.5/f0center(1))
                        for r = nn{instr}
                            F = 0;
                            G = 0;
                            Fsmooth = 0;
                            Gsmooth = 0;
                            if sum(h{instr}(r,:)) ~= 0
                                if k<=floor(0.5/f0center(r))
                                    for t=1:N
                                        kf0 = f0{instr}(t,r)*2*M;
                                        deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                                        fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                                        vloc = LambdaBetam2(deb:fin,t).*g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(t,r)) + len/2))*h{instr}(r,t);
                                        
                                        vloc1 = vloc';
                                        F = F + vloc1*Lambda(deb:fin,t);
                                        G = G + vloc1*V(deb:fin,t);
                                        
                                        Fsmooth = Fsmooth + 4*ak{instr}(k,r).^2;
                                        Gsmooth = Gsmooth + 2*ak{instr}(k,r)*ak{instr}(k+1,r);
                                        if k>1
                                            Gsmooth = Gsmooth + 2*ak{instr}(k,r)*ak{instr}(k-1,r);
                                        end
                                    end
                                end
                                if F>0
                                    ak{instr}(k,r) = ak{instr}(k,r).*(G + muSmoothW*Gsmooth)/(F + muSmoothW*Fsmooth + eps);% %+muSmoothW*Gsmooth%+muSmoothW*Fsmooth
                                end
                            end
                        end
                    end
                    if sum(isnan(ak{instr}(:)))>0
                        disp('ak is NaN')
                        return
                    end
                    
                    
                else % non time-varying fundamental frequency
                    
                    for k = 1:floor(0.5/f0center(1))
                        for r = nn{instr}
                            F = 0;
                            G = 0;
                            Fsmooth = 0;
                            Gsmooth = 0;
                            if sum(h{instr}(r,:)) ~= 0
                                kf0 = f0{instr}(r)*2*M;
                                deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                                fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                                if k<=floor(0.5/f0center(r))
                                    for t=1:N
                                        vloc = LambdaBetam2(deb:fin,t).*g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(r)) + len/2))*h{instr}(r,t);
                                        
                                        vloc1 = vloc';
                                        F = F + vloc1*Lambda(deb:fin,t);
                                        G = G + vloc1*V(deb:fin,t);
                                        
                                        Fsmooth = Fsmooth + 4*ak{instr}(k,r).^2;
                                        Gsmooth = Gsmooth + 2*ak{instr}(k,r)*ak{instr}(k+1,r);
                                        if k>1
                                            Gsmooth = Gsmooth + 2*ak{instr}(k,r)*ak{instr}(k-1,r);
                                        end
                                    end
                                end
                                if F > 0
                                    ak{instr}(k,r) = ak{instr}(k,r).*(G + muSmoothW*Gsmooth)/(F + muSmoothW*Fsmooth + eps);% %+muSmoothW*Gsmooth%+muSmoothW*Fsmooth
                                end
                            end
                        end
                    end
                    if sum(isnan(ak{instr}(:)))>0
                        disp('ak is NaN')
                        return
                    end
                    
                end
                
            end
            
        else  % same harmonic distribution for every atom
            
            muSmoothW = mu.muSmoothW;
            LambdaBetam2 = Lambda.^(beta-2);
            onesVect = ones(M,1);
            
            for instr = 1:nbInstr
                
                if timeVaryingF0(instr) % time-varying fundamental frequency
                    for k=1:floor(0.5/f0center(1))
                        F = 0;
                        G = 0;
                        Fsmooth = 0;
                        Gsmooth = 0;
                        
                        for r = nn{instr}
                            if k<=floor(0.5/f0center(r))
                                for t=1:N
                                    if sum(Lambda(:,t)) ~= 0
                                        kf0 = f0{instr}(t,r)*2*M;
                                        deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                                        fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                                        vloc = LambdaBetam2(deb:fin,t).*g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(t,r)) + len/2))*h{instr}(r,t);
                                        
                                        
                                        if isnan(vloc(:))
                                            return
                                        end
                                        vloc1 = vloc';
                                        F = F + vloc1*Lambda(deb:fin,t);
                                        G = G + vloc1*V(deb:fin,t);
                                        
                                        Fsmooth = Fsmooth + 4*ak{instr}(k,r).^2;
                                        Gsmooth = Gsmooth + 2*ak{instr}(k,r)*ak{instr}(k+1,r);
                                        if k>1
                                            Gsmooth = Gsmooth + 2*ak{instr}(k,r)*ak{instr}(k-1,r);
                                        end
                                    end
                                end
                            end
                        end
                        if isnan(F)||isnan(G)
                            return
                        end
                        if F~=0
                            ak{instr}(k,:) = ak{instr}(k,:).* G / (F + eps);
                        end
                        if sum(isnan(ak{instr}(:)))>0
                            disp('ak is NaN')
                            return
                        end
                    end
                    
                else % non time-varying fundamental frequency
                    
                    toBecomputed = zeros(R,1);
                    for r = nn{instr}
                        toBecomputed(r) = sum(h{instr}(r,:))>hThreshold*sum(h{instr}(r,:)>0);
                    end
                    
                    toBecomputedrt = (h{instr}>hThreshold);
                    
                    for k=1:floor(0.5/f0center(1))
                        F = 0;
                        G = 0;
                        Fsmooth = 0;
                        Gsmooth = 0;
                        
                        for r = nn{instr}
                            if k<=floor(0.5/f0center(r))
                                kf0 = f0{instr}(r)*2*M;
                                deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                                fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                                
                                if toBecomputed(r)
                                    peakHarmonic = g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(r)) + len/2));
                                    for t=1:N
                                        if toBecomputedrt(r,t)
                                            vloc =LambdaBetam2(deb:fin,t).*peakHarmonic*h{instr}(r,t);
                                            
                                            vloc1 = vloc';
                                            F = F + vloc1*Lambda(deb:fin,t);
                                            G = G + vloc1*V(deb:fin,t);
                                        end
                                        
                                    end
                                    Fsmooth = Fsmooth + 4*ak{instr}(k,r).^2;
                                    Gsmooth = Gsmooth + 2*ak{instr}(k,r)*ak{instr}(k+1,r);
                                    if k>1
                                        Gsmooth = Gsmooth + 2*ak{instr}(k,r)*ak{instr}(k-1,r);
                                    end
                                end
                            end
                        end
                        if F~=0
                            ak{instr}(k,:) = ak{instr}(k,:).* G / (F + eps) ;
                        end
                    end
                end
                
            end
            
        end
        
        
        
        
        % normalization of ak
        for instr = 1:nbInstr
            for r = nn{instr}
                sumak = sum(ak{instr}(1:floor(0.5/f0center(r)),r));
                ak{instr}(:,r) = ak{instr}(:,r)/sumak;
                h{instr}(r,:) = h{instr}(r,:)*sumak;
            end
        end
        
        % recomputation of w and Lambda
        Lambda(:) = 0;
        for instr = 1:nbInstr
            
            if timeVaryingF0(instr) % time-varying fundamental frequency
                
                for r = nn{instr}
                    w{instr,r}(:) = 0;
                    for t = 1:N
                        if h{instr}(r,t)>hThreshold
                            kf0 = f0{instr}(t,r)*2*M;
                            for k=1:floor(0.5/f0center(r))
                                deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                                fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                                w{instr,r}(deb:fin,t) = w{instr,r}(deb:fin,t) + ak{instr}(k,r) * g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(t,r)) + len/2));
                            end
                            w{instr,r}(:,t) = w{instr,r}(:,t) + wThreshold;
                            Lambda(:,t) = Lambda(:,t)+ w{instr,r}(:,t)*h{instr}(r,t);
                        end
                    end
                end
                
            else % non time-varying fundamental frequency
                
                for r = nn{instr}
                    w{instr,r}(:) = 0;
                    if sum(h{instr}(r,:))>hThreshold*sum(h{instr}(r,:)>0)
                        kf0 = f0{instr}(r)*2*M;
                        for k=1:floor(0.5/f0center(r))
                            deb = max(floor(k*sqrt(1+k^2*Binharm{instr}(r))*kf0-Deltamr(r)),1);
                            fin = min(ceil(k*sqrt(1+k^2*Binharm{instr}(r))*kf0+Deltamr(r)),M);
                            w{instr,r}(deb:fin) = w{instr,r}(deb:fin) + ak{instr}(k,r) * g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm{instr}(r))*f0{instr}(r)) + len/2));
                        end
                        w{instr,r} = w{instr,r} + wThreshold;
                    end
                end
                Lambda = Lambda + cell2mat(w(instr,:))*h{instr}(nn{instr},:);
                
            end
            
        end
        
        LambdaHarmo = Lambda;
        Lambda = Lambda + W*H + LambdaOnset;
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % non harmonic, non time-varying template update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    muSmoothWNonHarmonic = 100;
    
    indexInstrumentNoise = 1;
    for instr = 1:nbInstr
        if initVal.associatedTemplateNoise(instr)
            H(indexInstrumentNoise,:) = sum(h{instr});
            indexInstrumentNoise = indexInstrumentNoise + 1;
        end
    end
    Lambda = LambdaHarmo + W*H + LambdaOnset;
    
    
    % update of non harmonic patterns
    W = W.*((Lambda.^(beta-2).*V)*H' + 2*muSmoothWNonHarmonic*W.*(shiftUD(W,1)+shiftUD(W,-1)))./(Lambda.^(beta-1)*H' + 4*muSmoothWNonHarmonic*W.^2 + eps);
    
    % recomputation of Lambda
    Lambda = LambdaHarmo + W*H + LambdaOnset;
        
    
    
    
    
    
    % update of H for each value of t (which will be averaged)
    for t=0:Tonset-1
        Ht(t+1,:,:) = Honset.* (reshape(Wonset(t+1,:,:),M,Ronset)'*shiftLR(V./(Lambda + eps),t))./(reshape(Wonset(t+1,:,:),M,Ronset)'*One + eps);
    end
    
    % average along t
    Honset = reshape(mean(Ht),Ronset,N);
    
    
    % recomputation of Lambda
    LambdaOnset(:) = 0;
    for t=0:Tonset-1
        LambdaOnset = LambdaOnset + reshape(Wonset(t+1,:,:),M,Ronset)* shiftLR(Honset,-t);
    end
    
    Lambda = LambdaHarmo + W*H + LambdaOnset;
    
    % update of Wt
    for t=0:Tonset-1
        Wonset(t+1,:,:) = reshape(Wonset(t+1,:,:),M,Ronset).*(  ((V./(Lambda + eps)) * shiftLR(Honset,-t)') ./ (One*shiftLR(Honset,-t)' + eps));
    end
    
    % computation of Lambda
    LambdaOnset(:) = 0;
    for t=0:Tonset-1
        LambdaOnset = LambdaOnset + reshape(Wonset(t+1,:,:),M,Ronset)* shiftLR(Honset,-t);
    end
    
    Lambda = LambdaHarmo + W*H + LambdaOnset;
    
    % plotting results
    if plotResults
        
        % plot activations
        hplot = h{1};
        for instr = 2:nbInstr
            hplot = hplot + h{instr};
        end
    end
    
    
    disp(['iteration ' int2str(iter)])
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% annex functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stable] = isStable(A);

N = length(A)-1; % Order of A(z)
stable = 1;      % stable unless shown otherwise
A = A(:);        % make sure it's a column vector
for i=N:-1:1
    rci=A(i+1);
    if abs(rci) >= 1
        stable=0;
        return;
    end
    A = (A(1:i) - rci * A(i+1:-1:2))/(1-rci^2);
end


function [hn] = minphase(h)

r = cplxpair(roots(h));

% Find roots inside & on unit circle

ra = abs(r);

iru   = find(abs(ra-1) <= 1e-8); % indices for roots on uc
irin  = find(ra-1 < -1e-8);  % indices for roots inside uc
irout = find(ra-1 >  1e-8);  % indices for roots outside uc

% Map roots outside the unit circle to inside:

cr = [r(iru) ; r(irin) ; 1./conj(r(irout))];

hn = h(1)*prod(abs(r(irout)))*poly(cr);



function g = squaredWindow(f)

T = 1024;
a = 0.5;
b = 0.5;
w = 2*pi*f;
g = (2-2*cos(T*w)).*(T^2*w.^2*(b-a)+4*pi^2*a).^2./(T^2*w.^3 - 4*pi^2*w).^2;

g(abs(f)<10^-10) = a^2*T^2;
g(abs(f-1/T)<10^-10) = b^2*T^2/4;
g(abs(f+1/T)<10^-10) = b^2*T^2/4;


function dg = derivativeWindowOnf(f)

a = 0.5;
b = 0.5;
T = 1024;

derivat = 1/(4*pi^2)*( (2-2*cos(2*pi*T*f)).*2.*(f.^2*T^2*(b-a)+a)./(f.^2.*(f.^2*T^2-1).^2)   .*...
    ( 2*f*T^2*(b-a) - (2*f*T^2.*(f.^2*T^2*(b-a)+a))./(f.^2*T^2 - 1) -  (f.^2*T^2*(b-a)+a)./f) + ...
    4*pi*T*sin(2*pi*T*f) .* (T^2*f.^2*(b-a) + a).^2./(f.^2.*(T^2*f.^2-1).^2) );

dg = derivat./f;

f0 = 0.000001;
derivat0 = 1/(4*pi^2)*( (2-2*cos(2*pi*T*f0)).*2.*(f0.^2*T^2*(b-a)+a)./(f0.^2.*(f0.^2*T^2-1).^2)   .*...
    ( 2*f0*T^2*(b-a) - (2*f0*T^2.*(f0.^2*T^2*(b-a)+a))./(f0.^2*T^2 - 1) -  (f0.^2*T^2*(b-a)+a)./f0) + ...
    4*pi*T*sin(2*pi*T*f0) .* (T^2*f0.^2*(b-a) + a).^2./(f0.^2.*(T^2*f0.^2-1).^2) );

dg0 = derivat0/f0;
dg(f==0) = dg0;


function bD = betaDiv(V,Vh,beta)
if beta == 0
    bD = sum((V(:)./Vh(:))-log(V(:)./Vh(:)) - 1);
elseif beta == 1
    bD = sum(V(:).*(log(V(:))-log(Vh(:))) + Vh(:) - V(:));
else
    bD = sum(max(1/(beta*(beta-1))*(V(:).^beta + (beta-1)*Vh(:).^beta - beta*V(:).*Vh(:).^(beta-1)),0));
end





function bD = betaDivInharm(Binharm,f0,h,ak,V,VhatMinusVhatr,beta,g,Deltamr,f0center,freq,len,hThreshold,wThreshold,division,M,N,r)

wr = zeros(M,1);
Vhatr = zeros(M,N);

if sum(h(r,:))>hThreshold*sum(h(r,:)>0)
    kf0 = f0(r)*2*M;
    for k=1:floor(0.5/f0center(r))
        deb = max(floor(k*sqrt(1+k^2*Binharm)*kf0-Deltamr(r)),1);
        fin = min(ceil(k*sqrt(1+k^2*Binharm)*kf0+Deltamr(r)),M);
        wr(deb:fin) = wr(deb:fin) + ak(k,r) * g(round(division*(freq(deb:fin)-k*sqrt(1+k^2*Binharm)*f0(r)) + len/2));
    end
    wr = wr + wThreshold;
end


Vhatr = wr(:)*h(r,:);



Vh = VhatMinusVhatr + Vhatr;
if beta == 0
    bD = sum((V(:)./Vh(:))-log(V(:)./Vh(:)) - 1);
elseif beta == 1
    bD = sum(V(:).*(log(V(:))-log(Vh(:))) + Vh(:) - V(:));
else
    bD = sum(max(1/(beta*(beta-1))*(V(:).^beta + (beta-1)*Vh(:).^beta - beta*V(:).*Vh(:).^(beta-1)),0));
end
