function smoothedSignal = smooth(signal,N)
% smoothedSignal = smooth(signal,N)
% Filtering with N point moving average.
% Copyright (c) 2010 Romain Hennequin (http://www.romain-hennequin.fr/).

f = ones(N,1)/N;
smoothedSignal = filter(f,1,signal);
