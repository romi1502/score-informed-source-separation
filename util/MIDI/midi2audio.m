function [y,Fs] = midi2audio(input,Fs,synthtype,instrTrack)
% y = midi2audio(input, Fs, synthtype,instrTrack)
% y = midi2audio(input, Fs, synthtype)
% y = midi2audio(input, Fs)
% y = midi2audio(input)
%
% Convert midi structure to a digital waveform
%
% Inputs:
%  input - can be one of:
%    a structure: matlab midi structure (created by readmidi.m)
%    a string: a midi filename
%    other: a 'Notes' matrix (as ouput by midiInfo.m)
%
%  synthtype - string to choose synthesis method
%      passed to synth function in synth.m
%      current choices are: 'fm', 'sine' or 'saw'
%      default='fm'
%
%  Fs - sampling frequency in Hz (beware of aliasing!)
%       default =  44.1e3
%
%  instrTrack - cell array of string. These string are the name of the
%       instrument to synthesize for each track. (see the help of synth.m for a list of the supported instruments)

% Copyright (c) 2009 Ken Schutte
% more info at: http://www.kenschutte.com/midi

if (nargin<2)
  Fs=44.1e3;
end
if (nargin<3)
  synthtype='fm';
end
if (nargin<4 &&  synthtype == 'multi')
    error('need fourth argument for multi synthtype')
end

endtime = -1;
if (isstruct(input))
  [Notes,endtime] = midiInfo(input,0);
elseif (ischar(input))
  [Notes,endtime] = midiInfo(readmidi(input), 0);
else
  Notes = input;
end

% t2 = 6th col
if (endtime == -1)
  endtime = max(Notes(:,6));
end
if (length(endtime)>1)
  endtime = max(endtime);
end


y = zeros(1,ceil(endtime*Fs));

for i=1:size(Notes,1)

  f = midi2freq(Notes(i,3));
  dur = Notes(i,6) - Notes(i,5);
  velocity = Notes(i,4);
  amp = velocity/127;


  if strcmp(synthtype,'multi')    
    yt = synth(f, dur, amp, Fs, instrTrack{Notes(i,1)});
  else
    yt = synth(f, dur, amp, Fs, synthtype);  
  end

  n1 = floor(Notes(i,5)*Fs)+1;
  N = length(yt);  

  % ensure yt is [1,N]:
  y(n1:n1+N-1) = y(n1:n1+N-1) + reshape(yt,1,[]);

end