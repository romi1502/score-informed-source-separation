function Ashifted = shiftUD(A,k)
% B = shiftUD(A,k)
%
%   Shifts the values of the matrix A by k elements down (setting to
%   zero last values).
%   If k<0, shifts the values by -k up (setting to zero first values).
%
%           example :
%               >>A = [1 2;3 4; 5 6 ;7 8];
%               >>shiftUD(A,2)
%               ans =
%                   3     4     0     0
%                   7     8     0     0



Ashifted = circshift(A,[k 0]);
if k<=0
    Ashifted(max(end+k+1,1):end,:) = 0;
else 
    Ashifted(1:min(k,end),:) = 0;
end