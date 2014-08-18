function Ashifted = shiftLR(A,k)
% B = shiftLR(A,k)
%
%   Shifts the values of the matrix A by k elements to the left (setting to
%   zero last values).
%   If k<0, shifts the values by -k to the right (setting to zero first values).
%
%           example :
%               >>A = [1 2 3 4; 5 6 7 8];
%               >>shiftLR(A,2)
%               ans =
%                   3     4     0     0
%                   7     8     0     0

Ashifted = circshift(A,[0 -k]);
if k>=0
    Ashifted(:,max(end-k+1,1):end) = 0;
else
    Ashifted(:,1:min(-k,end)) = 0;
end