function xMono = toMono(x,method)
% xMono = toMono(x,method)
%   transforms x (a mono or stereo to signal) to a mono signal (if x is already mono it does nothing)
% input :
%       - x : signal to transform (can be a 2xn, nx2, 1xn or nx1 signal)
%       - method : method used for transformation. possible values are :
%                   - 'sum' : sums channels (default method)
%                   - '1' : keep only first channel
%                   - '2' : keep only second channel
%
% Copyright (c) 2010 Romain Hennequin (http://www.romain-hennequin.fr/).

if nargin == 1
    method = 'sum';
elseif nargin>2 || nargin == 0
    error('wrong number of input arguments, type help toMono');
end

xMono = x;

[nbChannel dimChannel] = min(size(x));
if nbChannel == 2
    if strcmp(method,'1')
        if dimChannel == 1
            xMono = x(1,:);
        else
            xMono = x(:,1);
        end
    elseif strcmp(method,'2')
        if dimChannel == 1
            xMono = x(2,:);
        else
            xMono = x(:,2);
        end
    else
        if strcmp(method,'sum') == 0

            warning('this method does not exist, type help toMono. default method used instead');
        end
        xMono = sum(x,dimChannel);
    end
elseif nbChannel>2
    error('too many channels : toMono convert only stereo signal to mono signal');
end
