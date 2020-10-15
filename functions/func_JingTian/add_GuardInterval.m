function Xg = add_GuardInterval(X, Ng, type_GI)
% add_GuardInterval add guard intervals to data blocks.
%   
%   X:       matrix whose columns representing different data blocks.
%   Ng:      number of samples per guard interval.
%   type_GI: string, either 'cp' or 'zp'.
%   Xg:      guarded output.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

if ~ismatrix(X)
    error('The first input should be a matrix.');
end

switch lower(type_GI)
    case 'cp'
        Xg = [X(end-Ng+1 : end, :); X];
    case 'zp'
        Xg = [X; zeros(Ng, size(X,2))];
    otherwise
        error('Invalid guard interval type.');
end

end

