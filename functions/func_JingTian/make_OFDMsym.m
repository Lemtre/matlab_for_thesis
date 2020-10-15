function x_ofdm = make_OFDMsym(ofdm, data, pilots, varargin)
% MAKE_OFDMSYM synthesize OFDM symbols using IFFT.
%   
%   sym:    data symbols.
%   pilots: pilots.
%   ofdm:   structure containing OFDM parameters,
%           ofdm.idd - positions of data symbols;
%           ofdm.idp - positions of data pilots;
%           ofdm.idz - positions of null carriers;
%           ofdm.Nc  - number of subcarriers per block (number of IFFT points).
%           ofdm.window - transmission window 
%   x_ofdm: OFDM symbol in baseband (without CP/ZP).
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

if nargin == 3
    shiftOp = true;
elseif nargin > 5
    error('The number of optional parameters should be no greater than two.');
elseif nargin == 4
    error('Optional parameters should always go by pairs.');
elseif ~strcmp(varargin{1}, 'shiftOp')
    prop = varargin{1};
    error(['There is no ' prop ' property of this function.']);
else
    shiftOp = varargin{2};
end

data = data(:);
pilots = pilots(:);
L = numel(ofdm.idp)+numel(ofdm.idd)+numel(ofdm.idz);
if numel(data)~=numel(ofdm.idd) || numel(pilots)~=numel(ofdm.idp)
    error('Data length mismatch.');
elseif numel(union(union(ofdm.idp, ofdm.idd), ofdm.idz))~=L
    error('Positions conflict.');
elseif L > ofdm.Nc
    error('IFFT size should be made greater.');
end
Xf = zeros(ofdm.Nc,1);
Xf(ofdm.idd) = data;
Xf(ofdm.idp) = pilots;
if true(shiftOp)
    Xf = ifftshift(Xf);
end
x_ofdm = ifft(Xf);%.*sqrt(ofdm.Nc);
if length(ofdm.window) > ofdm.Nc % windowing if specified 
    N_extra = length(ofdm.window) - ofdm.Nc;
    x_ofdm = [x_ofdm; x_ofdm(1:N_extra)];
    x_ofdm = x_ofdm .* ofdm.window;
end


end

