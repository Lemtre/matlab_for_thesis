function Xf = ser2para_rmvGI_fft(sig, Ncyc, Nblk, typeGI, varargin)
% SER2PARA_RMVDI_FFT convert serially block-transmitted signal into
% parallel, remove the guard invertals, and perform FFT to each block.
% 
%   sig:      received signal.
%   Ncyc:     number of samples per cycle.
%   Nblk:     number of blocks.
%   typeGI:   string, either 'cp' or 'zp'.
%   varargin: 0 or 2 optional inputs, specifying whether to perform
%             fftshift.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

if nargin == 4
    shiftOp = true;
elseif nargin > 6
    error('The number of optional parameters should be no greater than two.');
elseif nargin == 5
    error('Optional parameters should always go by pairs.');
elseif ~strcmp(varargin{1}, 'shiftOp')
    prop = varargin{1};
    error(['There is no ' prop ' property of this function.']);
else
    shiftOp = varargin{2};
end

x_para = reshape(sig, [], Nblk);
Ncol = size(x_para, 1);
Ng = Ncol - Ncyc;
switch lower(typeGI)
    case 'cp'
        x_data = x_para(Ng+1:end, :);
    case 'zp'
        x_data = x_para(1:end-Ng, :);
        x_data(1:Ng, :) = x_data(1:Ng, :) + x_para(end-Ng+1:end, :);
    otherwise
        error('Invalid guard interval.');
end
Xf = fft(x_data);
if true(shiftOp)
    Xf = fftshift(Xf, 1);
end

end