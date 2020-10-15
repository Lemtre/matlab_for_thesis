function cfr_at_pilots = chEst_ofdm_LS(ofdm_fft, pilots, idp, varargin)
% CHEST_LS perform least square channel estimation at pilot positions.
% 
%   ofdm_fft: baseband OFDM block after FFT.
%   pilots:   pilots.
%   idp:      pilot positions.
%   varargin: 0 or 2 optional inputs, specifying whether to perform
%             fftshift.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

if nargin == 3
    shiftOp = false;
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

ofdm_fft = ofdm_fft(:);
pilots = pilots(:);
if true(shiftOp)
    ofdm_fft = fftshift(ofdm_fft, 1);
end
cfr_at_pilots = ofdm_fft(idp)./pilots;

end