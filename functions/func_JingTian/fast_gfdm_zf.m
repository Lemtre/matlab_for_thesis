function x = fast_gfdm_zf(g, K, M, group_order, sig)
% FAST_GFDM_MOD Fast algorithm of GFDM zero-forcing (ZF) demodulation based on
%   Po-Chih Chen's paper, DOI: 10.1109/TSP.2017.2718971.
% 
%   g:           prototype filter.
%   K:           number of subcarriers.
%   M:           number of subsymbols.
%   group_order: data grouped in subcarriers ('sc') or subsymbols ('ss'). 
%   sig:         received signal.
%   x:           ZF demodulated GFDM signal.
% 
%   Data belonging to the same subsymbols are grouped together by default
%       if 'group_order' is omitted. 
% 
% Author: Jing Tian
% Last modified: 9th April, 2020.
if nargin < 5
    group_order = 'ss';
    warning('The output data is subsymbol-grouped by default.');
end
g = g(:);
N = length(g);
if N~=K*M
    error('Filter length mismatch.');
elseif numel(sig)~=N;
    error('Filter length and signal length mismatch.');
end
g = g./norm(g);

G = fft(reshape(g, K, M).');
Dm = fft(reshape(sig(:), K, M).');
Dg = Dm./G;
Dk = fft(Dg.');
xss = ifft(Dk.').';
xss = xss./sqrt(N*M);

if strcmpi(group_order, 'sc')
    xsc = xss.';
    x = xsc(:);
else
    x = xss(:);
end

end

