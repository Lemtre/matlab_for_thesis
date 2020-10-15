function gfdm_sig = fast_gfdm_mod(g, K, M, group_order, data)
% FAST_GFDM_MOD Fast algorithm of GFDM modulation based on
%   Po-Chih Chen's paper, DOI: 10.1109/TSP.2017.2718971.
% 
%   g:           prototype filter.
%   K:           number of subcarriers.
%   M:           number of subsymbols.
%   group_order: data grouped in subcarriers ('sc') or subsymbols ('ss'). 
%   data:        baseband complex symbols.
%   gfdm_sig:    modulated GFDM signal.
% 
% Author: Jing Tian
% Last modified: 9th April, 2020.

if nargin < 5
    group_order = 'ss';
    warning('The input data is treated as subsymbol-grouped by default.');
end
g = g(:);
N = length(g);
if N~=K*M
    error('Filter length mismatch.');
elseif numel(data)~=N;
    error('Filter length and data length mismatch.');
end
g = g./norm(g); % normalize the prototype filter g

if strcmpi(group_order, 'sc')
    data = reshape(data(:), M, K).';
end
G = fft(reshape(g, K, M).');
data_ss = reshape(data(:), K, M);
Dk = ifft(data_ss);
Dm = fft(Dk.');
Dg = Dm.*G;
Dout = ifft(Dg).';
gfdm_sig = K*Dout(:); % This last step of normalization is necessary for OQAM modulation.

end

