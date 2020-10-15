function freqDom_gfdm_sig = mod_FreqDom_GFDM(g, K, M, group_order, data)
% MOD_FREQDOM_GFDM frequency domain GFDM modulation.
% 
%   g:           the prototype filter in time domain;
%   K:           number of subcarriers;
%   M:           number of subsymbols;
%   group_order: group data into subcarriers ('sc') or subsymbols ('ss')
%   data:        baseband data symbols;
%   freqDom_gfdm_sig: frequency domain GFDM signal.
%
% Author: Jing Tian
% Last modified: June 26th, 2020.
if nargin < 5
    group_order = 'sc';
    warning('The input data is treated as subsymbol-grouped by default.');
end
g = g(:);
N = length(g);
if N~=K*M
    error('Filter length mismatch.');
elseif numel(data)~=N;
    error('Filter length and data length mismatch.');
end

gf = fft(g)./sqrt(N);
if strcmpi(group_order, 'ss')
    data = reshape(data(:), K, M).';
    data = data(:);
end
data_conj = conj(data);
gf_conj = conj(gf);
freq_gfdm_sig_conj = fast_gfdm_mod(gf_conj, M, K, 'ss', data_conj);
freqDom_gfdm_sig = conj(freq_gfdm_sig_conj);

end