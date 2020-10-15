function freq_gfdm_mf = mf_FreqDom_GFDM(g, K, M, group_order, sig)
% MF_FREQDOM_GFDM frequency domain matched filtering demodulation of
%                 GFDM signal (fast algorithm).
% 
%   g:            the prototype filter in time domain;
%   K:            number of subcarriers;
%   M:            number of subsymbols;
%   group_order:  group data into subcarriers ('sc') or subsymbols ('ss')
%   sig:          received signal.
%   freq_gfdm_mf: matched filtering output.
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
elseif numel(sig)~=N;
    error('Filter length and data length mismatch.');
end

gf = fft(g);
gf = gf./norm(gf);
gf_conj = conj(gf);
sig_conj = conj(sig);
freq_gfdm_mf = fast_gfdm_mf(gf_conj, M, K, 'ss', sig_conj);
freq_gfdm_mf = conj(freq_gfdm_mf);

if strcmpi(group_order, 'ss')
    freq_gfdm_mf = reshape(freq_gfdm_mf, M, K).';
    freq_gfdm_mf = freq_gfdm_mf(:);
end
        
end