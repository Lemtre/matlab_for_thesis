function freq_gfdm_zf = zf_FreqDom_GFDM(g, K, M, group_order, sig)
% MOD_FREQDOM_GFDM frequency domain zero-forcing (ZF) GFDM demodulation.
% 
%   g:           the prototype filter in time domain;
%   K:           number of subcarriers;
%   M:           number of subsymbols;
%   group_order: group data into subcarriers ('sc') or subsymbols ('ss')
%   sig:         received signal;
%   freq_gfdm_zf: frequency domain ZF demodulation output.
%
% Author: Jing Tian
% Last modified: June 26th, 2020.
if nargin < 5
    group_order = 'ss';
    warning('The output data is treated as subsymbol-grouped by default.');
end
g = g(:);
N = length(g);
if N~=K*M
    error('Filter length mismatch.');
elseif numel(sig)~=N;
    error('Filter length and data length mismatch.');
end
sig_K_M = reshape(sig, M, K).';
Dm = ifft(sig_K_M);

gf = fft(g);
gf = gf./norm(gf);
G = conj(fft(reshape(gf, M, K)'));

Dg = Dm./G;
Dk = fft(Dg);
xss = ifft(Dk.').';

if strcmpi(group_order, 'sc')
    xsc = xss.';
    freq_gfdm_zf = xsc(:);
else
    freq_gfdm_zf = xss(:);
end
        
end