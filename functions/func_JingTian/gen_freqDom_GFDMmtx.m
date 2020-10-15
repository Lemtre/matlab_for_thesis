function G = gen_freqDom_GFDMmtx(g, K, M, group_order)
% GEN_FREQDOM_GFDM_OQAMMTX generate the frequency domain GFDM modulation
%   matrix for a given prototype filter "g";
% 
%   g:           the prototype filter in time domain;
%   K:           number of subcarriers;
%   M:           number of subsymbols;
%   group_order: group data into subcarriers ('sc') or subsymbols ('ss')
%   G:           frequency domain modulation matrix.
%
% Author: Jing Tian
% Last modified: June 26th, 2020.
if nargin < 4
    group_order = 'sc';
end
g = g(:);
N = length(g);
if N~=K*M
    error('Filter length mismatch!');
end

gf = fft(g);
gf = gf./norm(gf);
shift = 0:M:N-1;
Smtx = arrayfun(@(s) circshift(gf, s), shift, 'Uniformoutput', false);
mm = 0:M-1;
nn = 0:N-1;
freq = exp(-1j*2*pi*nn.'*mm/M);
G = cell2mat(cellfun(@(gf_s) diag(gf_s)*freq, Smtx, 'UniformOutput', false));

if strcmpi(group_order, 'ss')
    idx = 1:N;
    idx = reshape(idx, M, K).';
    idx_ss = idx(:);
    G = G(:,idx_ss);
end

end