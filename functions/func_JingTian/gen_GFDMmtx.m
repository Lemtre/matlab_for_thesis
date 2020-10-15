function A = gen_GFDMmtx(g, K, M, group_order)
% GEN_GFDMMTX generate the GFDM modulation matrix for a given prototype
%   filter "g";
%   g:           prototype filter in time domain.
%   K:           number of subcarriers;
%   M:           number of subsymbols;
%   group_order: grouping the data vector into subsymbols (ss) or subcarriers (sc)
%   A:           the GFDM modulation matrix.
% 
% Author: Jing Tian
% Last modified: Dec 11th, 2019.
if nargin < 4
    group_order = 'sc';
end
g = g(:);
N = length(g);
if N~=K*M
    error('Filter length mismatch!');
end
g = g./norm(g);

shift = 0:K:N-1;
Smtx = arrayfun(@(s) circshift(g, s), shift, 'Uniformoutput', false);
kk = 0:K-1;
nn = 0:N-1;
freq = exp(1j*2*pi*nn.'*kk/K);
A = cell2mat(cellfun(@(gf_s) diag(gf_s)*freq, Smtx, 'UniformOutput', false));

if strcmpi(group_order, 'sc')
    idx = 1:N;
    idx = reshape(idx, K, M).';
    idx_sc = idx(:);
    A = A(:,idx_sc);
end

end