function [x_oqam, x_i, x_q] = mod_gfdm_oqam(g, K, M,...
    group_order, oqam_method, data_cplx)
% MOD_GFDM_OQAM modulation of GFDM/OQAM signal (fast algorithm).
% 
%   g:           the prototype filter in time domain;
%   K:           number of subcarriers;
%   M:           number of subsymbols;
%   group_order: group data into subcarriers ('sc') or subsymbols ('ss')
%   oqam_method: 'SMT' or 'CMT'.
%   data_cplx:   complex baseband signal.
%   x_oqam:      the final GFDM/OQAM signal in time-domain;
%   x_i:         the in-phase (I) complex signal in time domain;
%   x_q:         the quadrature (Q) complex signal in time domain;
%
% Author: Jing Tian
% Last modified: June 26th, 2020.

g = g(:);
N = length(g);
if N~=K*M
    error('Filter length mismatch.');
elseif mod(K,2)~=0 || mod(M,2)~=0 
    error('Numbers of subsymbols and subcarriers must be even.');
end
g = g./norm(g);

if strcmpi(group_order, 'sc')
    data_cplx = reshape(reshape(data_cplx(:,1),M,K).',[],1);
elseif ~strcmpi(group_order, 'ss')
    error('Invalid data format.');
end
[d_i, d_q] = deal(real(data_cplx), imag(data_cplx));


switch lower(oqam_method)
    case 'smt'
        kk = 0:K-1;
        phase_i = repmat((1j).^kk.', M, 1);
        phase_q = 1j*phase_i;
        d_i = d_i .* phase_i;
        d_q = d_q .* phase_q;
        x_i = fast_gfdm_mod(g, K, M, 'ss', d_i);
        g_q = circshift(g, K/2);
        x_q = fast_gfdm_mod(g_q, K, M, 'ss', d_q);
    case 'cmt'
        mm = 0:M-1;
        phase_i = repmat((1j).^mm, K, 1);
        phase_i = phase_i(:);
        phase_q = 1j*phase_i;
        d_i = d_i .* phase_i;
        d_q = d_q .* phase_q;
        x_i = fast_gfdm_mod(g, K, M, 'ss', d_i);
        x_q = fast_gfdm_mod(g, K, M, 'ss', d_q);
        nn = 0 : N-1;
        x_q = exp(1j*pi*nn.'/K) .* x_q;
    otherwise
        error('Invalid method.');
end
x_oqam = x_i + x_q;

end

