function [x_oqam, x_i, x_q] = mod_FreqDom_GFDM_OQAM...
    (g, K, M, group_order, oqam_method, data_cplx)
% MOD_FREQDOM_GFDM_OQAM frequency domain modulation of GFDM/OQAM signal.
% 
%   g:           the prototype filter in time domain;
%   K:           number of subcarriers;
%   M:           number of subsymbols;
%   group_order: group data into subcarriers ('sc') or subsymbols ('ss')
%   oqam_method: 'SMT' or 'CMT'.
%   data_cplx:   complex baseband signal.
%   x_oqam:      DFT of the final GFDM/OQAM signal;
%   x_i:         the in-phase (I) complex signal in frequency domain;
%   x_q:         the quadrature (Q) complex signal in frequency domain;
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
if strcmpi(group_order, 'ss')
    data_cplx = reshape(reshape(data_cplx(:,1),K,M).',[],1);
elseif ~strcmpi(group_order, 'sc')
    error('Invalid data format.');
end
[d_i, d_q] = deal(real(data_cplx), imag(data_cplx));

gf = fft(g);
gf = gf./norm(gf);
switch lower(oqam_method)
    case 'smt'
        kk = 0:K-1;
        phase_i = repmat((1j).^kk, M, 1);
        phase_q = repmat((1j).^(3*kk+1), M, 1);
        d_i = d_i .* phase_i(:);
        d_q = d_q .* phase_q(:);
        x_i = negFreq_gfdm_mod(gf, M, K, d_i);
        x_q = negFreq_gfdm_mod(gf, M, K, d_q);
        nn = 0:N-1;
        x_q = exp(-1j*pi*nn.'/M) .* x_q;
    case 'cmt'
        mm = 0:M-1;
        phase_i = repmat((1j).^mm.', K, 1);
        phase_q = repmat((1j).^(3*mm.'+1), K, 1);
        d_i = d_i .* phase_i;
        d_q = d_q .* phase_q;
        x_i = negFreq_gfdm_mod(gf, M, K, d_i);
        gf_q = circshift(gf, M/2);
        x_q = negFreq_gfdm_mod(gf_q, M, K, d_q);
    otherwise
        error('Invalid method.');
end
x_oqam = x_i + x_q;

end

function x = negFreq_gfdm_mod(g, Nsc, Nss, d)

d_conj = conj(d);
g_conj = conj(g);
x_conj = fast_gfdm_mod(g_conj, Nsc, Nss, 'ss', d_conj);
x = conj(x_conj);

end

