function [mf_out, mf_i, mf_q] = mf_gfdm_oqam...
    (g, K, M, group_order, oqam_method, sig)
% MF_GFDM_OQAM matched filtering (MF) demodulation of GFDM/OQAM signal.
% 
%   g:           the prototype filter in time domain;
%   K:           number of subcarriers;
%   M:           number of subsymbols;
%   group_order: group data into subcarriers ('sc') or subsymbols ('ss')
%   oqam_method: 'SMT' or 'CMT'.
%   sig:         received signal.
%   mf_output:   take the real parts of I/Q outputs and combine into a
%                complex signal;
%   mf_i:        the in-phase (I) complex output;
%   mf_q:        the quadrature (Q) complex output;
%
% Author: Jing Tian
% Last modified: June 26th, 2020.
sig = sig(:);
g = g(:);
N = length(g);
if N~=K*M
    error('Filter length mismatch.');
elseif N~=length(sig)
    error('Filter length and signal length mismatch.');
elseif mod(K,2)~=0 || mod(M,2)~=0 
    error('Numbers of subsymbols and subcarriers must be even.');
end

switch lower(oqam_method)
    case 'smt'
        kk = 0:K-1;
        phase_i = repmat((1j).^kk.', M, 1);
        phase_q = 1j*phase_i;
        mf_i = fast_gfdm_mf(g, K, M, 'ss', sig);
        mf_i = mf_i .* conj(phase_i);
        g_q = circshift(g, K/2);
        mf_q = fast_gfdm_mf(g_q, K, M, 'ss', sig);
        mf_q = mf_q .* conj(phase_q);
    case 'cmt'
        mm = 0:M-1;
        phase_i = repmat((1j).^mm, K, 1);
        phase_i = phase_i(:);
        phase_q = 1j*phase_i;
        mf_i = fast_gfdm_mf(g, K, M, 'ss', sig);
        mf_i = mf_i .* conj(phase_i(:));
        nn = 0 : N-1;
        mf_q = exp(-1j*pi*nn.'/K) .* sig;
        mf_q = fast_gfdm_mf(g, K, M, 'ss', mf_q);
        mf_q = mf_q .* conj(phase_q(:));
    otherwise
        error('Invalid method.');
end

if strcmpi(group_order, 'sc')
    [mf_i, mf_q] = deal(reshape(reshape(mf_i,K,M).',[],1),...
        reshape(reshape(mf_q,K,M).',[],1));
end
mf_out = real(mf_i) + 1j*real(mf_q);

end

