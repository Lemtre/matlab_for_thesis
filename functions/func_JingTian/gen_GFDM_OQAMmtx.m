function [Ai, Aq] = gen_GFDM_OQAMmtx(g, K, M, group_order, oqam_method)
% GEN_GFDM_OQAMMTX generate the GFDM/OQAM modulation matrix for a given prototype
%   filter "g";
%   g:           prototype filter in time domain.
%   K:           number of subcarriers;
%   M:           number of subsymbols;
%   group_order: group the data vector into subsymbols (ss) or subcarriers (sc)
%   oqam_method: 'SMT' or 'CMT'.
%   Ai:          the in-phase modulation matrix.
%   Aq:          the quadrature modulation matrix.
% 
% Author: Jing Tian
% Last modified: Dec 11th, 2019.
    if nargin < 5
        oqam_method = 'SMT';
    end
    
    g = g(:);
    N = length(g);
    if N~=K*M
        error('Filter length mismatch!');
    end
    g = g./norm(g);
   
    switch lower(oqam_method)
        case 'smt'
            time_shift = @(s_) circshift(g, s_);
            shift_i = 0:K:N-1;
            shift_q = shift_i + floor(K/2);
            Gi = cell2mat(arrayfun(time_shift, shift_i, 'UniformOutput', false));
            Gq = cell2mat(arrayfun(time_shift, shift_q, 'UniformOutput', false));
            freq = 0: 1/K :1-1/K;
            E = num2cell(exp(1j*2*pi*(0:N-1)'*freq), 1);
            pha = num2cell((1j).^(0:K-1), 1);
            Ai = cell2mat(cellfun(@(pha_k, E_k) pha_k*diag(E_k)*Gi, pha, E, 'UniformOutput', false));
            Aq = cell2mat(cellfun(@(pha_k, E_k) (1j)*pha_k*diag(E_k)*Gq, pha, E, 'UniformOutput', false));
        case 'cmt'
            shift = 0:K:N-1;
            pha = (1j).^(0:M-1);
            G = cell2mat(arrayfun(@(s_, pha_m) circshift(g, s_)*pha_m, shift, pha, 'UniformOutput', false));
            freq_i = 0: 1/K :1-1/K;
            freq_q = freq_i + 1/K/2;
            Ei = num2cell(exp(1j*2*pi*(0:N-1)'*freq_i), 1);
            Eq = num2cell(exp(1j*2*pi*(0:N-1)'*freq_q), 1);
            Ai = cell2mat(cellfun(@(E_k) diag(E_k)*G, Ei, 'UniformOutput', false));
            Aq = cell2mat(cellfun(@(E_k) diag(E_k)*G*(1j), Eq, 'UniformOutput', false));
        otherwise
            error('Non-specific method');
    end
    
    if strcmpi(group_order, 'ss')
        idx = 1:N;
        idx = reshape(idx, M, K).';
        idx_ss = idx(:);
        [Ai, Aq] = deal(Ai(:,idx_ss), Aq(:,idx_ss));
    end

end