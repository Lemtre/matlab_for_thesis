function [Gi, Gq] = gen_freqDom_GFDM_OQAMmtx(g, K, M, group_order, oqam_method)
% GEN_FREQDOM_GFDM_OQAMMTX generate the frequency domain GFDM/OQAM modulation
%   matrix for a given prototype filter "g";
% 
%   g:           the prototype filter in time domain;
%   K:           number of subcarriers;
%   M:           number of subsymbols;
%   group_order: group data into subcarriers ('sc') or subsymbols ('ss')
%   oqam_method: the realization of OQAM in time domain, i.e., 'SMT' or 'CMT'.
%   Gi:          the in-phase frequency domain modulation matrix.
%   Gq:          the quadrature frequency domain modulation matrix.
% 
% Author: Jing Tian
% Last modified: June 26th, 2020.
    if nargin < 5
        oqam_method = 'SMT';
    end
    
    g = g(:);
    N = length(g);
    if N~=K*M
        error('Filter length mismatch!');
    elseif mod(K,2)~=0 || mod(M,2)~=0
        error('The numbers of subsymbols and subcarriers must be even!');
    end
    gf = fft(g);
    gf = gf./norm(gf);
    
    mm = 0:M-1;
    kk = 0:K-1;
    nn = 0:N-1;
    eqvfreq = exp(-1j*2*pi*nn.'*mm/M);

    time_shift = @(s_) circshift(gf, s_);
    shift = 0:M:N-1;
    G = arrayfun(time_shift, shift, 'UniformOutput', false);
   
    switch lower(oqam_method)
        case 'smt'
            eqvfreq_c = exp(-1j*2*pi*nn.'*(mm+1/2)/M);
            pha_i = num2cell(1j.^kk);
            pha_q = num2cell(1j.^(3*kk+1));
            Gi = cell2mat(cellfun(@(gf_k, pha_k) diag(gf_k)*eqvfreq*pha_k,...
                G, pha_i, 'UniformOutput', false));
            Gq = cell2mat(cellfun(@(gf_k, pha_k) diag(gf_k)*eqvfreq_c*pha_k,...
                G, pha_q, 'UniformOutput', false));
        case 'cmt'
            G_c = arrayfun(time_shift, shift+M/2, 'UniformOutput', false);
            pha_i = (1j).^mm;
            pha_q = (1j).^(3*mm+1);
            freq_pha_i = eqvfreq .* repmat(pha_i, N, 1);
            freq_pha_q = eqvfreq .* repmat(pha_q, N, 1);
            Gi = cell2mat(cellfun(@(gf_k) diag(gf_k)*freq_pha_i, G, 'UniformOutput', false));
            Gq = cell2mat(cellfun(@(gf_k) diag(gf_k)*freq_pha_q, G_c, 'UniformOutput', false));
        otherwise
            error('Non-specific method');
    end
    
    if strcmpi(group_order, 'ss')
        idx = 1:N;
        idx = reshape(idx, M, K).';
        idx_ss = idx(:);
        [Gi, Gq] = deal(Gi(:,idx_ss), Gq(:,idx_ss));
    end

end