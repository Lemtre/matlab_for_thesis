function [coeMtx, id_AD, id_data] = gen_ADcoeMtx_OQAM...
    (g, K, M, group_order, oqam_method, id_pilot, posAD)
% GEN_ADCOEMTX_OQAM generate the coefficient matrix for auxiliary data (AD) 
%                   aided CFBMC/OQAM.
% 
%   g:           prototype filter.
%   K:           number of subcarriers.
%   M:           number of subsymbols.
%   group_order: data group in subcarriers ('sc') or subsymbols ('ss').
%   oqam_method: OQAM method, 'smt' or 'cmt'.
%   id_pilot:    indices of pilots (which rows in the frequency domain 
%                modulation matrix, 1-based).
%   posAD:       a cell array containing AD positions. Each cell is a two-column
%                matrix whose first column corresponds to subcarrier indices
%                and second column subsymbol indices (0-based).
%   coeMtx:      output coefficient matrix.
%   id_AD:       indices of AD (which columns in the time/frequency domain
%                modulation matrix, 1-based).
%   id_data:     indices of free data (which columns in the time/frequency 
%                domain modulation matrix, 1-based).
% 
% Author: Jing Tian
% Last modified: 9th April, 2020. 

if ~iscell(posAD) || length(posAD)~=2
    error('The #7 input should be a 2 by 1 cell array containing AD positions.');
elseif ~strcmpi(oqam_method,'smt') && ~strcmpi(oqam_method,'cmt')
    error([oqam_method ' is an invalid OQAM method.']);
end

switch lower(group_order)
    case 'ss'
        get_ADid = @(posADmtx) posADmtx(:,2)*K + posADmtx(:,1); 
    case 'sc'
        get_ADid = @(posADmtx) posADmtx(:,1)*M + posADmtx(:,2); 
    otherwise
        error('Invalid grouping.');
end
N = K*M;
id_AD_i = get_ADid(posAD{1});
id_AD_q = get_ADid(posAD{2});
id_AD = union(id_AD_i, id_AD_q+N)+1;
id_data = setdiff((1:2*N)', id_AD);
if mod(length(id_data), 2)~=0
    warning('The number of free, real data is supposed to be even. The last data position has been discarded.');
    id_data = id_data(1:end-1);
end

[Gi, Gq] = gen_freqDom_GFDM_OQAMmtx(g, K, M, group_order, oqam_method);
G = [Gi Gq];
G_real = [real(G); imag(G)];
id_row = union(id_pilot, id_pilot+N);
Phi = G_real(id_row, id_AD);
coeMtx = Phi'/(Phi*Phi');

end