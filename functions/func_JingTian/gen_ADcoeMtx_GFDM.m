function [coeMtx, id_AD, id_data] = gen_ADcoeMtx_GFDM...
    (g, K, M, group_order, id_pilot, posAD)
% GEN_ADCOEMTX_GFDM generate the coefficient matrix for auxiliary data (AD) 
%                   aided GFDM.
% 
%   g:           prototype filter.
%   K:           number of subcarriers.
%   M:           number of subsymbols.
%   group_order: data grouped in subcarriers ('sc') or subsymbols ('ss'). 
%   id_pilot:    indices of pilots (which rows in the frequency domain 
%                modulation matrix, 1-based).
%   posAD:       AD positions (0-based), two-column matrix.
%   coeMtx:      output coefficient matrix.
%   id_AD:       indices of AD (which columns in the time/frequency domain
%                modulation matrix, 1-based).
%   id_data:     indices of free data (which columns in the time/frequency 
%                domain modulation matrix, 1-based).
% 
% Author: Jing Tian
% Last modified: 9th April, 2020. 

if ~ismatrix(posAD) || size(posAD, 2)~=2
    error('The #6 input should be a two-column matrix containing AD positions.');
end

switch lower(group_order)
    case 'ss'
        id_AD = posAD(:,2)*K + posAD(:,1); 
    case 'sc'
        id_AD = posAD(:,1)*M + posAD(:,2);        
    otherwise
    error('Invalid grouping.');
end
id_AD = sort(id_AD(:))+1;
G = gen_freqDom_GFDMmtx(g, K, M, group_order);
Phi = G(id_pilot, id_AD);
coeMtx = Phi'/(Phi*Phi');
id_data = setdiff((1:K*M).', id_AD);

end