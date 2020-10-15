function x_gfdm_oqam = make_GFDM_OQAMsym_AD(oqamPara, data, pilots)
%  MAKE_GFDMSYM_AD make GFDM/OQAM symbol with auxiliary-data(AD)-aided 
%                  interference-free pilots.  
% 
%   oqamPara: GFDM/OQAM parameters.
%   data:     informative data.
%   pilots:   pilots symbols.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

K = oqamPara.Nsc; % number of subcarrier
M = oqamPara.Nss; % number of subsymbols
g = oqamPara.filter; % prototype filter
grp = oqamPara.group_order; % whether to group subcarriers or subsymbols together
method = oqamPara.oqam_method;  % 'smt' or 'cmt'
coeMtx = oqamPara.ADcoeMtx; % AD coefficient matrix
id_pilot = oqamPara.id_pilot; % pilot indices
id_AD = oqamPara.id_AD;  % AD indices 
id_data = oqamPara.id_data; % free data indicies

N = K*M;
data = data(:);

d0 = zeros(2*N, 1);
d0(id_data) = [real(data); imag(data)];
D0 = reshape(d0, [], 2);
D0_cplx = complex(D0(:,1), D0(:,2));
df = mod_FreqDom_GFDM_OQAM(g, K, M, grp, method, D0_cplx);
delta = pilots - df(id_pilot);
dp = coeMtx * [real(delta); imag(delta)];
d0(id_AD) = dp;
D_cplx = complex(d0(1:N), d0(N+1:end));
x_gfdm_oqam = mod_gfdm_oqam(g, K, M, grp, method, D_cplx);

end