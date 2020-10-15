function x_gfdm = make_GFDMsym_AD(gfdmPara, data, pilots)
%  MAKE_GFDMSYM_AD make GFDM symbol with auxiliary-data(AD)-aided 
%                  interference-free pilots.  
% 
%   gfdmPara: GFDM parameters.
%   data:     informative data.
%   pilots:   pilots symbols.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

K = gfdmPara.Nsc; % number of subcarrier
M = gfdmPara.Nss; % number of subsymbols
g = gfdmPara.filter; % prototype filter
grp = gfdmPara.group_order; % whether to group subcarriers or subsymbols together
coeMtx = gfdmPara.ADcoeMtx; % AD coefficient matrix
id_pilot = gfdmPara.id_pilot; % pilot indices
id_AD = gfdmPara.id_AD;  % AD indices 
id_data = gfdmPara.id_data; % free data indicies

N = K*M;
data = data(:);

d0 = zeros(N,1);
d0(id_data) = data;
df = mod_FreqDom_GFDM(g, K, M, grp, d0);
delta = pilots - df(id_pilot);
dp = coeMtx * delta; 

d0(id_AD) = dp;
x_gfdm = fast_gfdm_mod(g, K, M, grp, d0);

end