function g = gen_Pei_HGLvec0(N, m)
% GEN_Pei_HGLVEC0 generate 0th-order HGL DFT eigenvector using the
%   closed-form expression in Pei Soo-Chang's paper,?DOI 10.1109/TSP.2016.2540601.
% 
%   N: dimension of the DFT matrix.
%   m: bandwidth.
%   g: 0th-order discrete HGL vector.
% 
% Author: Jing Tian
% Last modified: Dec 11th, 2019.
L = floor(N/2);
if mod(L,2)~=0
    warning('The length constraint is violated: FLOOR(N/2) should be even.');
end
nn = (0:N-1)';
ss = (m+1):L;
g = prod(bsxfun(@minus, cos(2*pi*nn/N), cos(2*pi*ss/N)), 2);
g = g./norm(g);

end