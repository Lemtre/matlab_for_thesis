function A = paf_eqv(x, y, sigma, nu)
% PAF_EQV calculate the doubly periodic (cross) ambiguity function (PAF)
%         between two periodic discrete signals.
% 
%   x, y:  periodic discrete signals of the same length;
%   sigma: range of normalized delay;
%   nu:    range of normalized frequency;
%   A:     output PAF in the desired range.
%  
% Author: Jing Tian
% Last modified: 7th Oct, 2019.

x = x(:);
y = y(:);
nu = nu(:);
sigma = sigma(:);

Nx = length(x);
Ny = length(y);
if Nx ~= Ny
    error('The two input sequences should have the same length.');
end
Ys = per_interp_delay(y, sigma);

Gamma_nu = exp(1j*2*pi*nu*(0:Nx-1));

x_Gamma_nu = Gamma_nu * diag(x');
A = x_Gamma_nu * Ys;

end