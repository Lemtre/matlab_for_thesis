function hn_eval = hermite_func(n, t_var)
% HERMITE_FUNC generate Hermite functions of arbitrary order.
% 
%   n:       order.
%   t_var:   time.
%   hn_eval: order-n Hermite function within the range of t_var.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

syms t;
x = exp(-t^2);
y = @(nn)diff(x, nn);
Hn_poly = @(nn)(-1)^nn*exp(t^2)*y(nn);
hn = @(nn,tt)subs(Hn_poly(nn), tt).*exp(-tt.^2/2)/sqrt((2^nn*factorial(nn)*sqrt(pi)));
hn_eval = hn(n, t_var);
hn_eval = double(hn_eval(:));

end
