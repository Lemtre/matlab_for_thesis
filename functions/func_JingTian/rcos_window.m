function [rcoswin, beta_out] = rcos_window(T, beta_in)
% RCOS_WINDOW generate raised cosine window
% 
%   T:        length of original window length.
%   beta_in:  the desired roll-off factor.
%   rcoswin:  the raised cosine window.
%   beta_out: the possibly modified roll-off factor.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.
    
    Lp = T * beta_in;
    L = Lp + T;
    beta_out = beta_in;
    if Lp ~= fix(Lp)
        warning('\beta has been adjusted to make the length of raised cosine window integer');
        Lp = ceil(Lp);
        L = Lp + T;
        beta_out = Lp/T;
    end
    rcoswin = zeros(L,1);
    rcoswin(Lp+1:end-Lp) = 1;
    t1 = 1:Lp;
    t2 = T+1:L;
    rcoswin(t1) = 0.5*(1+cos(pi/(beta_out*T)*(abs((t1)-(1+beta_out)/2*T)-(1-beta_out)/2*T)));
    rcoswin(t2) = 0.5*(1+cos(pi/(beta_out*T)*(abs((t2)-(1+beta_out)/2*T)-(1-beta_out)/2*T)));
    
end