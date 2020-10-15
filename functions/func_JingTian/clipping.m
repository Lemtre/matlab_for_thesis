function sig = clipping(sig, ratio)
% CLIPPING clip the signal using a simple model called soft limiter. 
% 
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.
sig = sig(:);
L = length(sig);
P = sig'*sig/L;
A = ratio*sqrt(P);
sig(sig>A) = A;
sig(sig<-A) = -A;

end