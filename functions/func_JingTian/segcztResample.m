function z = segcztResample(x, M, w0, a0)
% SEGCZTRESAMPLE resample in frequency domain using chirp-Z transform (CZT)
% Perform segment-wise CZT if the input sequence is too long.
% 
%   x:  input signal.
%   a0: starting point on the Z plane.
%   w0: spiral radius and angle.
%   M:  length of the output sequence.
% 
% Author: Jing Tian
% Last modified: April 26th, 2019

x = x(:);
N = length(x);
f = N./factor(N);
[~,id] = min(abs(f-M));
% L = round(N/32);
L = f(id);
if mod(N,L) ~= 0
    error('The input sequence cannot be divided into integer segments of length L!');
end
S = round(N/L);
segX = reshape(x, L, []);
segcztX = czt(segX, M, w0, a0);
mm = 0:M-1;
ss = 0:S-1;
Q = mm.'*ss*L;
Ws = w0.^Q;
As = repmat(a0.^(-ss*L), M, 1);
z = sum(segcztX .* Ws .* As, 2);
z = z(:);

end