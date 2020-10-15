function y = per_interp_delay(x, tau)
% PER_INTERP_DELAY interpolation of periodic signals at arbitrary delays.
%   x:   uniform-spaced sampled discrete sequence of a periodic signal with
%        period T.
%   tau: delays normalized with period T.
%   y:   uniform-spaced sampled discrete sequence
% 
% Author: Jing Tian
% Last modified: 9th Oct, 2019.
tau = tau(:);
xf = fft(x(:));
N = length(x);
if mod(N,2)==0
    Gamma_sigma = fftshift([(exp(1j*2*pi*N/2*tau)+exp(-1j*2*pi*N/2*tau))/2,...
        exp(-1j*2*pi*tau*(-N/2+1:N/2-1))], 2);
else
    Gamma_sigma = ifftshift(exp(-1j*2*pi*tau*(-(N-1)/2:(N-1)/2)), 2);
end
xfs = diag(xf) * Gamma_sigma.';
y = ifft(xfs);

if isreal(x)
    y = real(y);
end
% if size(x,1) < size(x,2)
%     y = y.';
% end

end