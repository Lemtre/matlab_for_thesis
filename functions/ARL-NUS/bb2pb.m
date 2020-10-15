% BB2PB Convert complex baseband signal to passband
%
% y = bb2pb(x,Fd,Fc,Fs)
%   x is the complex baseband signal (column vector)
%   Fd is the baseband sampling rate
%   Fc is the carrier frequency
%   Fs is the passband sampling rate
%   y is the passband signal
%
% Author: Mandar Chitre
% Last modified: Dec 8, 2009

function y = bb2pb(x,Fd,Fc,Fs)

x = resample(x,Fs,Fd);
t = (0:length(x)-1)'/Fs;
carrier = exp(2i*pi*Fc*t);
y = real(x.*carrier);
