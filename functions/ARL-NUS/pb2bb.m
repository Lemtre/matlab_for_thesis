% PB2BB Convert passband signal to complex baseband
%
% x = pb2bb(y,Fd,Fc,Fs)
%   y is the passband signal (column vector)
%   Fd is the baseband sampling rate
%   Fc is the carrier frequency
%   Fs is the passband sampling rate
%   x is the complex baseband signal
%
% Author: Mandar Chitre
% Last modified: Dec 8, 2009

function x = pb2bb(y,Fd,Fc,Fs)

t = (0:length(y)-1)'/Fs;
osc = exp(2i*pi*Fc*t);
x = 2*conj(osc).*y;
Hb = fir1(256,2.1*Fd/Fs);
x = filtfilt(Hb,1,x);
x = resample(x,2*Fd,Fs);
x = x(1:2:end);
