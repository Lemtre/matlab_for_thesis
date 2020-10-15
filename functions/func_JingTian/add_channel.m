function rx = add_channel(tx, fs, chan)
% ADD_CHANNEL add channel effects to a signal vector
% 
%   tx:   transmitted signal.
%   fs:   sampling frequency.
%   chan: channel delay-Doppler-amplitude triplets.
%   rx:   received signal passed through channel
% 
% Example: 
%   chan.PathDelays = [0 0.0001 0.0002 0.0004].';
%   chan.PathGains = [0.2354 0.6331 0.7370 0.7079].';
%   chan.DopplerFactor = 1e-4*[0.3960 -1.4573 -0.8847 -0.9607].';
%   fs = 96e3; t = 0:1/fs:1; x = sin(2*pi*6e3*t);
%   y = add_channel(x, fs, chan);
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

tx = tx(:);
L_tx = length(tx);
L_max = ceil(L_tx * (1+max(chan.DopplerFactor)));

add_dpl = @(dpl)interp1(0:L_tx-1, tx, (0:L_tx*(1+dpl)).'/(1+dpl), 'PCHIP', 'extrap');
tx_scaled = arrayfun(@(a,d) a*add_dpl(d), chan.PathGains, chan.DopplerFactor,...
    'Uniformoutput', false);

delay_in_samples = round(chan.PathDelays * fs);
L_delay = delay_in_samples(end);
delay_in_samples = num2cell(delay_in_samples(:));
L = L_max + L_delay;
tx_multipath = cellfun(@(x,tau) circshift([x; zeros(L-length(x),1)], tau),...
    tx_scaled, delay_in_samples, 'Uniformoutput', false);
tx_multipath = tx_multipath(:).';
rx = sum(cell2mat(tx_multipath), 2);

end