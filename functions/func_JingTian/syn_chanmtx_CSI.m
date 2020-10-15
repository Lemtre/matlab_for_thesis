function H = syn_chanmtx_CSI(chan, ofdm, D)
% SYN_CHANNELMTX synthesize the channel matrix with full CSI.
% 
%   chan: structure containing channel state information (CSI).
%   ofdm: structure containing OFDM system parameters.
%   H:    frequency-domain channel matrix.
% 
% Example:
%   chan.PathDelays = [0 0.001 0.0027 0.0042].';
%   chan.PathGains = [0.2354 0.6331 0.7370 0.7079].';
%   chan.DopplerFactor = 1e-4*[0.3960 -1.4573 -0.8847 -0.9607].';
%   ofdm.Nc = 1024;
%   ofdm.ovrsamp = 2;
%   ofdm.window = rcos_window(ofdm.Nc,1/4);
%   ofdm.bw = 6e3;
%   ofdm.fc = 12e3;
%   D = 5;
%   H = syn_channelmtx(chan, ofdm, D);
%   figure; imagesc(abs(H));
% 
% Authour: Jing Tian
% Last modified: 10th Sept, 2020.

amp = chan.PathGains(:).';
tau = chan.PathDelays(:).';
alpha = chan.DopplerFactor(:).';
xi = num2cell(arrayfun(@(a,t,d)(1+d)*a*exp(-1j*2*pi*ofdm.fc*t),...
    amp, tau, alpha), 1);
Lambda = arrayfun(@(t) compute_Lambda_mtx(t,ofdm), tau,...
    'UniformOutput', false);
Gamma = arrayfun(@(d) compute_Gamma_mtx(d,ofdm,D), alpha,...
    'UniformOutput', false);
Hp = cellfun(@(a,b,c) a*b*c, xi, Lambda, Gamma, 'UniformOutput', false);
size_H = size(Hp{1});
Hp = reshape(full(cell2mat(Hp)), [size_H numel(amp)]);
H = sparse(sum(Hp, 3));

end