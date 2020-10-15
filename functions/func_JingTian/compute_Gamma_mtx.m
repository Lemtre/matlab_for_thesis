function Gamma_p = compute_Gamma_mtx(alpha_p, ofdm, D)
% COMPUTE_GAMMA_MTX compute the Gamma_p matrix for a given Doppler, alpha_p.
% 
%   alpha_p: relative Doppler factor.
%   ofdm:    structure containing OFDM system parameters.
%   D:       the depth of ICI.
%   Gamma_p: the output Gamma_p matrix.
% 
% Example:
%   ofdm.Nc = 1024;
%   ofdm.ovrsamp = 2;
%   ofdm.window = rcos_window(ofdm.Nc,1/4);
%   ofdm.bw = 6e3;
%   ofdm.fc = 12e3;
%   alpha = 7e-5;
%   D = 5;
%   Gamma = compute_Gamma_mtx(alpha, ofdm, D);
%   figure; imagesc(abs(Gamma));
% 
% Authour: Jing Tian
% Last modified: 10th Sept, 2020.

M = ofdm.ovrsamp;
K = ofdm.Nc;
g = ofdm.window(:);
fc = ofdm.fc;
bw = ofdm.bw;
if mod(K,2)~=0
    error('The number of subcarriers is supposed to be even.');
end

L_czt_out = 2*M*D+1;
W = exp(-1j*2*pi*(1+alpha_p)/M/K);
kk = -K/2 : K/2-1;
Ak = exp(1j*2*pi*(fc/bw*alpha_p+(kk*alpha_p-D*(1+alpha_p))/K));
gamma_k = arrayfun(@(ai) segcztResample(g, L_czt_out, W, ai), Ak,...
    'UniformOutput', false);
gamma_k = num2cell([cell2mat(gamma_k); zeros(M*K-L_czt_out,K)], 1);
shift = num2cell(M*(0:K-1)-D*M, 1);
gamma_k_shift = cellfun(@(gi,si) circshift(gi,si), gamma_k, shift,...
    'UniformOutput', false);
Gamma_p = sparse(cell2mat(gamma_k_shift))/K;

end

