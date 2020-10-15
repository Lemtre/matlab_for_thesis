function Lambda_p = compute_Lambda_mtx(tau_p, ofdm)
% COMPUTE_LAMBDA_MTX compute the Lambda_p matrix for a given delay, tau_p.
% 
%   tau_p:   actual delay (in seconds).
%   ofdm:    structure containing OFDM system parameters.
%   Gamma_p: the output Gamma_p matrix.
% 
% Example:
%   ofdm.Nc = 1024;
%   ofdm.ovrsamp = 2;
%   ofdm.bw = 6e3;
%   tau = 7e-3;
%   Lambda = compute_Lambda_mtx(tau, ofdm);
%   figure; imagesc(abs(Lambda));
% 
% Authour: Jing Tian
% Last modified: 10th Sept, 2020.

M = ofdm.ovrsamp;
K = ofdm.Nc;
T = K/ofdm.bw;
if mod(K,2)~=0
    error('The number of subcarriers is supposed to be even.');
end
mm = -M*K/2 : M*K/2-1;
Lambda_p = sparse(diag(exp(-1j*2*pi*mm*tau_p/M/T)));

end