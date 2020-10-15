function [Hp, delay_set, Doppler_set] = compute_chanmtxset(ofdm,...
    ovrsamp_in_delay, maxDopplerFactor, D)
% COMPUTE_CHANMTXSET compute channel matrices set for the BP method.
% 
%   ofdm:             structure containing OFDM parameters.
%   ovrsamp_in_delay: over-sampling rate along the delay axis.
%   maxDopplerFactor: maximum relative Doppler factor.
%   D:                ICI depth.     
%   Hp:               cell array of channel matrices.
%   delay_set:        set of delay grids.
%   Doppler_set:      set of Doppler grids.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

fs = ovrsamp_in_delay * ofdm.bw;
N_tau = floor(ofdm.Tg*fs)+1;
delay_set = linspace(0, ofdm.Tg, N_tau).';

delta_f = 0.05*ofdm.bw/ofdm.Nc; % Doppler searching step: 0.05 of subcarrier spacing.
N_alpha = 2*ceil(maxDopplerFactor*ofdm.fc/delta_f)+1;
Doppler_set = linspace(-maxDopplerFactor, maxDopplerFactor, N_alpha).';

ii = 1:length(Doppler_set);
jj = 1:length(delay_set);
[X,Y] = meshgrid(ii, jj);
id_2D = [X(:) Y(:)];

fprintf(1, '-------------------------------------------------------------\n');
t0 = tic;
fprintf(1, 'Computing Gamma matrices...\n');
tic;
Gamma_set = arrayfun(@(d) compute_Gamma_mtx(d, ofdm, D), Doppler_set,...
    'UniformOutput', false);
toc;
fprintf(1, 'Computing Lambda matrices...\n');
tic;
Lambda_set = arrayfun(@(d) compute_Lambda_mtx(d, ofdm), delay_set,...
    'UniformOutput', false);
toc;

fprintf(1, 'Computing channel matrices at interwoven grids...\n');
tic;
Hp = cellfun(@(q) Lambda_set{q(2)}*Gamma_set{q(1)}, num2cell(id_2D,2),...
    'UniformOutput', false);
toc;
t1 = toc(t0);
fprintf(1, 'Channel matrices computation done!\n');
fprintf(1, 'Total elapsed time: %0.1f seconds\n', t1);
fprintf(1, '-------------------------------------------------------------\n');


end