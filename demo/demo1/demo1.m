clear
close all
dbstop if error
%% set parameters
set_para_BP; % set parameters
%% generate channel
ch = gen_Channel(ofdm, chanPara);
% figure; stem(ch.PathDelays, ch.PathGains); xlim([0 ofdm.Tg]);
%% BP-OFDM-Tx
% -------------------------------------------------------------------------
tx_bit = gen_binbits(modType, Nd, 1); % generate 0-1 bits
% -------------------------------------------------------------------------
[tx_sym, ~] = bbSymMap(tx_bit, modType); % constellation mapping
% -------------------------------------------------------------------------
pilots = exp(1j*2*pi*rand(Np,1)); % generate pilot symbols 
% pilots = bbSymMap(randi([0 1], Np, 2), modType);
% -------------------------------------------------------------------------
xbb = make_OFDMsym(ofdm, tx_sym, pilots, 'shiftOp',1); % make OFDM symbol
% -------------------------------------------------------------------------
xbb_g = add_GuardInterval(xbb, ofdm.Ng, 'zp'); % add ZP
% -------------------------------------------------------------------------
xpb = bb2pb(xbb_g, ofdm.bw, ofdm.fc, fs_tx); % convert to passband
% -------------------------------------------------------------------------
ypb = add_channel(xpb, fs_tx, ch); % pass through channel
ypb = awgn(ypb, snr, 'measured');
figure; plot(ypb); xlim([-100 length(ypb)+100]); title('Received signal');
% -------------------------------------------------------------------------
%% OFDM-Rx
% -------------------------------------------------------------------------
zbb = pb2bb(ypb, ofdm.bw, ofdm.fc, fs_tx); % convert to baseband
% -------------------------------------------------------------------------
% transform into frequency domain
z_1blk = zbb(1 : length(ofdm.window)+ofdm.Ng); % disgard samples outside ZP
if ofdm.ovrsamp == 1 % overlap-add 
    Zf = ser2para_rmvGI_fft(z_1blk, ofdm.Nc, 1, 'zp', 'shiftOp', 1);
else % or oversample
    Zf = fftshift(fft(z_1blk, ofdm.ovrsamp*ofdm.Nc), 1);
end
Zf_0ovrsamp = Zf(1: ofdm.ovrsamp : end);
% -------------------------------------------------------------------------
% estimate noise variance via null carriers
Z0 = Zf(idz_ovrsamp);
pn = (Z0'*Z0)/length(idz_ovrsamp);
%% MMSE equalization with known CSI
% -------------------------------------------------------------------------
% reconstruct channel matrix using known CSI
H_CSI = syn_chanmtx_CSI(ch, ofdm, D); % reconstruct channel matrix H using CSI
% -------------------------------------------------------------------------
% MMSE equalization using known CSI
s_CSI = mmse_eqlz(Zf, H_CSI, pn);
% -------------------------------------------------------------------------
% demodulation
rx_sym_CSI = s_CSI(ofdm.idd);
[~, rx_bit_CSI] = bbSymDemap(modType, rx_sym_CSI);
[Nerr_CSI, ber_CSI] = biterr(tx_bit, rx_bit_CSI);

%% MMSE equalization with BP channel estimation
% -------------------------------------------------------------------------
% compute channel matrices
[Hp,delay_set,Doppler_set] = compute_chanmtxset(ofdm, 2, maxDopplerFactor, D);
% -------------------------------------------------------------------------
% compute dictionary matrix
A  = compute_dicmtx(ofdm, pilots, false, Hp);
% -------------------------------------------------------------------------
% BP 2-D channel estimation 
mu = 50; % l1-regularizer coefficient
obz = Zf(ofdm.idobz); % observations
fprintf(1, '-------------------------------------------------------------\n');
fprintf(1, 'BP optimization begins...\n');
[xi,~,~,timeBPopt] = SpaRSA(obz, A, mu, 'Verbose', false);
fprintf(1, 'Total elapsed time: %0.1f seconds\n', timeBPopt(end));
fprintf(1, 'BP optimization complete!\n');
fprintf(1, '-------------------------------------------------------------\n');
% -------------------------------------------------------------------------
% channel matrix reconstruction using BP estimation result
[H_BP, id_xi] = syn_chanmtx_BP(xi, Hp);
% -------------------------------------------------------------------------
% MMSE equalization using BP result
s_BP = mmse_eqlz(Zf, H_BP, pn);
% -------------------------------------------------------------------------
% demodulation
rx_sym_BP = s_BP(ofdm.idd);
[~, rx_bit_BP] = bbSymDemap(modType, rx_sym_BP);
[Nerr_BP, ber_BP] = biterr(tx_bit, rx_bit_BP);
%% show results
show_results_BP;