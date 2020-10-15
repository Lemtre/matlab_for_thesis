clear
close all
dbstop if error
%% set parameters
set_para_OFDMchEst; % set parameters
%% generate channel
ch = gen_Channel(ofdm, chanPara);
% figure; stem(ch.PathDelays, ch.PathGains); xlim([0 ofdm.Tg]);
%% OFDM-Tx
% -------------------------------------------------------------------------
tx_bit = gen_binbits(modType, Nd, 1); % generate 0-1 bits
% -------------------------------------------------------------------------
[tx_sym, ~] = bbSymMap(tx_bit, modType); % constellation mapping
% -------------------------------------------------------------------------
pilots = exp(1j*2*pi*rand(Np,1)); % generate pilot symbols 
% -------------------------------------------------------------------------
xbb = make_OFDMsym(ofdm, tx_sym, pilots, 'shiftOp',1); % make OFDM symbol
% -------------------------------------------------------------------------
xbb_g = add_GuardInterval(xbb, ofdm.Ng, typeGI); % add ZP
% -------------------------------------------------------------------------
xpb = bb2pb(xbb_g, ofdm.bw, ofdm.fc, fs_tx); % convert to passband
% -------------------------------------------------------------------------
ypb = clipping(xpb, 6); % clip to reduce PAPR
% -------------------------------------------------------------------------
ypb = add_channel(ypb, fs_tx, ch); % pass through channel
ypb = awgn(ypb, snr, 'measured');
figure; plot(ypb); xlim([-100 length(ypb)+100]); title('Received signal');
%% OFDM-Rx
% -------------------------------------------------------------------------
zbb = pb2bb(ypb, ofdm.bw, ofdm.fc, fs_tx); % convert to baseband
% -------------------------------------------------------------------------
% transform into frequency domain
z_1blk = zbb(1 : length(ofdm.window)+ofdm.Ng); % disgard samples outside ZP
Zf = ser2para_rmvGI_fft(z_1blk, ofdm.Nc, 1, typeGI, 'shiftOp', 1);
%% Noise power estimation
Z0 = Zf(ofdm.idz);
pn = Z0(:)'*Z0(:)/Nz;
%% Channel estimation
% -------------------------------------------------------------------------
% LS channel estimation at pilots
h_LS = chEst_ofdm_LS(Zf, pilots, ofdm.idp);
% -------------------------------------------------------------------------
% Conventional 1-D interpolation
cfr_interp = ofdm_ch_interp1(h_LS, ofdm.idp, ofdm.Nc, 'pchip');
% -------------------------------------------------------------------------
% Transform domain method: IDFT
[cfr_TD, ~] = chEst_TransDom_LP_IDFT(h_LS, ofdm.idp, ofdm.Nc, ofdm.Ng, []);
% -------------------------------------------------------------------------
%% Equalization and demodulation
% -------------------------------------------------------------------------
% Conventional 1-D interpolation + zero-forcing (ZF) equalization
s1 = Zf./cfr_interp;
rx_sym1 = s1(ofdm.idd);
% scatterplot(rx_sym1); title('Conventional 1-D interpolation');
[~, rx_bit1] = bbSymDemap(modType, rx_sym1);
[~, ber1] = biterr(tx_bit, rx_bit1);
% -------------------------------------------------------------------------
% Transform domain method + zero-forcing (ZF) equalization
s2 = Zf./cfr_TD;
rx_sym2 = s2(ofdm.idd);
% scatterplot(rx_sym2); title('Transform domain method: IDFT');
[~, rx_bit2] = bbSymDemap(modType, rx_sym2);
[~, ber2] = biterr(tx_bit, rx_bit2);
% -------------------------------------------------------------------------
% Conventional 1-D interpolation + MMSE equalization
s3 = mmse_eqlz(Zf, diag(cfr_interp), pn);
rx_sym3 = s3(ofdm.idd);
% scatterplot(rx_sym3); title('Transform domain method: IDFT');
[~, rx_bit3] = bbSymDemap(modType, rx_sym3);
[~, ber3] = biterr(tx_bit, rx_bit3);
% -------------------------------------------------------------------------
% Transform domain method + MMSE equalization
s4 = mmse_eqlz(Zf, diag(cfr_TD), pn);
rx_sym4 = s4(ofdm.idd);
% scatterplot(rx_sym4); title('Transform domain method: IDFT');
[~, rx_bit4] = bbSymDemap(modType, rx_sym4);
[~, ber4] = biterr(tx_bit, rx_bit4);
%% Show results
show_results_OFDMchEst;
