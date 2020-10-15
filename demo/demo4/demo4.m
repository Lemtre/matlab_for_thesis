clear
close all
dbstop if error
%% set parameters
set_para_CFBMC_AD; % set parameters
%% generate channel
ch = gen_Channel(sysPara, chanPara);
% figure; stem(ch.PathDelays, ch.PathGains); xlim([0 ofdm.Tg]);
%% generate baseband symbols for GFDM
% -------------------------------------------------------------------------
% compute AD coefficient matrix and the positions of AD and data.
[gfdmPara.ADcoeMtx, gfdmPara.id_AD, gfdmPara.id_data] = ...
    gen_ADcoeMtx_GFDM(g_gfdm, K, M, group_order, idp, posAD_gfdm);
Nd_gfdm = length(gfdmPara.id_data);
eta_gfdm = Nd_gfdm/(N+Ng); % spectral efficiency
% -------------------------------------------------------------------------
tx_bit_gfdm = gen_binbits(modType, Nd_gfdm, 1); % generate 0-1 bits
% -------------------------------------------------------------------------
[tx_sym_gfdm, ~] = bbSymMap(tx_bit_gfdm, modType); % constellation mapping
% -------------------------------------------------------------------------
pilots = exp(1j*2*pi*rand(Np,1)); % generate pilot symbols 
%% GFDM: Tx
% -------------------------------------------------------------------------
xbb_gfdm = make_GFDMsym_AD(gfdmPara, tx_sym_gfdm, pilots); % make GFDM symbol
% -------------------------------------------------------------------------
xbb_gfdm_g = add_GuardInterval(xbb_gfdm, Ng, typeGI); % add CP
% -------------------------------------------------------------------------
xpb_gfdm = bb2pb(xbb_gfdm_g, bw, fc, fs); % convert to passband
% figure; plot(xpb_gfdm); xlim([-100 length(xpb_gfdm)+100]); title('GFDM signal');
% -------------------------------------------------------------------------
ypb_gfdm = add_channel(xpb_gfdm, fs, ch); % pass through channel
ypb_gfdm = awgn(ypb_gfdm, snr, 'measured');
%% GFDM: Rx
% -------------------------------------------------------------------------
zbb_gfdm = pb2bb(ypb_gfdm, bw, fc, fs); % convert to baseband
% -------------------------------------------------------------------------
% transform into frequeNy domain
z_1blk_gfdm = zbb_gfdm(1 : N+Ng); % disgard samples outside guard interval
Zf_gfdm = ser2para_rmvGI_fft(z_1blk_gfdm, N, 1, typeGI, 'shiftOp', false);
%% GFDM: channel estimation
% -------------------------------------------------------------------------
% LS channel estimation at pilots
h_LS_gfdm = chEst_ofdm_LS(Zf_gfdm, pilots, idp);
% -------------------------------------------------------------------------
% Transform domain channel estimation
[cfr_gfdm, ~] = chEst_TransDom_LP_IDFT(h_LS_gfdm, idp, N, Ng, []);
% -------------------------------------------------------------------------
%% GFDM: equalization
% -------------------------------------------------------------------------
% Zero-forcing (ZF) channel equalization by default in this section.
% Beware of the difference between channel effect compensation and symbol
% demodulation, which can also be named "ZF".
s_gfdm = Zf_gfdm./cfr_gfdm;
% -------------------------------------------------------------------------
% matched filtering (MF) equalization/demodulation 
mf_gfdm = mf_FreqDom_GFDM(g_gfdm, K, M, group_order, s_gfdm);
rx_sym_gfdm_mf = mf_gfdm(gfdmPara.id_data);
% scatterplot(rx_sym_gfdm_mf); title('MF GFDM');
[~, rx_bit_gfdm_mf] = bbSymDemap(modType, rx_sym_gfdm_mf);
[~, ber_gfdm_mf] = biterr(tx_bit_gfdm, rx_bit_gfdm_mf);
% -------------------------------------------------------------------------
% zero-forcing (ZF) equalization/demodulation
zf_gfdm = zf_FreqDom_GFDM(g_gfdm, K, M, group_order, s_gfdm);
rx_sym_gfdm_zf = mf_gfdm(gfdmPara.id_data);
% scatterplot(rx_sym_gfdm_zf); title('ZF GFDM');
[~, rx_bit_gfdm_zf] = bbSymDemap(modType, rx_sym_gfdm_zf);
[~, ber_gfdm_zf] = biterr(tx_bit_gfdm, rx_bit_gfdm_zf);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% generate baseband symbols for GFDM/OQAM
% -------------------------------------------------------------------------
% compute AD coefficient matrix and the positions of AD and data.
[oqamPara.ADcoeMtx, oqamPara.id_AD, oqamPara.id_data] = ...
    gen_ADcoeMtx_OQAM(g_oqam, K, M, group_order, oqam_method, idp, posAD_oqam);
Nd_oqam = length(oqamPara.id_data)/2;
eta_oqam = Nd_oqam/(N+Ng); % spectral efficiency
% -------------------------------------------------------------------------
tx_bit_oqam = gen_binbits(modType, Nd_oqam, 1); % generate 0-1 bits
% -------------------------------------------------------------------------
[tx_sym_oqam, ~] = bbSymMap(tx_bit_oqam, modType); % constellation mapping
% -------------------------------------------------------------------------
%% GFDM/OQAM: Tx
% -------------------------------------------------------------------------
xbb_oqam = make_GFDM_OQAMsym_AD(oqamPara, tx_sym_oqam, pilots); % make GFDM symbol
% -------------------------------------------------------------------------
xbb_oqam_g = add_GuardInterval(xbb_oqam, Ng, typeGI); % add CP
% -------------------------------------------------------------------------
xpb_oqam = bb2pb(xbb_oqam_g, bw, fc, fs); % convert to passband
% figure; plot(xpb_oqam); xlim([-100 length(xpb_oqam)+100]); title('GFDM signal');
% -------------------------------------------------------------------------
ypb_oqam = add_channel(xpb_oqam, fs, ch); % pass through channel
ypb_oqam = awgn(ypb_oqam, snr, 'measured');
%% GFDM/OQAM: Rx
% -------------------------------------------------------------------------
zbb_oqam = pb2bb(ypb_oqam, bw, fc, fs); % convert to baseband
% -------------------------------------------------------------------------
% transform into frequeNy domain
z_1blk_oqam = zbb_oqam(1 : N+Ng); % disgard samples outside guard interval
Zf_oqam = ser2para_rmvGI_fft(z_1blk_oqam, N, 1, typeGI, 'shiftOp', false);
%% GFDM/OQAM: channel estimation
% -------------------------------------------------------------------------
% LS channel estimation at pilots
h_LS_oqam = chEst_ofdm_LS(Zf_oqam, pilots, idp);
% -------------------------------------------------------------------------
% Transform domain channel estimation
[cfr_oqam, ~] = chEst_TransDom_LP_IDFT(h_LS_oqam, idp, N, Ng, []);
% -------------------------------------------------------------------------
%% GFDM/OQAM: equalization
% -------------------------------------------------------------------------
% Zero-forcing (ZF) channel equalization by default in this section.
S_oqam = Zf_oqam./cfr_oqam;
% -------------------------------------------------------------------------
% matched filtering (MF) equalization/demodulation 
[~, mf_oqam_i, mf_oqam_q] = mf_FreqDom_GFDM_OQAM...
    (g_oqam, K, M, group_order, oqam_method, S_oqam);
mf_oqam = real([mf_oqam_i(:); mf_oqam_q(:)]);
mf_oqam_data = mf_oqam(oqamPara.id_data);
mf_oqam_data  =reshape(mf_oqam_data, [], 2);
rx_sym_oqam = complex(mf_oqam_data(:,1), mf_oqam_data(:,2));
% scatterplot(rx_sym_oqam); title('MF OQAM');
[~, rx_bit_oqam] = bbSymDemap(modType, rx_sym_oqam);
[~, ber_oqam] = biterr(tx_bit_oqam, rx_bit_oqam);
%% show results
show_results_CFBMC_AD;