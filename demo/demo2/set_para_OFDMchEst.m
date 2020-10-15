%% Set OFDM parameters
%% Basic settings
modType = 'QPSK';
typeGI = 'cp';
%% Set bandwidth, center frequency, and guard time, etc.
ofdm.Nc = 1024; % total number of subcarriers
ofdm.bw = 6e3;
ofdm.fc = 12e3;
ofdm.Ng = round(ofdm.Nc/8);
ofdm.Tg = ofdm.Ng/ofdm.bw;
fs_tx = 96e3; % sampling frequency of the transmitted signal in passband
%% channel parameters
chanPara.Npa = 7; % the number of paths,
chanPara.interarrv_mean = 2e-3; % the mean of inter-arrvals (in seconds)
chanPara.maxAvgPathGainLoss = -20; % average path gain loss (in dB)
chanPara.TimeVar = 1; % channel is time-varying
chanPara.UniDopplerFactor = 0; % path-specific Doppler
chanPara.sigma_v = 0.1; % motion speed (m/s)

snr = 15; 
%% Set positions of data, pilots, and null subcarriers 
% each OFDM block has 12 null subcarriers on both band edges, and 24 null
% carriers randomly distributed within the band.
Nz_edge = 12;
B_coh = 1/ofdm.Tg;
F = ofdm.bw/ofdm.Nc;
P = floor(B_coh/F);
% P = 4;
ofdm.idp = (Nz_edge+1: P : ofdm.Nc-Nz_edge).';
Np = length(ofdm.idp);

id_d_z = setdiff((Nz_edge+1:ofdm.Nc-Nz_edge).', ofdm.idp);
idz_inband = randperm(length(id_d_z), 24);
idz_inband = id_d_z(idz_inband);
ofdm.idz = union([1:Nz_edge ofdm.Nc-Nz_edge+1:ofdm.Nc].', idz_inband);
Nz = length(ofdm.idz);

ofdm.idd = setdiff((1:ofdm.Nc).', union(ofdm.idp, ofdm.idz));
Nd = length(ofdm.idd);
%% Assign transmission window
beta = 1/8;
ofdm.window = rcos_window(ofdm.Nc, beta); % raised cosine window
% ofdm.window = ones(ofdm.Nc,1); % traditional rectangular window
if length(ofdm.window)>ofdm.Nc && ~strcmpi(typeGI,'zp')
    warning('Window applied. Guard type changed to ZP');
    typeGI = 'zp';
end


