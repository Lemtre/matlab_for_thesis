%% Set OFDM parameters for BP method
%% Basic settings
modType = 'QPSK';
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
chanPara.sigma_v = 0.2; % motion speed (m/s)
maxDopplerFactor = sqrt(3)*chanPara.sigma_v/1500; % max relative Doppler Factor 

snr = 20; 
%% Set positions of data, pilots, and null carriers 
ofdm.ovrsamp = 2; % upsampling factor

Nez = 24; % number of null carriers at the band edges
N_inband = ofdm.Nc-Nez;
id_seg = 0:20:N_inband-1;
idz = [1;7];
ofdm.idz = bsxfun(@plus, idz, id_seg);
ofdm.idz = union(ofdm.idz(:), ofdm.Nc-Nez+1 : ofdm.Nc);

idz_ovrsamp = bsxfun(@plus, ofdm.ovrsamp*(ofdm.idz-1), 1:ofdm.ovrsamp);
idz_ovrsamp = sort(idz_ovrsamp(:)); % positions of null carriers if upsampled

idp = 2:6;
ofdm.idp = bsxfun(@plus, idp(:), id_seg);
ofdm.idp = sort(ofdm.idp(:));
Np = numel(ofdm.idp);

idd = 8:20;
ofdm.idd = bsxfun(@plus, idd(:), id_seg);
ofdm.idd = sort(ofdm.idd(:));
Nd = numel(ofdm.idd);


ofdm.idobz = bsxfun(@plus, ofdm.ovrsamp*(ofdm.idp-1), 1:ofdm.ovrsamp);
ofdm.idobz = sort(ofdm.idobz(:)); % positions of observations
%% Assign transmission window
beta = 1/8;
ofdm.window = rcos_window(ofdm.Nc, beta); % raised cosine window
% ofdm.window = ones(ofdm.Nc,1); % traditional rectangular window
%% ICI depth
D = 5;


