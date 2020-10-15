# MATLAB code for doctoral thesis 
Part of MATLAB code for doctoral thesis, "Multicarrier Communications in Doubly Dispersive Underwater Acoustic Channels", by Jing Tian.
## func_JingTian folder under functions
FUNC_JINGTIAN contains some commonly needed funtctions in multicarrier communications:
* OFDM modulation/demodulation modules, including several data-aided channel estimation methods and differential demodulation;
* fast implementation of GFDM and C-FBMC modulation/demodulation, including zero-forcing (ZF) and matched filtering (MF) equalization in time and frequency domain respectively;
* functions that calculate the dictionaries for compressed sensing approach used in OFDM channel estimation;
* functions that add channel effect, applying broadband Doppler distortion;
* resample a signal in frequency domain with arbitrary precision via chirp-Z transform (CZT);
* calculation of cross ambiguity functions between signals synthesized with Hermite functions.
## demo folder
There are four demos as listed below
* demo1: sparse 2-D channel estimation for OFDM using a compressed sensing approach in form of basis-pursuit (BP) algorithm;
* demo2: comparison between basic data-aided OFDM channel estimation algorithms, including conventional frequency domain interpolation and IDFT based transform domain method;
* demo3: DFT-eigenvector-based prototype filter synthesis for circular filterbank multicarriers (C-FBMC/OQAM); 
* demo4: auxiliary date aided interference-free pilot synthesis for generalized frequency division multiplexing (GFDM).
## Usage
* step1: run "Run_me_first.m";
* step2: run files in the DEMO folder named as "demo + number";
## Reference
See "Read_me.pdf".
