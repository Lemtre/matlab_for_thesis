%% Table of performance results
fprintf('\n BER performance:\n\n');
fprintf('            |    known CSI    |      BP      \n')
fprintf('------------+-----------------+--------------\n')
fprintf('Uncoded BER |      %5.4f     |     %5.4f    \n\n', ber_CSI, ber_BP );
%% estimated 2-D sparse channel contour plot
ch_BPest = abs(xi).^2;
ch_BPest = ch_BPest./max(ch_BPest);
ch_2d = reshape(ch_BPest, length(delay_set), []);
ch_2d = 10*log10(ch_2d);
c_vec = linspace(0,-20,41);

figure; set(gcf,'unit','centimeters','position',[27 42 30 30]);
subplot(2,2,1); 
[~,h] = contour(Doppler_set, delay_set, ch_2d, c_vec);
h.LineWidth = 2.7;
xlim([-5e-4 5e-4]);
xlabel('Relative Doppler factor');
ylim([-5e-3 ofdm.Tg+5e-3]);
ylabel('Delay (ms)');
grid on;
title('2-D channel estimation');
c = colorbar;
caxis([-20 0]);
xlabel(c,'dB');
%% scatterplots
subplot(2,2,2); 
plot(real(Zf_0ovrsamp(ofdm.idd)),imag(Zf_0ovrsamp(ofdm.idd)),...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',3);
% scatterplot(Zf_0ovrsamp(ofdm.idd)); 
grid on; axis equal;
title('Unequalized scatterplot');

subplot(2,2,3); 
plot(real(s_CSI(ofdm.idd)),imag(s_CSI(ofdm.idd)),...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',3);
% scatterplot(s_CSI(ofdm.idd)); 
grid on; axis equal;
title('MMSE equalization with known CSI');

subplot(2,2,4); 
plot(real(s_BP(ofdm.idd)),imag(s_BP(ofdm.idd)),...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',3);
% scatterplot(s_BP(ofdm.idd)); 
grid on; axis equal;
title('MMSE equalization with BP estimation.');


