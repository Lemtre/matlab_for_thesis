%% Table of performance results
fprintf('\n BER performance:\n\n');
fprintf('            |    GFDM, MF    |    GFDM, ZF    |    OQAM, MF    \n')
fprintf('------------+----------------+----------------+----------------\n')
fprintf('Uncoded BER |     %5.4f     |     %5.4f     |     %5.4f    \n\n',...
    ber_gfdm_mf, ber_gfdm_zf, ber_oqam);

%% scatterplots
figure; set(gcf,'unit','centimeters','position',[27 42 40 10]);
subplot(1,3,1); 
plot(real(rx_sym_gfdm_mf), imag(rx_sym_gfdm_mf),...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',2);
grid on; axis equal;
title('GFDM + MF');

subplot(1,3,2); 
plot(real(rx_sym_gfdm_zf), imag(rx_sym_gfdm_zf),...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',2);
grid on; axis equal;
title('GFDM + ZF');

subplot(1,3,3); 
plot(real(rx_sym_oqam),imag(rx_sym_oqam),...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',2);
grid on; axis equal;
title('OQAM + MF');

