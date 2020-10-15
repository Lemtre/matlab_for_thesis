%% Table of performance results
fprintf('\n BER performance:\n\n');
fprintf('Uncoded BER | Conventional interpolation | Transform domain method \n')
fprintf('------------+----------------------------+-------------------------\n')
fprintf('    ZF      |          %5.4f            |     %5.4f    \n', ber1, ber2 );
fprintf('------------+----------------------------+-------------------------\n')
fprintf('   MMSE     |          %5.4f            |     %5.4f    \n\n', ber3, ber4 );

%% scatterplots
figure; set(gcf,'unit','centimeters','position',[27 42 30 30]);
subplot(2,2,1); 
plot(real(rx_sym1), imag(rx_sym1),...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',2);
grid on; axis equal;
title('Interp. + ZF');

subplot(2,2,2); 
plot(real(rx_sym2), imag(rx_sym2),...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',2);
grid on; axis equal;
title('TransDom. + ZF');

subplot(2,2,3); 
plot(real(rx_sym3),imag(rx_sym3),...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',2);
% scatterplot(s_BP(ofdm.idd)); 
grid on; axis equal;
title('Interp. + MMSE');

subplot(2,2,4); 
plot(real(rx_sym4),imag(rx_sym4),...
    'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',2);
% scatterplot(s_BP(ofdm.idd)); 
grid on; axis equal;
title('TransDom. + MMSE');
