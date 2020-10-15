function [A, At, Af] = plot_result_filtdes(hh, K, M, ratio_null, method)

%% set range and resolution
switch lower(method)
    case 'smt'
        T_max = 2;
        F_max = 4;
        T0 = 1/M;
        F0 = 2/K;
        x_ticklabel = {'-1','0','1'};
        y_ticklabel = {'-2','0','2'};
    case 'cmt'
        T_max = 4;
        F_max = 2;
        T0 = 2/M;
        F0 = 1/K;
        x_ticklabel = {'-2','0','2'};
        y_ticklabel = {'-1','0','1'};
    otherwise
        error('Invalid method.');
end
ovr = 100;
T_max = min(floor(M/2), T_max);
F_max = min(floor(K/2), F_max);

L_T = 2*ovr*T_max+1;
L_F = 2*ovr*F_max+1;
nu = linspace(-F_max/K, F_max/K, L_F);
sigma = linspace(-T_max/M, T_max/M, L_T);

x_tick = [-T0 0 T0];
y_tick = [-F0 0 F0];

%% waveform in time domain and its frequency response
hh = bsxfun(@rdivide, hh, max(abs(hh)));
hf = fftshift(fft(fftshift(hh)));
hf = bsxfun(@rdivide, hf, max(abs(hf)));

figure; set(gcf,'unit','centimeters','position',[27 42 45 30]);
subplot(2,3,1);
stem(real(hh)); 
axis tight; title('Synthesized filter(s)');

subplot(2,3,2);
stem(real(hf)); 
axis tight; title('DFT of the pulse(s)');
%% Calculate paf
if size(hh,2)==1
    [x, y] = deal(hh);
else
    x = hh(:,1);
    y = hh(:,2);
end
A = paf_eqv(x, y, sigma, nu);
A = A./max(abs(A(:))); % normalization
A = 20*log10(abs(A));

At = A(ceil(L_F/2), :);

Af = A(:, ceil(L_T/2));
%% plot paf
subplot(2,3,3);
contourf(sigma, nu, A, 100); 
c = colorbar;
caxis([-60 0]);
xlabel(c, 'dB');
axis tight; axis square; grid on;
set(gca,'TickLabelInterpreter', 'latex');

set(gca,'XTick',x_tick);
set(gca,'xticklabel',x_ticklabel);
set(gca,'YTick',y_tick);
set(gca,'Yticklabel',y_ticklabel);
xlabel('Normalized delay'); ylabel('Normalized frequency');
title('$20\log_{10}(|A(\sigma,\nu)|)$','Interpreter','Latex');


%% plot PAF's projection along the time axis
subplot(2,3,4);
plot(sigma, At); 
axis tight; grid on; hold on;

delta_T = ratio_null/M/2;
tt1 = T0-delta_T;
tt2 = T0+delta_T;
line([tt1 tt1], [min(At) 0], 'LineStyle', '-.');
line([tt2 tt2], [min(At) 0], 'LineStyle', '-.');
line([T0 T0], [min(At) 0], 'LineStyle', ':');
line([2*T0 2*T0], [min(At) 0], 'LineStyle', ':');
line([-tt1 -tt1], [min(At) 0], 'LineStyle', '-.');
line([-tt2 -tt2], [min(At) 0], 'LineStyle', '-.');
line([-T0 -T0], [min(At) 0], 'LineStyle', ':');
line([-2*T0 -2*T0], [min(At) 0], 'LineStyle', ':');
ylim([-80 0]); 
ylabel('$20\log_{10}(|A(\sigma,0)|)$','Interpreter','Latex');

set(gca, 'XTick', x_tick);
set(gca, 'XTickLabel', x_ticklabel);
xlabel('Normalized delay');
%% plot PAF's projection along the frequency axis
subplot(2,3,5);
plot(nu, Af); 
axis tight; grid on; hold on;

delta_F = ratio_null/K/2;
ff1 = F0-delta_F;
ff2 = F0+delta_F;
line([ff1 ff1], [min(Af) 0], 'LineStyle', '-.');
line([ff2 ff2], [min(Af) 0], 'LineStyle', '-.');
line([F0 F0], [min(Af) 0], 'LineStyle', ':');
line([2*F0 2*F0], [min(At) 0], 'LineStyle', ':');
line([-ff1 -ff1], [min(Af) 0], 'LineStyle', '-.');
line([-ff2 -ff2], [min(Af) 0], 'LineStyle', '-.');
line([-F0 -F0], [min(Af) 0], 'LineStyle', ':');
line([-2*F0 -2*F0], [min(At) 0], 'LineStyle', ':');
ylim([-80 0]);
ylabel('$20\log_{10}(|A(0,\nu)|)$','Interpreter','Latex');

set(gca, 'XTick', y_tick);
set(gca, 'XTickLabel', y_ticklabel);
xlabel('Normalized frequency');


end
