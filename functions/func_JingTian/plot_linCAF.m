function [z, H, AXES, COLORMAP] = plot_linCAF(tx, rx, T, rangeT, rangeF, varargin)
% PLOT_LINCAF plot the linear cross-ambiguity function between the synthesis 
% filter and the analysis filter.
% 
%   T: actual time duration;
%   rangeT: interested range of time;
%   rangeF: interested range of frequency;
%   varargin: - upsampF: frequency domain oversampling factor;
%             - plot_opt: plot style options
%   varargout: if output is non-void, pass the complex-valued CAF via "z"
%              and the handle of the plot via "h".
% 
% Author: Jing Tian
% Last modified: 20th April, 2019   
    default.markpeak = 0;
    default.marknull = 0;
    default.symmetry = 0;
    default.zscale = 'log';
    default.contour = 1;
    default.cshading = 0;
    default.showcolorbar = 0;
    switch length(varargin)
        case 0
            upsampF = 4;
            plot_opt = default;
        case 1
            if isstruct(varargin)
                plot_opt = varargin{1};
                upsampF = 4;
            else
                upsampF = varargin{1};
                plot_opt = default;
            end
        case 2
            upsampF = varargin{1};
            plot_opt = varargin{2};
    end
    
    L1 = length(tx);
    L2 = length(rx);
    L = L1 + L2 -1;
    N = round(T*upsampF);
    if N < L
        Nfft = 2.^ceil(log2(L));
        upsampF = L/T;
    else
        Nfft = N;
    end
    M = 2*round(rangeF*upsampF)+1;
    
    z = lin_caf(tx, rx);
    zs = abs(z);
    
    hold on
    if strcmp(plot_opt.zscale, 'log') 
        zs = 20*log10(zs);
        if plot_opt.cshading == 0
            ZS = linspace(-20,0,50);
            ZS1 = linspace(-40,-20,27);
            ZS2 = linspace(-60,-40,42);
            ZS = unique([ZS ZS1 ZS2],'sorted');
%             ZS = linspace(-60,0,100);
        elseif plot_opt.cshading == 1
            ZS = linspace(-60,0,60);
        end
    elseif strcmp(plot_opt.zscale, 'lin')
            plot_opt.cshading = 0;        
            ZS = linspace(1e-3,1,100);
    end
    
    if plot_opt.contour == 1
        if plot_opt.cshading == 1
            contourf(zs,ZS);
        else
            contour(zs,ZS); brighten(-0.15)
        end
        box on; grid on; %set(gca,'xgrid','off');
    else
        surfc(zs); shading interp; alpha(0.5);
        zlim([ZS(1) ZS(end)]);
    end
    caxis([ZS(1) ZS(end)]); 
   
    if plot_opt.markpeak == 1
        hl = line([L1 L2],[ceil(M/2) ceil(M/2)]);
        hl.LineStyle = '-';
        hl.Color = 'r';
        hl.LineWidth = 1.2;
        hs = scatter(round((L1+L2)/2), ceil(M/2), 'pentagram');
        hs.MarkerFaceColor = 'k';
        hs.LineWidth = 0.3;
        hs.SizeData = 63;
    end
    
    if plot_opt.marknull == 1
        xnull = ceil(L/2);
        ynull = 1:upsampF:M;
        ynull(ceil(length(ynull)/2)) = [];
        [xn,yn] = meshgrid(xnull,ynull);
        s = scatter(xn(:),yn(:),'o');
        s.MarkerEdgeColor = 'k';
        s.MarkerFaceColor = 'g';
        s.LineWidth = 1;
    end
    
    if plot_opt.showcolorbar == 1
        colorbar;
    end
    
    if nargout > 4
        error('Too many output arguments!');
    end
    
    %% Axes settings
    t_cent = (L+1)/2;
    tt = [t_cent-rangeT*T  t_cent+rangeT*T];
    xlim(tt);
    if plot_opt.symmetry == 1
        xticklabel = cell(2*floor(rangeT)+1,1);
        for tk = -floor(rangeT) : floor(rangeT);
            xticklabel{tk+floor(rangeT)+1} = num2str(tk);
        end
        xtick = t_cent-floor(rangeT)*T : T: t_cent+floor(rangeT)*T;
    else
        xtick = unique([L1-T L1 L2 L2+T]);
        if L2>L1
            xticklabel = {'$T$','$0$','$T_{\rm g}$','$T+T_{\rm g}$'};
        else
            xticklabel = {'$-T$','$0$','$T$'};
        end
    end
    set(gca,'XTick',xtick);
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'xticklabel',xticklabel);
    
    ylim([1 M]);
    if rangeF > 4
        yedge = min(rangeF*upsampF, floor(M/2));
        ytick = unique([0:2*upsampF:yedge, 0:-2*upsampF:-yedge], 'sorted') + ceil(M/2);
        Num_ytick = length(ytick);
        yticklabel = cell(Num_ytick,1);
        for fk = 1 : Num_ytick
            yticklabel{fk} = 2*fk-Num_ytick-1;
        end
    else
        ytick = 1:upsampF:M;
        yticklabel = cell(2*floor(rangeF)+1,1);
        for fk = -floor(rangeF) : floor(rangeF)
            yticklabel{fk+floor(rangeF)+1} = num2str(fk);
        end
    end
    set(gca,'YTick',ytick);
    set(gca,'yticklabel',yticklabel);
    H = gcf;  
    AXES = gca;
    COLORMAP = colormap;
    %% calculate CAF
    function z = lin_caf(tx, rx)       
        E = sqrt(max(conv(tx, rx)));
        tx = tx ./ E;
        rx = rx ./ E;

        tx0 = [zeros(L-L1,1); tx];
        rx0 = [rx; zeros(L-L2,1)];       
        w = exp(-1j*2*pi/Nfft);
        a = exp(-2j*pi*floor(M/2)/Nfft);
        z = NaN(M,L);
        for ii = 1 : L
            temp = tx0 .* circshift(rx0,ii);
            z(:,ii) = czt(temp, M, w, a);
        end
    end

end