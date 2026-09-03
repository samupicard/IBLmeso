function fig_leg = plot_hsv_fraction_legend(fraclim, vlim, labelstr)

    if nargin < 3
        labelstr = 'frac. sig.';
    end

    n = 250;

    [F, V] = meshgrid(linspace(fraclim(1), fraclim(2), n), ...
                      linspace(vlim(1),   vlim(2),   n));

    Fnorm = min(max((F - fraclim(1)) ./ diff(fraclim), 0), 1);
    Vnorm = min(max((V - vlim(1)) ./ diff(vlim), 0), 1);

    H = 0.70 * (1 - Fnorm);  % blue -> cyan/green/yellow/red
    S = ones(size(H));

    RGB = hsv2rgb(cat(3, H, S, Vnorm));

    fig_leg = figure('Position',[1500,-100,180,100], ...
                     'Name','HSV fraction-density legend');

    image(linspace(fraclim(1), fraclim(2), n), ...
          linspace(vlim(1),   vlim(2),   n), ...
          RGB);

    set(gca,'YDir','normal', 'XColor', 'w', 'YColor', 'w')
    %axis square
    box on 
    set(gcf,'Color','k')

    xlabel(labelstr, 'FontSize', 12,'Color',[1 1 1])
    %ylabel('nr. ROIs', 'FontSize', 12)
    ylabel('nr. sess', 'FontSize', 12,'Color',[1 1 1])
    %title('Hue = fraction, value = neuron count', 'FontSize', 13)
end