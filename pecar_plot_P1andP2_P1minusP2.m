%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vals = 1:2;

% P1 and P2
n_obs_data_filename = [save_loc sprintf('%iobs_P1_P2_Delta', n_obs)];
load(n_obs_data_filename)

plot_order = [2, 1, 4, 3];
figure('Position', get(groot, 'ScreenSize'));
plotn = 1;
for val = vals
    subplot(2, size(vals, 2), plot_order(plotn)); plotn = plotn + 1; hold on;

    Pind = 1;
    for Pall = {P1_all, P2_all}
        if Pind == 1; shift = -5;
        else shift = 5; end
        if val == 1; cols = [230, 60, 23; 230, 138, 23] / 256.;
        else cols = [51, 78, 198; 51, 180, 255] / 256.; end
        Pall = Pall{1};
        meantoplot = mean(Pall(:, val, :), 3);
        % compute SEM
        SEM = std(Pall(:, val, :), [], 3) ./ sqrt(n_obs);
        % compute 95% confidence interval
        perc = .05; ts = tinv([perc/2., 1 - (perc/2.)], n_obs - 1);
        CI = [meantoplot + ts(1) * SEM, meantoplot + ts(2) * SEM];
        plot(delays + shift, meantoplot, 'o-', 'Color', cols(Pind, :), ...
            'MarkerEdgeColor', cols(Pind, :), 'MarkerFaceColor', [1, 1, 1],...
            'LineWidth', 3, 'MarkerSize', 10);
        for delind = 1:length(delays)
            plot([delays(delind) + shift, delays(delind) + shift],...
                [CI(delind, 1), CI(delind, 2)],...
                'Color', cols(Pind, :), 'LineWidth', 3);
        end
        Pind = Pind + 1;
    end

    if val == 2; title('Valid'); elseif val == 1; title('Invalid'); end

    set(gca, 'YTick', -.1:.2:.7, 'FontSize', 13, 'LineWidth', 2',...
        'Fontname', 'Ariel')
    set(gca, 'XTick', 0:100:500, 'FontSize', 13, 'LineWidth', 2',...
        'Fontname', 'Ariel')

    
    ylabel('Probe report probabilities', 'FontSize', 12, 'Fontname', 'Ariel')
    xlabel('Time from search task onset [ms]', 'FontSize', 12, 'Fontname', 'Ariel')

    ylim([0 .7]); xlim([0 540])
end

suptitle(sprintf('P1 & P2 - %i subj', n_obs));


%% P1 minus P2
Pdiff = P1_all - P2_all;

plotindiv = false;
for val = vals
    subplot(2, size(vals, 2), plot_order(plotn)); plotn = plotn + 1; hold on;
    if plotindiv;
        for obs_i = 1:n_obs
            plot(delays, Pdiff(:, val, obs_i), 'ko-', 'LineWidth', 1,...
             'MarkerFaceColor', [1 1 1], 'MarkerSize', 6, 'Color', [0 0 0])
        end
    end

    if val == 1; cols = [230, 60, 23; 230, 138, 23] / 256.;
    else cols = [51, 78, 198; 51, 180, 255] / 256.; end

    meantoplot = mean(Pdiff(:, val, :), 3);
    % compute SEM
    SEM = std(Pdiff(:, val, :), [], 3) ./ sqrt(n_obs);
    % compute 95% confidence interval
    perc = .05; ts = tinv([perc/2., 1 - (perc/2.)], n_obs - 1);
    CI = [meantoplot + ts(1) * SEM, meantoplot + ts(2) * SEM];
    for delind = 1:length(delays)
        plot([delays(delind), delays(delind)],...
            [CI(delind, 1), CI(delind, 2)],...
            'Color', 'k', 'LineWidth', 3);
    end
    plot(delays, meantoplot, 'o-', 'Color', 'k',...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1, 1, 1],...
        'LineWidth', 3, 'MarkerSize', 10);

    set(gca, 'YTick', -.3:.1:.7, 'FontSize', 13, 'LineWidth', 2',...
        'Fontname', 'Ariel')
    set(gca, 'XTick', 0:100:500, 'FontSize', 13, 'LineWidth', 2',...
        'Fontname', 'Ariel')

    if val == 2; title('Valid'); elseif val == 1; title('Invalid'); end
    ylabel('difference P1-P2','FontSize',12, 'Fontname', 'Ariel')
    xlabel('Time from search task onset [ms]', 'FontSize', 12, 'Fontname', 'Ariel')
    
    xlim([0 540]); ylim([-.35 .48]);
end
suptitle(sprintf('P1 minus P2 - %i subj', n_obs));

