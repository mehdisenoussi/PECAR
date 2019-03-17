pecar_loc = './';
save_loc = [pecar_loc, 'results/congruTrick/'];

delays = 40:40:520;
n_del = size(delays, 2);

n_obs = 11;

byvalidity = true; bycongru = true; onlycorrect = true;

vals = 1:2; txtval = '_byvalidity';
congrus = [1]; txtcongru = '_bycongrus';
txtcorrect = '_onlycorrect';

n_congrus = size(congrus, 2);

%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%

vals = 2:-1:1;

n_obs_data_filename = [save_loc, sprintf('%iobs_P1_P2_Delta%s%s%s_congruTrick',...
    n_obs, txtval, txtcongru, txtcorrect)];
load(n_obs_data_filename)

lineW = 2;
markersize = 8;
perc = .05; ts = tinv(1-perc/2., n_obs - 1);

% plot the "probe on same quadrant" congruency condition (merged congruenc
% 2 and 3)
cong = 1;

% P1 minus P2
Pdiff = P1_all - P2_all;

plotn = 1;
ylims = [-.4, .65];
figure;
for val = vals
    subplot(2, 2, plotn); plotn = plotn+1; hold on;        

    avgtoplot = nanmean(Pdiff(:, val, cong, :), 4);
    % compute SEM
    SEM = nanstd(Pdiff(:, val, cong, :), [], 4) ./ sqrt(n_obs);

    % compute 95% confidence interval
    % CI = ts * SEM;
    CI = SEM;

    for delind = 1:length(delays)
        ys = [avgtoplot(delind) - CI(delind),...
                avgtoplot(delind) + CI(delind)];
        plot([delays(delind), delays(delind)], ys, 'Color',...
                'k', 'LineWidth', lineW);
    end
    plot(delays, avgtoplot, 'o-', 'Color', 'k', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1, 1, 1],...
        'LineWidth', lineW, 'MarkerSize', markersize);

    set(gca, 'YTick', -.4:.2:.4, 'FontSize', 13, 'LineWidth', 2', 'Fontname', 'Ariel')
    set(gca, 'XTick', 0:100:500, 'FontSize', 13, 'LineWidth', 2', 'Fontname', 'Ariel')
    
    if val == 2
        ylabel('P1-P2 difference','FontSize',12,'Fontname','Ariel')
    end

    xlabel('Time from search task onset [ms]','FontSize',12,'Fontname','Ariel')
    
    ylim(ylims); xlim([0, 540])
    %axis square
    grid; plot([0 600], [0 0], 'k--')
    
    if cong == 1 && val==2; title('Valid'); elseif cong == 1 && val==1; title('Invalid'); end
end


%%

load([save_loc, sprintf('fft_Pdiff_byObs_p_pad_2s_valid%s_congru1_%isubjs_congruTrick.mat',...
    txtcorrect, n_obs)])
load([save_loc, sprintf('fft_Pdiff_byObs_p_pad_2s_invalid%s_congru1_%isubjs_congruTrick.mat',...
    txtcorrect, n_obs)])
fft_p_all = cat(3, fft_Pdiff_byObs_p_pad_invalid, fft_Pdiff_byObs_p_pad_valid);

repeatnumber = 100000;

ci = 1-(.05/6.);
ci2 = 1-.05;

xfreq = .5:.5:12;
n_freqs = size(xfreq, 2);

% load custom black to green colormap
load('custom_black2green_cmap.mat')

Pdiff = P1_all - P2_all;
pdiff_pad_pecar = [];
for val = 1:2
    for obs_i = 1:11
        avg_pdiff_pecar_obs_i = mean(Pdiff(:, val, cong, obs_i), 1);
        pdiff_pad_pecar(val, cong, :, obs_i) =...
         [repmat(avg_pdiff_pecar_obs_i, 1, 18), Pdiff(:, val, cong, obs_i)',...
                 repmat(avg_pdiff_pecar_obs_i, 1, 19)];
    end
end

fft_pdiff_pad_pecar = fft(pdiff_pad_pecar, [], 3);
a_fft_pdiff_pad = abs(fft_pdiff_pad_pecar(:, :, 2:size(xfreq, 2)+1, :));

ymin = .4; ymax = 2.1;
xmin = 0; xmax = 12.5;
for val = vals
    subplot(2, 2, plotn); plotn = plotn + 1; hold on;

    if val==1; col = [1, .4, .4];
    else col = [.4, .4, 1]; end

    allftimefunction = sort(squeeze(fft_p_all(:, :, val, cong)), 1);
    percentile = floor(ci*repeatnumber);
    upperlim = allftimefunction(percentile, :);
    percentile2 = floor(ci2*repeatnumber);
    upperlim2 = allftimefunction(percentile2, :);
    expected = mean(fft_p_all(:, :, val, cong), 1);

    surf(xfreq, 1:100000, allftimefunction);
    hold on; view(0, 0); colormap(C); set(gcf, 'Renderer', 'Zbuffer');
    shading interp;

    plot3(xfreq, zeros(size(xfreq, 2)), expected, 'k-', 'LineWidth', 1.5);
    plot3(xfreq, zeros(size(xfreq, 2)), upperlim, 'k--', 'LineWidth', 1.5);
    plot3(xfreq, zeros(size(xfreq, 2)), upperlim2, '--', 'Color', [.4, .4, .4], 'LineWidth', 1.5);
    avgtoplot = squeeze(nanmean(a_fft_pdiff_pad(val, cong, :, :), 4));
    plot3(xfreq, zeros(size(xfreq, 2)), avgtoplot,...
        'ko-', 'LineWidth', lineW, 'Color', 'k', 'MarkerFaceColor', [1, 1, 1],...
        'MarkerSize', markersize)
    zlim([ymin, ymax]);
    xlim([xmin, xmax])

    grid

    %axis square
    set(gca, 'ZTick', .0:.5:2.5, 'FontSize', 13, 'LineWidth', 2', 'Fontname', 'Ariel')
    set(gca, 'XTick', xfreq(2:2:end), 'FontSize', 13,'LineWidth', 2', 'Fontname', 'Ariel')
    zlabel('Amplitude (a.u)');
    xlabel('Frequencies (Hz)');
end

if onlycorrect; trialmask = 'only correct trials';
else trialmask = 'all trials'; end
% this is commented out because it isn't always available depending on the
% matlab version and toolboxes you have.
% suptitle(sprintf(['Probability estimates results for probes on same quadrants - %i subj\n',...
%     '"Probes on distractor quadrant" and "Probes on target quadrant combined"\n',...
%     'p<0.05 (bonferroni corrected) - mask: %s'], n_obs, trialmask));



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now plot the two "Probe on same quadrant" congruencies in invalid trials separately
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_loc = [pecar_loc, 'results/'];

n_obs_data_filename = [save_loc, sprintf('%iobs_P1_P2_Delta%s%s%s',...
    n_obs, txtval, txtcongru, txtcorrect)];
load(n_obs_data_filename)

congrus = [2, 3];

n_congrus = size(congrus, 2);

%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%

vals = [1];

n_obs_data_filename = [save_loc, sprintf('%iobs_P1_P2_Delta%s%s%s',...
    n_obs, txtval, txtcongru, txtcorrect)];
load(n_obs_data_filename)

lineW = 2;
markersize = 8;
perc = .05; ts = tinv(1-perc/2., n_obs - 1);

% plot the two "probe on same quadrant" congruencies separately
cong = [2, 3];

% P1 minus P2
Pdiff = P1_all - P2_all;

plotn = 1;
ylims = [-.5, .7];
figure;
for val = vals
    for cong = congrus
        subplot(2, 2, plotn); plotn = plotn+1; hold on;        

        avgtoplot = nanmean(Pdiff(:, val, cong, :), 4);
        % compute SEM
        SEM = nanstd(Pdiff(:, val, cong, :), [], 4) ./ sqrt(n_obs);

        % compute 95% confidence interval
        % CI = ts * SEM;
        CI = SEM;

        for delind = 1:length(delays)
            ys = [avgtoplot(delind) - CI(delind),...
                    avgtoplot(delind) + CI(delind)];
            plot([delays(delind), delays(delind)], ys, 'Color',...
                    'k', 'LineWidth', lineW);
        end
        plot(delays, avgtoplot, 'o-', 'Color', 'k', ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1, 1, 1],...
            'LineWidth', lineW, 'MarkerSize', markersize);

        set(gca, 'YTick', -.4:.2:.4, 'FontSize', 13, 'LineWidth', 2', 'Fontname', 'Ariel')
        set(gca, 'XTick', 0:100:500, 'FontSize', 13, 'LineWidth', 2', 'Fontname', 'Ariel')

        if val == 2
            ylabel('P1-P2 difference','FontSize',12,'Fontname','Ariel')
        end

        xlabel('Time from search task onset [ms]','FontSize',12,'Fontname','Ariel')
        ylabel('P1-P2')
        
        ylim(ylims); xlim([0, 540])

        grid; plot([0 600], [0 0], 'k--')

        if cong == 2
            title('Probe on Distractor Side');
        elseif cong == 3
            title('Probe on Target Side');
        end
    end
end


% Plot spectra

load([save_loc, sprintf('fft_Pdiff_byObs_p_pad_2s_invalid%s_congru2_%isubjs.mat', txtcorrect, n_obs)])
fft_p_all = fft_Pdiff_byObs_p_pad_invalid;
load([save_loc, sprintf('fft_Pdiff_byObs_p_pad_2s_invalid%s_congru3_%isubjs.mat', txtcorrect, n_obs)])
fft_p_all = cat(3, fft_p_all, cat(3, fft_Pdiff_byObs_p_pad_invalid));

repeatnumber = 100000;

ci = 1-(.05/6.);
ci2 = 1-.05;

xfreq = .5:.5:12;
n_freqs = size(xfreq, 2);

% load custom black to green colormap
load('custom_black2green_cmap.mat')

pdiff_pad_pecar = [];
for val = vals
    for cong = congrus
        for obs_i = 1:11
            avg_pdiff_pecar_obs_i = mean(Pdiff(:, val, cong, obs_i), 1);
            pdiff_pad_pecar(val, cong, :, obs_i) =...
             [repmat(avg_pdiff_pecar_obs_i, 1, 18), Pdiff(:, val, cong, obs_i)',...
                     repmat(avg_pdiff_pecar_obs_i, 1, 19)];
        end
    end
end

fft_pdiff_pad_pecar = fft(pdiff_pad_pecar, [], 3);
a_fft_pdiff_pad = abs(fft_pdiff_pad_pecar(:, :, 2:size(xfreq, 2)+1, :));

ymin = .4; ymax = 2.7;
xmin = 0; xmax = 12.5;
for val = vals
    for cong = congrus
        subplot(2, 2, plotn); plotn = plotn + 1; hold on;

        if val==1; col = [1, .4, .4];
        else col = [.4, .4, 1]; end

        allftimefunction = sort(squeeze(fft_p_all(:, :, cong-1)), 1);
        percentile = floor(ci*repeatnumber);
        upperlim = allftimefunction(percentile, :);
        percentile2 = floor(ci2*repeatnumber);
        upperlim2 = allftimefunction(percentile2, :);
        expected = mean(fft_p_all(:, :, cong-1), 1);

        surf(xfreq, 1:100000, allftimefunction);
        hold on; view(0, 0); colormap(C); set(gcf, 'Renderer', 'Zbuffer');
        shading interp;

        plot3(xfreq, zeros(size(xfreq, 2)), expected, 'k-', 'LineWidth', 1.5);
        plot3(xfreq, zeros(size(xfreq, 2)), upperlim, 'k--', 'LineWidth', 1.5);
        plot3(xfreq, zeros(size(xfreq, 2)), upperlim2, '--', 'Color', [.4, .4, .4], 'LineWidth', 1.5);
        avgtoplot = squeeze(nanmean(a_fft_pdiff_pad(val, cong, :, :), 4));
        plot3(xfreq, zeros(size(xfreq, 2)), avgtoplot,...
            'ko-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', [1, 1, 1],...
            'MarkerSize', markersize)
        zlim([ymin, ymax]);
        xlim([xmin, xmax])

        grid

        %axis square
        set(gca, 'ZTick', .0:.5:2.5, 'FontSize', 13, 'LineWidth', 2', 'Fontname', 'Ariel')
        set(gca, 'XTick', xfreq(2:2:end), 'FontSize', 13,'LineWidth', 2', 'Fontname', 'Ariel')
        zlabel('Amplitude (a.u)');
        xlabel('Frequencies (Hz)');
    end
end

if onlycorrect; trialmask = 'only correct trials';
else trialmask = 'all trials'; end

% this is commented out because it isn't always available depending on the
% matlab version and toolboxes you have.
% suptitle(sprintf(['Probability estimates results for probes on same quadrants - %i subj\n'...
%     'p<0.05 (bonferroni corrected) - mask: %s'], n_obs, trialmask));













