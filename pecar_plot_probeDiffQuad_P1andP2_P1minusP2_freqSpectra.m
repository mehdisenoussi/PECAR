pecar_loc = './';
save_loc = [pecar_loc, 'results/'];

delays = 40:40:520;
n_del = size(delays, 2);

n_obs = 11;

byvalidity = true; bycongru = true; onlycorrect = true;

txtval = ''; txtcongru = ''; txtcorrect = '';
if byvalidity; vals = 1:2; txtval = '_byvalidity'; else vals = [1]; end
if bycongru; congrus = 1:3; txtcongru = '_bycongrus'; else congrus = [1]; end
if onlycorrect; txtcorrect = '_onlycorrect'; end

n_congrus = size(congrus, 2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P1 and P2
ylims = [-.1, .77];

vals = 2:-1:1;

n_obs_data_filename = [save_loc, sprintf('%iobs_P1_P2_Delta%s%s%s',...
    n_obs, txtval, txtcongru, txtcorrect)];
load(n_obs_data_filename)

lineW = 2;
markersize = 8;
perc = .05; ts = tinv(1-perc/2., n_obs - 1);

% plot the "probe on different quadrant" congruency condition
cong = 1;

figure()
plotn = 1;
for val = vals
    if val==1; cols = [230, 60, 23; 230, 138, 23] / 256.;
    else cols = [51, 78, 198; 51, 180, 255] / 256.; end
    
    subplot(3, 2, plotn); grid;
    plotn = plotn + 1; hold on;

    avgstoplot = [nanmean(P1_all(:, val, cong, :), 4),...
        nanmean(P2_all(:, val, cong, :), 4)];
    barDir = avgstoplot(:, 1) > avgstoplot(:, 2);
    
    for Pind = 1:2
        plot(delays, avgstoplot(:, Pind), 'o-', 'Color', cols(Pind, :), ...
        'MarkerEdgeColor', cols(Pind, :), 'MarkerFaceColor', [1, 1, 1],...
        'LineWidth', lineW, 'MarkerSize', markersize);
    end
    legend('p1', 'p2', 'Location', 'NorthWest')
    
    Pind = 1;
    for Pall = {P1_all, P2_all}
        Pall = Pall{1};

        % compute SEM
        SEM = nanstd(Pall(:, val, cong, :), [], 4) ./ sqrt(n_obs);
        % compute 95% confidence interval
        % CI = ts * SEM;
        CI = SEM;

        if Pind == 2; barDir = ~barDir; end
        for delind = 1:length(delays)
            if barDir(delind)
                ys = [avgstoplot(delind, Pind),...
                    avgstoplot(delind, Pind) + CI(delind)];
            else
                ys = [avgstoplot(delind, Pind) - CI(delind),...
                    avgstoplot(delind, Pind)];
            end
            plot([delays(delind), delays(delind)], ys, 'Color',...
                cols(Pind, :), 'LineWidth', lineW);
        end
        Pind = Pind + 1;
    end
    
    if cong == 1
        if val == 2; title('Valid');
        elseif val == 1; title('Invalid');
        end
    end
    
    
    set(gca, 'YTick', -.1:.2:.75, 'FontSize', 13, 'LineWidth', 2', 'Fontname', 'Ariel')
    set(gca, 'XTick', 0:100:500, 'FontSize', 13, 'LineWidth', 2', 'Fontname', 'Ariel')

    xlabel('Time from search task onset [ms]','FontSize',12,'Fontname','Ariel')
    
    if val == 2
        ylabel(sprintf('Probe report probabilities'), 'FontSize', 12,...
            'Fontname', 'Ariel')
    end

    ylim(ylims); xlim([0, 540])
    
end


%% P1 minus P2
Pdiff = P1_all - P2_all;

ylims = [-.5, .7];

for val = vals
    subplot(3, 2, plotn); plotn = plotn+1; hold on;        

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
end


%%

load([save_loc sprintf('fft_Pdiff_byObs_p_pad_2s_valid%s_congru1_%isubjs.mat',...
    txtcorrect, n_obs)])
load([save_loc sprintf('fft_Pdiff_byObs_p_pad_2s_invalid%s_congru1_%isubjs.mat',...
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

ymin = .4; ymax = 2.5;
xmin = 0; xmax = 12.5;
for val = vals
    subplot(3, 2, plotn); plotn = plotn + 1; hold on;

    if val==1; col = [1, .4, .4];
    else col = [.4, .4, 1]; end

    allftimefunction = sort(squeeze(fft_p_all(:, :, val, cong)), 1);
    percentile = floor(ci*repeatnumber);
    upperlim = allftimefunction(percentile, :);
    percentile2 = floor(ci2*repeatnumber);
    upperlim2 = allftimefunction(percentile2, :);
    expected = mean(fft_p_all(:, :, val, cong), 1);

    surf(xfreq, 1:100000, allftimefunction);
    hold on; view(0,0); colormap(C); set(gcf, 'Renderer', 'Zbuffer');
    shading interp;

    plot3(xfreq, zeros(size(xfreq, 2)), expected, 'k-', 'LineWidth', 1.5);
    plot3(xfreq, zeros(size(xfreq, 2)), upperlim, 'k--', 'LineWidth', 1.5);
    plot3(xfreq, zeros(size(xfreq, 2)), upperlim2, '--', 'Color', [.4, .4, .4], 'LineWidth', 1.5);
    avgtoplot = squeeze(nanmean(a_fft_pdiff_pad(val, cong, :, :), 4));
    plot3(xfreq, zeros(size(xfreq, 2)), avgtoplot,...
        'ko-', 'LineWidth', 2, 'MarkerFaceColor', [1, 1, 1],...
        'MarkerSize', markersize, 'Color', [.2, .2, .2])
    zlim([ymin, ymax]);
    xlim([xmin, xmax])

    grid

    set(gca, 'ZTick', .0:.5:2.5, 'FontSize', 13, 'LineWidth', 2', 'Fontname', 'Ariel')
    set(gca, 'XTick', xfreq(2:2:end), 'FontSize', 13,'LineWidth', 2', 'Fontname', 'Ariel')
    zlabel('Amplitude (a.u)');
    xlabel('Frequencies (Hz)');
end

if onlycorrect; trialmask = 'only correct trials';
else trialmask = 'all trials'; end

% this is commented out because it isn't always available depending on the
% matlab version and toolboxes you have.
% suptitle(sprintf(['Probability estimates results for probes on different quadrants - %i subj\n'...
%     'p<0.05 (bonferroni corrected) - mask: %s'], n_obs, trialmask));
