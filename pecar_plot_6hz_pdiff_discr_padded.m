% the PECAR data are sampled every 40ms for 480ms. So we add 18 time points
% (which reprensents 720ms) before and 19 time points after (760ms),
% yielding a 2030ms segment
n_obs = 11;
load([save_loc, sprintf('%iobs_P1_P2_Delta', n_obs)])

pdiff_pecar = P1_all - P2_all;
discr_pecar = Delta_all;
n_obs = size(Pdiff, 3);

discr_pad_pecar = []; pdiff_pad_pecar = [];
for val = 1:2
    for obs_i = 1:n_obs
        avg_pdiff_pecar_obs_i = mean(pdiff_pecar(:, val, obs_i), 1);
        pdiff_pad_pecar(val, :, obs_i) = [repmat(avg_pdiff_pecar_obs_i, 1, 18) ...
                        pdiff_pecar(:, val, obs_i)' ...
                        repmat(avg_pdiff_pecar_obs_i, 1, 19)];

        avg_d_pecar_obs_i = mean(discr_pecar(:, val, obs_i), 1);
        discr_pad_pecar(val, :, obs_i) = [repmat(avg_d_pecar_obs_i, 1, 18) ...
                        discr_pecar(:, val, obs_i)' ...
                        repmat(avg_d_pecar_obs_i, 1, 19)];
    end
end

sampled_padfreqs_pecar = .5:.5:12;

fft_pdiff_pad_pecar = fft(pdiff_pad_pecar, [], 2);
fft_discr_pad_pecar = fft(discr_pad_pecar, [], 2);

a_fft_pdiff_pad_pecar = abs(fft_pdiff_pad_pecar(:,...
    2:size(sampled_padfreqs_pecar, 2) + 1, :));
a_fft_discr_pad_pecar = abs(fft_discr_pad_pecar(:,...
    2:size(sampled_padfreqs_pecar, 2) + 1, :));


%% script to plot 6hz amplitude in PDiff and discriminant for valid/invalid

figure; hold on
freqind = 12;

% P1-P2 difference
subplot(1, 2, 1); hold on;
inds = [2, 1];
for obs_i = 1:n_obs
    plot([1, 2], [a_fft_pdiff_pad_pecar(2, freqind, obs_i),...
        a_fft_pdiff_pad_pecar(1, freqind, obs_i)], '-', 'Color', [.7, .7, .7]);
end

for val = 2:-1:1
    for obs_i = 1:n_obs
        plot(inds(val), a_fft_pdiff_pad_pecar(val, freqind, obs_i), 'ko',...
            'MarkerFaceColor', [1, 1, 1], 'MarkerSize', 6, 'Color', [.4, .4, .4])
    end
    
    meantoplot = squeeze(mean(a_fft_pdiff_pad_pecar(val, freqind, :), 3));
    % compute SEM
    SEM = std(a_fft_pdiff_pad_pecar(val, freqind, :), [], 3) / sqrt(n_obs);
    % compute 95% confidence interval
    perc = .05;
    ts = tinv([perc/2., 1 - (perc/2.)], 11-1);
    CI = meantoplot + ts * SEM;
    
    plot(inds(val), meantoplot, 'ko', 'MarkerEdgeColor', [0 0 0],...
            'MarkerFaceColor', [0 0 0], 'LineWidth', 1.5, 'MarkerSize', 10);
    plot([inds(val), inds(val)], CI, 'k', 'LineWidth', 2)
    
end
xlim([.75 2.25]); ylim([0 2.2])
[h_pdiff,p_pdiff,cip_diff,stat_pdiff] = ttest(squeeze(a_fft_pdiff_pad_pecar(1, freqind, :)),...
    squeeze(a_fft_pdiff_pad_pecar(2, freqind, :)));
text(.8, 2.15, 'Pdiff diff between valid and invalid:', 'fontsize', 14)
text(.8, 2, sprintf('t=%.2f, p=%.4f\n',stat_pdiff.tstat, p_pdiff/2.), 'fontsize', 14)
set(gca,'xtick',1:2);
set(gca,'xticklabel',{'Valid','Invalid'}, 'fontsize', 14);




% Discriminant
subplot(1, 2, 2); hold on;
inds = [2, 1];
for obs_i=1:n_obs
    plot([1,2],[a_fft_discr_pad_pecar(2, freqind, obs_i)...
         a_fft_discr_pad_pecar(1, freqind, obs_i)],'-', 'Color', [.7 .7 .7]);
end

for val = 2:-1:1
    for obs_i = 1:n_obs
        plot(inds(val), a_fft_discr_pad_pecar(val, freqind, obs_i), 'ko',...
            'MarkerFaceColor', [1 1 1], 'MarkerSize', 6, 'Color', [.4 .4 .4])
    end
    
    meantoplot = squeeze(mean(a_fft_discr_pad_pecar(val, freqind, :), 3));
    % compute SEM
    SEM = std(a_fft_discr_pad_pecar(val, freqind, :), [], 3) / sqrt(n_obs);
    % compute 95% confidence interval
    perc = .05;
    ts = tinv([perc/2., 1 - (perc/2.)], 11-1);
    CI = meantoplot + ts * SEM;
    
    plot(inds(val), meantoplot, 'ko', 'MarkerEdgeColor', [0 0 0],...
            'MarkerFaceColor', [0 0 0], 'LineWidth', 1.5, 'MarkerSize', 10);
    plot([inds(val), inds(val)], CI, 'k', 'LineWidth', 2)
end

xlim([.75 2.25]); ylim([0 .9])
[h_discr, p_discr, ci_discr, stat_discr] = ttest(squeeze(...
    a_fft_discr_pad_pecar(1, freqind, :)),...
    squeeze(a_fft_discr_pad_pecar(2, freqind, :)));
text(.8, .88, 'Discr diff between valid and invalid:', 'fontsize', 14)
text(.8, .82, sprintf('t=%.2f, p=%.4f\n',stat_discr.tstat, p_discr/2.), 'fontsize', 14)
set(gca,'xtick',1:2);
set(gca,'xticklabel',{'Valid','Invalid'}, 'fontsize', 14);

if onlycorrect; txtcorrect_title = 'Only correct trials';
else txtcorrect_title = 'All trials'; end

suptitle(sprintf(['Amplitude of 6hz component for P1-P2 difference ' ...
    'and Discriminant\nin valid versus invalid conditions\n%s'], txtcorrect_title));


% print it in the command line also
fprintf(sprintf('Pdiff diff between valid and invalid: t=%.2f, p=%.4f\n',...
    stat_pdiff.tstat, p_pdiff))
fprintf(sprintf('Discr diff between valid and invalid: t=%.2f, p=%.4f\n',...
    stat_discr.tstat, p_discr))



