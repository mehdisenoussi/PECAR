% the PECAR data are sampled every 40ms for 480ms. So we add 18 time points
% (which reprensents 720ms) before and 19 time points after (760ms),
% yielding a 2030ms segment
pecar_loc='/Users/mehdisenoussi/work/lpp/pecar/';
data_loc = [pecar_loc, 'pecar_data/'];
save_loc = '~/Dropbox/postphd/lpp/pecar/code/last_version/clean_forgit/';
byvalidity = true; bycongru = true; onlycorrect = true;

txtval=''; txtcongru=''; txtcorrect='';
if byvalidity; vals = 1:2; txtval = '_byvalidity'; else vals = [1]; end
if bycongru; congrus = 1:3; txtcongru = '_bycongrus'; else congrus = [1]; end
if onlycorrect; txtcorrect = '_onlycorrect'; end

n_obs = 11;

n_obs_data_filename = ...
    [save_loc sprintf('%iobs_P1_P2_Delta%s%s%s',...
        n_obs, txtval, txtcongru, txtcorrect)];
load(n_obs_data_filename)


Pdiff = P1_all - P2_all;
pdiff_pecar = squeeze(Pdiff(:, :, 1, :));
discr_pecar = squeeze(Delta_all(:, :, 1, :));
n_obs = size(Pdiff, 4);

discr_pad_pecar = []; pdiff_pad_pecar = [];
for val = 1:2
    for obs_i = 1:11
        avg_pdiff_pecar_obs_i = nanmean(pdiff_pecar(:, val, obs_i), 1);
        pdiff_pad_pecar(val, :, obs_i) = [repmat(avg_pdiff_pecar_obs_i, 1, 18) ...
                        pdiff_pecar(:, val, obs_i)' ...
                        repmat(avg_pdiff_pecar_obs_i, 1, 19)];
    end
    
    for obs_i = 1:11
        avg_d_pecar_obs_i = nanmean(discr_pecar(:, val, obs_i), 1);
        discr_pad_pecar(val, :, obs_i) = [repmat(avg_d_pecar_obs_i, 1, 18) ...
                        discr_pecar(:, val, obs_i)' ...
                        repmat(avg_d_pecar_obs_i, 1, 19)];
    end
end

sampled_padfreqs_pecar= .5:.5:12;

fft_pdiff_pad_pecar = fft(pdiff_pad_pecar, [], 2);
fft_discr_pad_pecar = fft(discr_pad_pecar, [], 2);

a_fft_pdiff_pad_pecar = abs(fft_pdiff_pad_pecar(:, 2:size(sampled_padfreqs_pecar, 2) + 1, :));
a_fft_discr_pad_pecar = abs(fft_discr_pad_pecar(:, 2:size(sampled_padfreqs_pecar, 2) + 1, :));


%% script to plot 4 Hz amplitude in PDiff and discriminant for valid/invalid

figure; hold on
% index 8 is 4 Hz
freqind = 8;

% P1-P2 difference
subplot(1, 2, 1); hold on;
inds = [2,  1];
for obs_i=1:n_obs; plot([1,2],[a_fft_pdiff_pad_pecar(2, freqind, obs_i)...
         a_fft_pdiff_pad_pecar(1, freqind, obs_i)],'-', 'Color', [.7, .7, .7]); end
for val = 2:-1:1
    for obs_i=1:n_obs
        plot(inds(val), a_fft_pdiff_pad_pecar(val, freqind, obs_i), 'ko',...
            'MarkerFaceColor', [1, 1, 1], 'MarkerSize', 6, 'Color', [.4, .4, .4])
    end
    
    meantoplot = squeeze(nanmean(a_fft_pdiff_pad_pecar(val, freqind, :), 3));
    % compute SEM
    SEM = nanstd(a_fft_pdiff_pad_pecar(val, freqind, :), [], 3) / sqrt(n_obs);
    % compute 95% confidence interval
    perc = .05;
    ts = tinv([perc/2., 1 - (perc/2.)], 11-1);
    CI = meantoplot + ts * SEM;
    
    plot(inds(val), meantoplot, 'ko', 'MarkerEdgeColor', [0, 0, 0],...
            'MarkerFaceColor', [0, 0, 0], 'LineWidth', 1.5, 'MarkerSize', 10);
    plot([inds(val), inds(val)], CI, 'k', 'LineWidth', 2)
    
end
xlim([.75, 2.25]); ylim([0, 3.7])
[h_pdiff,p_pdiff,cip_diff,stat_pdiff] = ttest(squeeze(a_fft_pdiff_pad_pecar(1, freqind, :)),...
    squeeze(a_fft_pdiff_pad_pecar(2, freqind, :)));
text(.8, 3.6, 'Pdiff diff between valid and invalid:', 'fontsize', 14)
text(.8, 3.3, sprintf('t=%.2f, p=%.4f\n', stat_pdiff.tstat, p_pdiff/2.), 'fontsize', 14)
set(gca, 'xtick',1:2);
set(gca, 'xticklabel', {'Valid','Invalid'}, 'fontsize', 14);




% Discriminant
subplot(1, 2, 2); hold on;
inds = [2, 1];
for obs_i = 1:n_obs; plot([1, 2],[a_fft_discr_pad_pecar(2, freqind, obs_i)...
         a_fft_discr_pad_pecar(1, freqind, obs_i)],'-', 'Color', [.7, .7, .7]); end
for val = 2:-1:1
    for obs_i = 1:n_obs
        plot(inds(val), a_fft_discr_pad_pecar(val, freqind, obs_i), 'ko',...
            'MarkerFaceColor', [1, 1, 1], 'MarkerSize', 6, 'Color', [.4, .4, .4])
    end
    
    meantoplot = squeeze(nanmean(a_fft_discr_pad_pecar(val, freqind, :), 3));
    % compute SEM
    SEM = nanstd(a_fft_discr_pad_pecar(val, freqind, :), [], 3) / sqrt(n_obs);
    % compute 95% confidence interval
    perc = .05;
    ts = tinv([perc/2., 1 - (perc/2.)], 11-1);
    CI = meantoplot + ts * SEM;
    
    plot(inds(val), meantoplot, 'ko', 'MarkerEdgeColor', [0, 0, 0],...
            'MarkerFaceColor', [0, 0, 0], 'LineWidth', 1.5, 'MarkerSize', 10);
    plot([inds(val), inds(val)], CI, 'k', 'LineWidth', 2)
    

end
xlim([.75, 2.25]); ylim([0, 3])
[h_discr, p_discr, ci_discr, stat_discr] = ttest(squeeze(...
    a_fft_discr_pad_pecar(1, freqind, :)),...
    squeeze(a_fft_discr_pad_pecar(2, freqind, :)));
text(.8, 3, 'Discr diff between valid and invalid:', 'fontsize', 14)
text(.8, 2.8, sprintf('t=%.2f, p=%.4f\n',stat_discr.tstat, p_discr/2.), 'fontsize', 14)
set(gca,'xtick',1:2);
set(gca,'xticklabel',{'Valid','Invalid'}, 'fontsize', 14);

if onlycorrect; txtcorrect_title = 'Only correct trials';
else txtcorrect_title = 'All trials'; end

suptitle(sprintf(['Amplitude of %ihz component for P1-P2 difference ' ...
    'and Discriminant\nin valid versus invalid conditions\n%s'],...
    sampled_padfreqs_pecar(freqind), txtcorrect_title));


% print it in the command line also
fprintf(sprintf('Pdiff diff between valid and invalid: t=%.2f, p=%.4f\n',...
    stat_pdiff.tstat, p_pdiff))
fprintf(sprintf('Discr diff between valid and invalid: t=%.2f, p=%.4f\n',...
    stat_discr.tstat, p_discr))



