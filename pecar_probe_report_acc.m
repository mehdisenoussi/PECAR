%% This script computes accuracy for probe reports depending on probe positions
pecar_loc = '/Users/mehdisenoussi/work/lpp/pecar/';
data_loc = [pecar_loc, 'pecar_data/'];
save_loc = '~/Dropbox/postphd/lpp/pecar/code/last_version/clean_forgit/';
% load all variable containing the "raw" data
% (i.e. correct/incorrect, etc. for each trial not P1 and P2 estimates)
load([save_loc, '11obs_probe_grat_delay_respCue_validity_congruency_respGratCorr_infos.mat'])
n_obs = size(probe_info_all, 2);

% only trials for which the observer was correct on the 2AFC task?
only_correct = 1;
if only_correct
    txtcorrect = 'correctTrials';
else
    txtcorrect = 'allTrials';
end
% delas in milliseconds
delays = 40:40:520;

counter_ntrials = [];
probe_results_all = {};

% figure();

%filters = {'diff', 'diffNoTargetPosNoDistrPos'};
probes_side = 'diff';

%probes_side = probes_side{1};

% dimensions: observers, attended/unattended, valid/invalid, delay
probe_results = zeros(n_obs, 2, 2, 13);

for obs_i = 1:n_obs
    n_trials = size(respCue_all{obs_i}, 2);
    % Get the probe display side by checking whether its X coordinate
    % was above or below 0 (which represents the center of the screen)
    % Yields a n_trial-by-2 matrix where 1 is left side and 2 is right side
    probes_side_all = squeeze((probe_info_all{obs_i}(:, 3, :) > 0) + 1);

    % Whether each probe was on the attended side (according to the
    % response cue)
    probe_attended = repmat(respCue_all{obs_i}, 2, 1)' == probes_side_all;

    % Whether the probe response was correct
    probe_resp = squeeze(probe_info_all{obs_i}(:, 5, :));
    
    % only keep trials in which probes were not on the same side
    mask_probe_side = probes_side_all(:, 1) ~= probes_side_all(:, 2);

    % If mask on probe-grat location relation then all trials are taken
    % and the only necessary mask is the probe-grat location relation.
    % If it's probe-probe location relation then use the
    % mask_probe_side mask to select only the trials of interest.
    % Also create the (grating) correctness mask depending on
    % only_correct which represent whether to only process trials in
    % which observers correctly responded to the primary task
    if only_correct
        correctness = resp_corr_grat_all{obs_i}';
    else
        correctness = true(n_trials, 1);
    end

    trialmask = mask_probe_side & correctness;
    probe_attended = probe_attended(trialmask, :);
    
    probes_side_all = probes_side_all(trialmask, :);
    probe_resp = probe_resp(trialmask, :);
    delays_per_trial = delays_all{obs_i}(trialmask);
    validity_per_trial = validity_all{obs_i}(trialmask);

    for delay = 1:13
        % mask for this delay
        delaymask = delays_per_trial == delay;
        % mask for valid trials
        valmask = validity_per_trial == 2;

        % create valid & delay mask together and repeat matrix for the two probes
        allmask = delaymask & valmask;
        allmask = repmat(allmask, 2, 1)';
        
        % compute average accuracy for probes displayed at attended side
        if size(probe_resp(logical(probe_attended & allmask)), 1) ~= 0
            probe_results(obs_i, 1, 1, delay) =...
                mean(probe_resp(logical(probe_attended & allmask)));
        else
            probe_results(obs_i, 1, 1, delay) = 0;
        end
        
        % at unattended side
        if size(probe_resp(logical(probe_attended & allmask)), 1) ~= 0
            probe_results(obs_i, 2, 1, delay) =...
                mean(probe_resp(logical(~probe_attended & allmask)));
        else
            probe_results(obs_i, 2, 1, delay) = 0;
        end


        % create invalid & delay mask together and repeat matrix for the two probes
        allmask = delaymask & ~valmask;
        allmask = repmat(allmask, 2, 1)';
        % at attended side (because it's invalid trials the attended
        % probe was on the side not cued by the response cue, therefore
        % attended probes are ~probe_attended (i.e. probes cued by the
        % pre-cue)
        if size(probe_resp(logical(~probe_attended & allmask)), 1) ~= 0
            probe_results(obs_i, 1, 2, delay) =...
                mean(probe_resp(logical(~probe_attended & allmask)));
        else
            probe_results(obs_i, 1, 2, delay) = 0;
        end
        % at unattended side
        if size(probe_resp(logical(probe_attended & allmask)), 1) ~= 0
            probe_results(obs_i, 2, 2, delay) =...
                mean(probe_resp(logical(probe_attended & allmask)));
        else
            probe_results(obs_i, 2, 2, delay) = 0;
        end
    end
end

%% Compute attended-unattended surrogates

xfreq= .5:.5:12; n_freqs = length(xfreq);
% compute permutations
compute_diff_perm = true;
repeatnumber = 100000; ndel = 13;
if compute_diff_perm
    % here 1 is valid and 2 is invalid
    for valind = 1:2
        perfdiff_pad = [];
        fft_accDiff_byObs_eachObs_p_pad = zeros(repeatnumber, n_obs, n_freqs);
        for repeat = 1:repeatnumber
            if ~mod(repeat, 500); disp(['Repeat number: ', num2str(repeat)]); end
            accdiff_pad = [];
            
            for obs_i = 1:n_obs
                rand_del_sample = randsample(1:ndel, ndel);
                
                accdiff = squeeze(probe_results(obs_i, 1, valind, :)) -...
                            squeeze(probe_results(obs_i, 2, valind, :));
                accdiff = accdiff(rand_del_sample);

                meandiff = mean(accdiff);
                accdiff_pad(obs_i, :) =...
                     [repmat(meandiff, 1, 18), accdiff', repmat(meandiff, 1, 19)];
            end
            absAccDiffAvgPad = abs( fft(accdiff_pad, [], 2) );
            fft_accDiff_byObs_eachObs_p_pad(repeat, :, :) = absAccDiffAvgPad(:, 2:size(xfreq, 2) + 1);
        end

        % Compute surrogates for accuracyDiff by shuffling the accuracy by
        % delay independently for each "attention condition" (attended side or 
        % unattended side) and then re-computing a new accuracy difference
        if valind == 1
            txtval = ''; fft_accDiff_byObs_eachObs_p_pad_valid = fft_accDiff_byObs_eachObs_p_pad;
        else txtval = 'in'; fft_accDiff_byObs_eachObs_p_pad_invalid = fft_accDiff_byObs_eachObs_p_pad;
        end

        save([save_loc, sprintf('fft_accDiff_byObs_eachObs_p_pad_2s_%sfilter_%svalid_%isubjs_%s.mat',...
            probes_side, txtval, n_obs, txtcorrect)],...
            sprintf('fft_accDiff_byObs_eachObs_p_pad_%svalid', txtval), 'xfreq')
    end
end


%% Plot diff accuracy and spectra

if only_correct
    txtcorrect = 'correctTrials';
else
    txtcorrect = '';
end

load([save_loc, sprintf('fft_accDiff_byObs_eachObs_p_pad_2s_%sfilter_valid_%isubjs_%s.mat',...
            probes_side, n_obs, txtcorrect)])
load([save_loc, sprintf('fft_accDiff_byObs_eachObs_p_pad_2s_%sfilter_invalid_%isubjs_%s.mat',...
            probes_side, n_obs, txtcorrect)])
% here we put valid before invalid (whereas it's the inverse in P1-P2
% scripts, because in this script index 1 is valid and index 2 is invalid
fft_p_all = cat(3, squeeze(mean(fft_accDiff_byObs_eachObs_p_pad_valid, 2)),...
    squeeze(mean(fft_accDiff_byObs_eachObs_p_pad_invalid, 2)));
ci = 1 - (.05/6);
ci2 = 1 - .05;

% what percentile of confidence interval do you want to plot
perc = .05;

n_del = size(delays, 2);

% load custom black to green colormap
load('custom_black2green_cmap.mat')


yminspectra = [.12, .3];
ymaxspectra = [.78, 1.3];
ylims_diff = [-.35 .35];

plot_gradient = true;
plot_only_surrs_gradient = false;
plot_surrs_lvl = true;
figure('Position', [118 65 715 740]);
for valind = 1:2
    if valind == 1
        valtxt = '';
    elseif valind == 2
        valtxt = 'in';
    end
    
    % plot the difference
    subplot(2, 2, valind); hold on;
    temp = squeeze(probe_results(:, 1, valind, :)) - squeeze(probe_results(:, 2, valind, :));
    avgtoplot = mean(temp, 1);
    % compute SEM
    SEM = nanstd(temp, [], 1) / sqrt(n_obs);
    % compute confidence interval
    ts = tinv([perc/2., 1 - perc/2.], n_obs-1);
    CI = repmat(avgtoplot, 2, 1) + ts' * SEM;
    plot(delays, avgtoplot, 'ko-', 'MarkerEdgeColor', [0, 0, 0],...
        'MarkerFaceColor', [1, 1, 1], 'LineWidth', 2, 'MarkerSize', 10);
    for delind = 1:n_del
        plot([delays(delind), delays(delind)], [CI(1, delind); CI(2, delind)], 'k', 'LineWidth', 2)
    end
    if valind == 1
        title('Valid');
    else
        title('Invalid');
    end
    
    xlabel('Time from search task onset [ms]','FontSize',12,'Fontname','Ariel')
    ylabel(sprintf('Probe report accuracy');
    
    ylim(ylims_diff); xlim([0, 540])
    grid()


    % do FFT on the difference
    for obs_i = 1:n_obs
        %temp = squeeze(probe_results_all{target_distr_ind}(obs_i, attind, valind, :));
        temp = squeeze(probe_results(obs_i, 1, valind, :))...
            - squeeze(probe_results(obs_i, 2, valind, :));
        avg_perfdiff_pecar_obs_i = mean(temp);
        perfdiff_pad_pecar(obs_i, valind, :) =...
         [repmat(avg_perfdiff_pecar_obs_i, 1, 18), temp',...
         repmat(avg_perfdiff_pecar_obs_i, 1, 19)];
    end
    fft_perfdiff_pad_pecar = fft(perfdiff_pad_pecar, [], 3);
    a_fft_perfdiff_pad = abs(fft_perfdiff_pad_pecar(:, :, 2:size(xfreq, 2) + 1));
    
    % plot spectra and surrogates
    subplot(2, 2, valind + 2); hold on;
    
    % compute probability levels from surrogates
    allftimefunction = sort(squeeze( fft_p_all(:, :, valind)), 1);
    percentile = floor(ci*repeatnumber);
    upperlim = allftimefunction(percentile, :);
    percentile2 = floor(ci2*repeatnumber);
    upperlim2 = allftimefunction(percentile2, :);
    expected = nanmean(fft_p_all(:, :, valind), 1);
    
    surf(xfreq, 1:100000, allftimefunction);
    hold on; view(0,0); colormap(C); set(gcf, 'Renderer', 'Zbuffer');
    shading interp;
    if plot_surrs_lvl
        plot3(xfreq, zeros(size(xfreq, 2)), upperlim, 'k--', 'LineWidth', 1.5);
        plot3(xfreq, zeros(size(xfreq, 2)), upperlim2, '--', 'Color', [.4, .4, .4], 'LineWidth', 1.5);
        avgtoplot = squeeze(nanmean(a_fft_perfdiff_pad(:, valind, :), 1));

        plot3(xfreq, zeros(size(xfreq, 2)), avgtoplot,...
            'ko-', 'LineWidth', 2, 'MarkerFaceColor', [1, 1, 1],...
            'MarkerSize', 8, 'Color', [.2, .2, .2])            
    end
    zlim([yminspectra(valind), ymaxspectra(valind)]);
        
        
end

suptitle(sprintf('Probe report accuracy difference analysis for probe on %s side', probes_side));




