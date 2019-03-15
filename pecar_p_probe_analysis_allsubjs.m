addpath(genpath('./mgl-master/'))

pecar_loc = './';
data_loc = [pecar_loc, 'pecar_data/'];
save_loc = [pecar_loc, 'results/'];

delays = 40:40:520;
n_del = size(delays, 2);

observers = ['ym'; 'ac'; 'al'; 'sa'; 'el'; 'gm'; 'hs'; 'hw'; 'js'; 'ma'; 'nv'];
n_obs = size(observers, 1);

byvalidity = true; bycongru = true; onlycorrect = true;

vals = 1:2; txtval = '_byvalidity';
congrus = 1:3; txtcongru = '_bycongrus';
txtcorrect = '_onlycorrect';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% congruTrick is used to make aggregate results from congurency 2 and 3
% together to check for attentional distribution when both probes are on
% the same quadrant. This will create another probability estimates results
% file and save it in the folder "congruTrick"
% this is only useful to generale panels A and B of figure 5 in the article
congruTrick = false;
if congruTrick
    if ~exist('./congruTrick', 'dir')
        mkdir('congruTrick')
    end
    congrus = 1:2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_congrus = size(congrus, 2);

%% Get data from raw files, compute P1 and P2 for all subjects and save
save_data = true;
save_raw_data = true;

P1_all = zeros(13, size(vals, 2), n_congrus, n_obs);
P2_all = zeros(13, size(vals, 2), n_congrus, n_obs);
Delta_all = zeros(13, size(vals, 2), n_congrus, n_obs);
probe_info_all = {}; grat_info_all = {}; delays_all = {};
validity_all = {}; congruency_all = {}; resp_corr_grat_all = {};
fixbreaks_all = {}; respCue_all = {};

for obs_i = 1:n_obs
    obs = observers(obs_i, :);
    fprintf(sprintf('subject %s\n', obs));
    [P1_all(:, :, :, obs_i), P2_all(:, :, :, obs_i),...
        Delta_all(:, :, :, obs_i), probe_info_all{obs_i},...
        grat_info_all{obs_i}, validity_all{obs_i},...
        delays_all{obs_i}, respCue_all{obs_i},...
        congruency_all{obs_i}, resp_corr_grat_all{obs_i}, ...
        fixbreaks_all{obs_i}] = pecar_p_probe_analysis(data_loc, ...
            obs, byvalidity, bycongru, onlycorrect, congruTrick);
end

% if congruTrick (merging congruencies 2 and 3) then only keep this merged
% results from the probability estimations
if congruTrick
    P1_all = P1_all(:, :, 2, :);
    P2_all = P2_all(:, :, 2, :);
    Delta_all = Delta_all(:, :, 2, :);
end

% Save P1, P2 and Discriminant data according to the chosen
% "filters", i.e. by validity, by congruence, only correct
% trials, etc.
if save_data
    if ~congruTrick
        n_obs_data_filename = [save_loc,...
            sprintf('%iobs_P1_P2_Delta%s%s%s',...
            n_obs, txtval, txtcongru, txtcorrect)];
    else
        n_obs_data_filename=[save_loc,...
            sprintf('congruTrick/%iobs_P1_P2_Delta%s%s%s_congruTrick',...
            n_obs, txtval, txtcongru, txtcorrect)];
    end
    save(n_obs_data_filename, 'observers', 'P1_all','P2_all', 'Delta_all');
end

% Save all infos by trial order for all observers
if save_raw_data
    n_obs_data_filename = [save_loc, sprintf(...
        '%iobs_probe_grat_delay_respCue_validity_congruency_respGratCorr_infos', n_obs)];
    save(n_obs_data_filename, 'observers', 'probe_info_all',...
    'grat_info_all', 'delays_all', 'respCue_all', 'validity_all',...
    'congruency_all', 'resp_corr_grat_all');
end


%% Surrogates on FFT computed on individual observers
xfreq = .5:.5:12;
n_freqs = size(xfreq, 2);

% if congruTrick (merging congruencies 2 and 3) then we only left the
% merged probability estimates results, therefore the size of the
% congruency dimension is 1
if congruTrick
    congrus = [1];
end
repeatnumber = 100000; rng('shuffle');
for cong = congrus
    for valind = 1:2
        fft_Pdiff_p_pad = zeros(repeatnumber, n_freqs);
        for repeat = 1:repeatnumber
            if ~mod(repeat, 500); disp(['Repeat number: ' num2str(repeat)]); end

            p1 = zeros(18+n_del+19, n_obs); p2 = zeros(18+n_del+19, n_obs);
            rand_samp = randsample(1:n_del, n_del);
            for obs_i = 1:n_obs
                temp1 = P1_all(rand_samp, valind, cong, obs_i);
                p1(:, obs_i) = [repmat(nanmean(temp1), 1, 18)...
                                temp1' repmat(nanmean(temp1), 1, 19)];
                temp2 = P2_all(rand_samp, valind, cong, obs_i);
                p2(:, obs_i) = [repmat(nanmean(temp2), 1, 18)...
                                temp2' repmat(nanmean(temp2), 1, 19)];
            end
            pdiff = p1 - p2;
            abspdiffpad = abs( fft(pdiff, [], 1) );
            fft_Pdiff_p_pad(repeat, :) = nanmean(abspdiffpad(2:size(xfreq, 2)+1, :), 2);

        end

        % Compute surrogates for Pdiff by shuffling P1 and P2 independently and
        % then re-computing a new Pdiff
        if valind == 2
            txtval = ''; fft_Pdiff_byObs_p_pad_valid = fft_Pdiff_p_pad;
        else txtval = 'in'; fft_Pdiff_byObs_p_pad_invalid = fft_Pdiff_p_pad;
        end
        if ~congruTrick
            save([save_loc, sprintf(['fft_Pdiff_byObs_p_pad_2s_%svalid%s_'...
                'congru%i_%isubjs.mat'], txtval, txtcorrect, cong, n_obs)],...
                sprintf('fft_Pdiff_byObs_p_pad_%svalid', txtval), 'xfreq')
        else
            save([save_loc, sprintf(['congruTrick/fft_Pdiff_byObs_p_pad_2s_%svalid%s_'...
            'congru%i_%isubjs_congruTrick.mat'], txtval, txtcorrect, cong, n_obs)],...
            sprintf('fft_Pdiff_byObs_p_pad_%svalid', txtval), 'xfreq')
        end

    end
end

