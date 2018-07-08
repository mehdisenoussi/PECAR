pecar_loc = './'; % adapt to location of the scripts and data
addpath(genpath([pecar_loc, 'mgl-master/']))
data_loc = './data_pecar/';
save_loc = './processed_data/';

if ~exist(save_loc); mkdir(save_loc); end
    
delays = 40:40:520;

observers = ['ym'; 'ac'; 'al'; 'sa'; 'el'; 'gm'; 'hs'; 'hw'; 'js'; 'ma'; 'nv'];
n_obs = size(observers, 1);

%% Get data from raw files, compute P1, P2 and Discriminant for all observers and save

save_data = true;
save_raw_data = true;

vals = 1:2;
            
P1_all = zeros(13, size(vals, 2), n_obs);
P2_all = zeros(13, size(vals, 2), n_obs);
Delta_all = zeros(13, size(vals, 2), n_obs);
probe_info_all = {}; grat_info_all = {}; delays_all = {}; respCue_all = {};
validity_all = {}; congruency_all = {};

for obs_i = 1:n_obs
    obs = observers(obs_i, :);
    fprintf(sprintf('subject %s\n', obs));
    [P1_all(:, :, obs_i), P2_all(:, :, obs_i), Delta_all(:, :, obs_i),...
        probe_info_all{obs_i}, grat_info_all{obs_i}, validity_all{obs_i},...
        delays_all{obs_i}, respCue_all{obs_i}, congruency_all{obs_i}] =...
        pecar_p_probe_analysis(data_loc, obs);
end

% Save P1, P2 and Discriminant
if save_data
    n_obs_data_filename = [save_loc sprintf('%iobs_P1_P2_Delta', n_obs)];
    save(n_obs_data_filename, 'observers', 'P1_all','P2_all', 'Delta_all');
end

% Save all infos by trial order for all observers
if save_raw_data
    n_obs_data_filename = [save_loc sprintf(...
        '%iobs_probe_grat_delay_respCue_validity_congruency_infos', n_obs)];
    save(n_obs_data_filename, 'observers', 'probe_info_all',...
    'grat_info_all', 'delays_all', 'respCue_all', 'validity_all',...
    'congruency_all');
end   

%% Frequency amplitude of Pdiff + surrogates

Pdiff = P1_all - P2_all;
Pdiffavg = mean(Pdiff, 3);

% pad grand average Pdiff for each validity
avgpdiff1 = squeeze(mean(Pdiffavg, 1));
Pdiffavgpad = [];
for valind = 1:2
        % pad the grand average Pdiff with its average across delays to
        % have a 2 seconds signal
        Pdiffavgpad(:, valind) = ...
            [repmat(avgpdiff1(valind), 1, 18), Pdiffavg(:, valind)',...
             repmat(avgpdiff1(valind), 1, 19)];
end

fft_Pdiffavgpad = fft(Pdiffavgpad, 50, 1); % do the FFT 

xfreq = .5:.5:12; n_freqs = size(xfreq, 2); % the sampled frequencies
% discard the symmetrical part of the FFT
a_fft_Pdiffavgpad=abs(fft_Pdiffavgpad(2:size(xfreq, 2) + 1, :, :));


% Permutation to create surrogate frequency profiles
repeatnumber = 100000; ndel = 13;
for valind = 1:2
    % Compute surrogates for Pdiff by shuffling P1 and P2 independently and
    % then re-computing a new Pdiff
    fft_Pdiff1_ALL_p_pad = zeros(repeatnumber, n_freqs);
    fft_Pdiff2_ALL_p_pad = zeros(repeatnumber, n_freqs);
    for repeat = 1:repeatnumber
        if ~mod(repeat, 500); disp(['Repeat number: ' num2str(repeat)]); end

        p1 = zeros(ndel, n_obs); p2 = zeros(ndel, n_obs);

        for obs_i = 1:n_obs
            p1(:, obs_i) = P1_all(randsample(1:ndel, ndel), valind, obs_i);
            p2(:, obs_i) = P2_all(randsample(1:ndel, ndel), valind, obs_i);
        end        
        pdiff1 = p1-p2; pdiff1avg = mean(pdiff1, 2); avgpdiff1 = mean(pdiff1avg);
        pdiffavgpad1 = [repmat(avgpdiff1, 1, 18) pdiff1avg'...
                        repmat(avgpdiff1, 1, 19)];
        abspdiffavgpad1 = abs( fft(pdiffavgpad1) );
        fft_Pdiff1_ALL_p_pad(repeat, :) = abspdiffavgpad1(2:size(xfreq, 2) + 1);
    end

    if valind == 2
        txtval = ''; fft_Pdiff_ALL_p_pad_valid = fft_Pdiff1_ALL_p_pad;
    else txtval = 'in'; fft_Pdiff_ALL_p_pad_invalid = fft_Pdiff1_ALL_p_pad;
    end
    save([save_loc, sprintf('fft_Pdiff_ALL_p_pad_2s_%svalid_%isubjs.mat',...
     txtval, n_obs)], sprintf('fft_Pdiff_ALL_p_pad_%svalid', txtval), 'xfreq')
end









