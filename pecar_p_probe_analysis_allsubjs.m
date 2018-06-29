pecar_loc='/Users/mehdisenoussi/Dropbox/postphd/lpp/pecar/';
addpath(genpath([pecar_loc 'soft/mgl-master/']))
data_loc = '/Users/mehdisenoussi/pecar/pecar_data/';
save_loc = './processed_data/';
delays = 40:40:520;

observers = ['ym'; 'ac'; 'al'; 'sa'; 'el'; 'gm'; 'hs'; 'hw'; 'js'; 'ma'; 'nv'];
n_obs = size(observers, 1);

%% Get data from raw files, compute P1 and P2 for all subjects and save

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

%% make bootstrap for P1 minus P2 FFT amplitude
tic
repeatnumber = 100000;

n_obs_data_filename=[save_loc sprintf('%iobs_P1_P2_Delta', n_obs)];
load(n_obs_data_filename);
for val = vals
    pecar_p1p2_bootstrap(save_loc, P1_all, P2_all, repeatnumber, val);
end
toc











