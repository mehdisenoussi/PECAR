%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results on the 2AFC orientation discrimination task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pecar_loc = './'; % adapt to location of the scripts and data
addpath(genpath([pecar_loc, 'mgl-master/']))
data_loc = [pecar_loc, 'data_pecar/'];

% Load data
observers = ['ym'; 'ac'; 'al'; 'sa'; 'el'; 'gm'; 'hs'; 'hw'; 'js'; 'ma'; 'nv'];
n_obs = size(observers, 1);
binstep = .05; bins = 0:binstep:2; n_hist_vals = size(bins, 2) - 1;
respTime_allsubjs = zeros(n_obs, 2, n_hist_vals);
dprime_allsubjs = zeros(n_obs, 2);
accuracy_allsubjs = zeros(n_obs, 2);
respTime_median = zeros(n_obs, 2);

for obs_i = 1:size(observers, 1)
    obs = observers(obs_i, :);
    
    subjdata_loc = [data_loc, 'subj_', obs, '/'];
    files = dir(strrep(subjdata_loc, '\', filesep));
    run_n = 0;
    for file_n = 1:size(files, 1)
        filename = files(file_n).name;
        fileL = size(filename, 2);
        if fileL == 17 && strcmp(filename(fileL-4+1:fileL), '.mat') &&...
                isa(str2double(filename(1:6)), 'double')
            load(strrep([subjdata_loc, filename], '\', filesep))
            exp = getTaskParameters(myscreen, task);
            run_n = run_n + 1;

            nofixbreaktrials = ~logical(task{1, 1}.randVars.fixBreak);
            if exp.nTrials ~= size(nofixbreaktrials, 2)
                nofixbreaktrials(end:exp.nTrials) = false;
            end
            rt_mask = nofixbreaktrials | ~logical(isnan(exp.responseVolume));
            respTime = exp.reactionTime(nofixbreaktrials);
            validity = exp.randVars.validity(nofixbreaktrials) < 4;
            respCorr = task{1,1}.randVars.responseCorr;
            if size(respCorr, 2) == exp.nTrials
                respCorr = task{1, 1}.randVars.responseCorr(nofixbreaktrials);
            else
                respCorr(end:exp.nTrials) = false;
                respCorr = respCorr(nofixbreaktrials);
            end
            cue = exp.randVars.cueside(nofixbreaktrials);
            resp = exp.response(nofixbreaktrials);
            ngoodtrials = sum(nofixbreaktrials);
            a = logical(cue - 1); a(~logical(validity)) = ~a(~logical(validity));
            respCue = a+1;

            stimOris = [exp.randVars.stimOrientation1(1, nofixbreaktrials);...
                      exp.randVars.stimOrientation2(1, nofixbreaktrials)];

            targetOris = zeros(1, ngoodtrials);
            for j = 1:ngoodtrials
                targetOris(j) = stimOris(respCue(j) ,j);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % compute dprime
            hits = zeros(2, ngoodtrials);
            falsealarms = zeros(2, ngoodtrials); n_sig = [0, 0];
            n_nsig = [0, 0];
            valid_trials_ind = find(validity == 1);
            invalid_trials = find(validity == 0);
            trialsbyvalidity = {valid_trials_ind, invalid_trials};
            for valid_ind = 1:2
                for j = trialsbyvalidity{valid_ind}
                    if targetOris( j) == 1
                        n_sig(valid_ind) = n_sig(valid_ind) + 1;
                        if resp(j) == 1; hits(valid_ind, j) = 1; end
                    else
                        n_nsig(valid_ind) = n_nsig(valid_ind) + 1;
                        if resp(j) == 1; falsealarms(valid_ind, j) = 1; end
                    end
                end
            end

            if (sum(hits(1, :)) / n_sig(1)) ~= 0 && ((sum(hits(1, :)) / n_sig(1))) ~= 1
                valhits = sum(hits(1, :)) / n_sig(1);
            elseif (sum(hits(1, :)) / n_sig(1)) == 0
                valhits = (sum(hits(1, :)) + 1) / n_sig(1);
            elseif ((sum(hits(1, :)) / n_sig(1))) == 1
                valhits = (sum(hits(1, :)) - 1) / n_sig(1);
            end

            if (sum(hits(2, :)) / n_sig(2)) ~= 0 && (sum(hits(2, :)) / n_sig(2)) ~= 1
                invalhits = sum(hits(2, :)) / n_sig(2);
            elseif (sum(hits(2, :)) / n_sig(2)) == 0
                invalhits = (sum(hits(2, :)) + 1) / n_sig(2);
            elseif (sum(hits(2, :)) / n_sig(2)) == 1
                invalhits = (sum(hits(2, :)) - 1) / n_sig(2); 
            end

            if ((sum(falsealarms(1, :)) / n_nsig(1))) ~= 0 &&...
                    ((sum(falsealarms(1, :)) / n_nsig(1))) ~= 1
                valfas = sum(falsealarms(1, :)) / n_nsig(1);
            elseif ((sum(falsealarms(1, :)) / n_nsig(1))) == 0
                valfas = (sum(falsealarms(1, :)) + 1) / n_nsig(1);
            elseif ((sum(falsealarms(1, :)) / n_nsig(1))) == 1
                valfas = (sum(falsealarms(1, :)) - 1) / n_nsig(1);
            end

            if (sum(falsealarms(2, :)) / n_nsig(2)) ~= 0 &&...
                    (sum(falsealarms(2, :)) / n_nsig(2)) ~= 1
                invalfas = sum(falsealarms(2, :)) / n_nsig(2);
            elseif (sum(falsealarms(2, :)) / n_nsig(2)) == 0
                invalfas = (sum(falsealarms(2, :)) + 1) / n_nsig(2);
            elseif (sum(falsealarms(2, :)) / n_nsig(2)) == 1
                invalfas = (sum(falsealarms(2, :)) - 1) / n_nsig(2);
            end

            dprime_valid = norminv(valhits) - norminv(valfas);
            dprime_invalid=norminv(invalhits)-norminv(invalfas);


            respCorr_byblock{run_n} = respCorr;
            respTime_byblock{run_n} = respTime;
            validity_byblock{run_n} = validity;
            dprimevalid_byblock{run_n} = dprime_valid;
            dprimeinvalid_byblock{run_n} = dprime_invalid;
        end
    end
    n_run = run_n;

    %% Compute across run values
    for run_n = 1:(n_run - 1)
        respTime_all = [respTime_byblock{run_n}, respTime_byblock{run_n + 1}];
        validity_all = [validity_byblock{run_n}, validity_byblock{run_n + 1}];
    end
    respTime_allsubjs(obs_i, 1, :) = histcounts(...
        respTime_all(logical(validity_all)), bins, 'Normalization', 'probability');
    respTime_allsubjs(obs_i, 2, :) = histcounts(...
        respTime_all(~logical(validity_all)), bins, 'Normalization', 'probability');
    respTime_median(obs_i, :) =...
        [nanmedian(respTime_all(logical(validity_all))),...
        nanmedian(respTime_all(~logical(validity_all)))];
    dprimevalid_all = zeros(1, 4); dprimeinvalid_all = zeros(1, 4);
    respCorr_all = zeros(2, 4); validity_all = zeros(2, 4);
    for run_n = 1:n_run
        dprimevalid_all(run_n) = dprimevalid_byblock{run_n};
        dprimeinvalid_all(run_n) = dprimeinvalid_byblock{run_n};
        respCorr_all(:, run_n) = [mean(respCorr_byblock{run_n}(validity_byblock{run_n})),...
            mean(respCorr_byblock{run_n}(~validity_byblock{run_n}))];
    end
    dprime_allsubjs(obs_i, :) = [mean(dprimevalid_all), mean(dprimeinvalid_all)];
    accuracy_byvalid = [mean(respCorr_all(1, :)), mean(respCorr_all(2, :))];
    accuracy_allsubjs(obs_i, :) = accuracy_byvalid;
end

%% Saving the file

save_loc = [pecar_loc, 'results/'];

n_obs_data_filename = [save_loc,...
    sprintf('%iobs_respTime_dprime_acc_byvalidity', n_obs)];
if ~exist(n_obs_data_filename, 'file')
    save(n_obs_data_filename, 'observers', 'respTime_allsubjs',...
        'respTime_median', 'dprime_allsubjs', 'accuracy_allsubjs');
end

%%
n_obs_data_filename = [save_loc,...
    sprintf('%iobs_respTime_dprime_acc_byvalidity', n_obs)];
load(n_obs_data_filename)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot reaction times
figure; subplot(1, 2, 2); hold on;

% compute SEM
SEM = std(respTime_median(:, 1), [], 1) / sqrt(n_obs);
% compute 95% confidence interval
perc = .05;
ts = tinv([perc/2., 1 - (perc/2.)], n_obs - 1);
CI = mean(respTime_median(:, 1), 1) + ts * SEM;

meantoplot = mean(respTime_median(:, 1), 1);
plot(1, meantoplot, 'ko', 'MarkerEdgeColor', [0, 0, 0],...
        'MarkerFaceColor', [0, 0, 0], 'LineWidth', 1.5, 'MarkerSize', 10);
plot([1, 1], CI, 'k', 'LineWidth', 2)

% compute SEM
SEM = std(respTime_median(:, 2), [], 1) / sqrt(n_obs);
% compute 95% confidence interval
perc = .05;
ts = tinv([perc / 2., 1 - (perc / 2.)], n_obs - 1);
CI = mean(respTime_median(:, 2), 1) + ts * SEM;

meantoplot = mean(respTime_median(:, 2), 1);
plot(2, meantoplot, 'ko', 'MarkerEdgeColor', [0, 0, 0],...
        'MarkerFaceColor', [0, 0, 0], 'LineWidth', 1.5, 'MarkerSize', 10);
plot([2, 2], CI, 'k', 'LineWidth', 2)

xlabel('Conditions'); ylabel('Average median'); grid on;
title(sprintf('Average median reaction times\nvalid versus invalid condition'));
xlim([0, 3]);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d-prime by validity
subplot(1, 2, 1); hold on; 

% compute SEM
SEM = std(dprime_allsubjs(:,1) ,[], 1) / sqrt(n_obs);
% compute 95% confidence interval
perc = .05;
ts = tinv([perc / 2., 1 - (perc / 2.)], n_obs - 1);
CI = mean(dprime_allsubjs(:, 1), 1) + ts * SEM;

meantoplot = mean(dprime_allsubjs(:, 1), 1);
plot(1, meantoplot, 'ko', 'MarkerEdgeColor', [0, 0, 0],...
        'MarkerFaceColor', [0, 0, 0], 'LineWidth', 1.5, 'MarkerSize', 10);
plot([1, 1], CI, 'k', 'LineWidth', 2)
scatter(ones(1, n_obs) + .2, dprime_allsubjs(:, 1), circ_sz,...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [.5, .5, .5], 'LineWidth', 1.5)


% compute SEM
SEM = std(dprime_allsubjs(:, 2) , [], 1) / sqrt(n_obs);
% compute 95% confidence interval
perc = .05;
ts = tinv([perc / 2., 1 - (perc / 2.)], n_obs - 1);
CI = mean(dprime_allsubjs(:, 2), 1) + ts * SEM;

meantoplot = mean(dprime_allsubjs(:, 2), 1);
plot(2, meantoplot, 'ko', 'MarkerEdgeColor', [0, 0, 0],...
        'MarkerFaceColor', [0, 0, 0], 'LineWidth', 1.5, 'MarkerSize', 10);
plot([2, 2], CI, 'k', 'LineWidth', 2)
scatter(ones(1, n_obs) + .8, dprime_allsubjs(:, 2), circ_sz,...
    'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [.5, .5, .5], 'LineWidth', 1.5)

    
for obs_i = 1:n_obs; plot([1.2, 1.8], [mean(dprime_allsubjs(obs_i, 1), 1),...
        mean(dprime_allsubjs(obs_i, 2), 1)], 'k-'); end

title('D-prime by validity');
ylabel('d-prime'); grid on;
set(gca, 'XTick', [1, 2]); set(gca, 'XTickLabel', {'Valid', 'Invalid'});
set(gca, 'XLim', [.8, 2.2]);


