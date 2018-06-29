function [P1, P2, Delta, probe_info, grat_info, validity, delays,...
    respCue, congruency] = pecar_p_probe_analysis(data_loc, obs)
    %
    % Parameters
    %
    % data_loc = string, location of the raw observer's data
    % obs = e.g. 'ax'; (observer's initials)
    %
    % Outputs
    %
    % P1 = 2D (or 3D) P1 values for each delay*validity(*congruency)
    % P2 = 2D (or 3D) P2 values for each delay*validity(*congruency)
    % probe_info = info for probes for each trial
    %              dimensions: n_trials*5*2 = trial_n * (probe identity, nothing,
    %                          probePosX, probePosY, correctResp) * probe 1/2
    % grat_info = info for gratings for each trial
    %             dimensions: n_trials*3*2 = trial_n * (grating x position
    %                         in degrees of visual angle from fixation,
    %                         grating y position, tilt in degrees) * grating 1/2
    % validity = validity for each trial (1 = invalid, 2 = valid);
    % delays = delay of probe presentation for each trial (1 to 13);
    % respCue = response cue for each trial (1 = left or 2 = right)
    % congruency = congruency for each trial (1 = one probe on each hemifield,
    %              2 = both probes on response cue hemifield,
    %              3 = both probes on contra-response cue hemifield);
        

    %% Obtain pboth, pone and pnone for each run and concatenate over sessions
    probe_info = []; validity = []; delays = []; grat_info = [];
    respCue = []; pboth_all = []; pone_all = []; pnone_all = [];
    congruency = [];

    obsdata_loc=[data_loc 'subj_' obs '/'];
    files = dir(strrep(obsdata_loc,'/', filesep));
    for file_n = 1:size(files, 1)
        filename = files(file_n).name;
        fileL = size(filename,2);
        if fileL == 17 && strcmp(filename(fileL-4+1:fileL),'.mat') &&...
                isa(str2double(filename(1:6)),'double')
            [pboth, pone, pnone, probe_info_block, validity_block,...
                delays_block, grat_info_block, respCue_block,...
                congruency_block] = pecar_probe_analysis(obsdata_loc, filename);
            probe_info = cat(1, probe_info, probe_info_block);
            validity = cat(2, validity, validity_block);
            delays = cat(2, delays, delays_block);
            grat_info = cat(1, grat_info, grat_info_block);
            respCue = cat(2, respCue, respCue_block);
            congruency = cat(2, congruency, congruency_block);
            pboth_all = cat(2, pboth_all, pboth);
            pone_all = cat(2, pone_all, pone);
            pnone_all = cat(2, pnone_all, pnone);
        end
    end

    %% Make 2D matrix with delay*validity and computes P1 and P2
    val_uniq = [2, 1]; n_vals = 2;

    pboth_2D = zeros(13, n_vals); pbothstd_2D = zeros(13, n_vals);
    pnone_2D = zeros(13, n_vals); pnonestd_2D = zeros(13, n_vals);
    pone_2D = zeros(13, n_vals); ponestd_2D = zeros(13, n_vals);
    for delay = unique(delays)
        for val = val_uniq
            mask = delays == delay & validity == val;

            probe_data = sum(probe_info(mask, 5, :), 3);
            pboth_2D(delay, val) = mean(probe_data == 2);
            pbothstd_2D(delay, val) = std(probe_data == 2) / sqrt(sum(mask));
            pone_2D(delay, val) = mean(probe_data == 1);
            ponestd_2D(delay, val) = std(probe_data == 1) / sqrt(sum(mask));
            pnone_2D(delay, val) = mean(probe_data == 0);
            pnonestd_2D(delay, val) = std(probe_data == 0) / sqrt(sum(mask));
        end
    end
    [P1, P2, Delta] = quadratic_analysis(pboth_2D, pnone_2D);
    
end
