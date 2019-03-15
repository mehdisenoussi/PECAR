function [P1, P2, Delta, probe_info, grat_info, validity, delays, respCue,...
    congruency, resp_corr, fixbreaks_all] = pecar_p_probe_analysis(...
        data_loc, obs, byvalidity, bycongru, onlycorrect, congruTrick)
    %
    % Parameters
    %
    % obs = e.g. 'ax'; (observer's initials)
    % byvalidity = true/false; (separate analysis for valid and invalid trials)
    % bycongru = true/false; (separate analysis by trial congruencies condition)
    % onlycorrect = true/false; (analyze only trials with correct response to grating task)
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
    probe_info=[]; validity=[]; delays=[];grat_info=[]; respCue=[];
    pboth_all=[]; pone_all=[]; pnone_all=[]; congruency=[]; resp_corr=[];
    fixbreaks_all=[];

    obsdata_loc=[data_loc 'subj_' obs '/'];
    files = dir(strrep(obsdata_loc,'/',filesep));
    for file_n = 1:size(files,1)
        filename = files(file_n).name;
        fileL = size(filename,2);
        if fileL == 17 && strcmp(filename(fileL-4+1:fileL),'.mat') && isa(str2double(filename(1:6)),'double')
            [pboth,pone,pnone,probe_info_block,validity_block,delays_block,...
                grat_info_block, respCue_block, congruency_block,...
                resp_corr_block, fixbreaks_block] = pecar_probe_analysis(...
                    obsdata_loc, filename, onlycorrect);
            probe_info=cat(1,probe_info,probe_info_block);
            validity=cat(2,validity,validity_block);
            delays=cat(2,delays,delays_block);
            grat_info=cat(1,grat_info,grat_info_block);
            respCue=cat(2,respCue,respCue_block);
            congruency=cat(2,congruency,congruency_block);
            resp_corr=cat(2,resp_corr,resp_corr_block);
            pboth_all=cat(2,pboth_all,pboth);
            pone_all=cat(2,pone_all,pone);
            pnone_all=cat(2,pnone_all,pnone);
            fixbreaks_all=cat(1,fixbreaks_all,fixbreaks_block);
        end
    end
    
    if congruTrick
        congruency(congruency == 3) = 2;
    end
    
    %% Make 2D (or 3D) matrix with delay*validity(*congruency) and computes P1 and P2
    if bycongru; cong_uniq=unique(congruency); n_congs=size(cong_uniq,2);
    else cong_uniq=[1]; n_congs=1; end
    if byvalidity; val_uniq=[2 1]; n_vals=size(val_uniq,2);
    else val_uniq=[1]; n_vals=1; end

    pboth_ND=zeros(13,n_vals,n_congs); pbothstd_ND=zeros(13,n_vals,n_congs);
    pnone_ND=zeros(13,n_vals,n_congs); pnonestd_ND=zeros(13,n_vals,n_congs);
    pone_ND=zeros(13,n_vals,n_congs); ponestd_ND=zeros(13,n_vals,n_congs);
    for delay = unique(delays)
        for val = val_uniq
            for cong = cong_uniq
                if bycongru && byvalidity; mask=delays==delay & validity==val & congruency==cong;
                elseif bycongru; mask=delays==delay & congruency==cong;
                elseif byvalidity; mask=delays==delay & validity==val;
                else mask=delays==delay; end
                probe_data=sum(probe_info(mask,5,:),3);
                pboth_ND(delay, val, cong)=mean(probe_data==2);
                pbothstd_ND(delay, val, cong)=std(probe_data==2)/sqrt(sum(mask));
                pone_ND(delay, val, cong)=mean(probe_data==1);
                ponestd_ND(delay, val, cong)=std(probe_data==1)/sqrt(sum(mask));
                pnone_ND(delay, val, cong)=mean(probe_data==0);
                pnonestd_ND(delay, val, cong)=std(probe_data==0)/sqrt(sum(mask));
            end
        end
    end
    [P1, P2, Delta] = quadratic_analysis(pboth_ND, pnone_ND);

end    

