function [pboth,pone,pnone,probe_info,validity,delays,grat_info,respCue,...
    congruency] = pecar_probe_analysis(obsdata_loc, file, onlycorrect, probeGratPos)
    %% This program analyzes the probe task for individual stim files
    %
    % Parameters
    %
    % obsdata_loc = directory containing the data file
    % file = '150716_stim01.mat'; (name of data file)
    % onlycorrect = whether to only analyze trials where the response to the
    %               grating task was correct
    % probeGratPos = 'All'/'NoOverlap'/'OneSame'/'BothSame'/'TargetSame'/'DistrSame'
    %                 (only analyses trials in which there was a specific
    %                 probe-grating position relationship)
    %
    % Outputs : All outputs contains data for one session (one raw data file)
    %
    % pboth = whether both probes were reported correctly for each trial
    % pone = whether only one probe was reported correctly for each trial
    % pboth = whether no probe was reported correctly for each trial
    % probe_info = info for probes for each trial
    %              dimensions: n_trials*5*2 = trial_n * (probe identity, nothing,
    %                          probePosX, probePosY, correctResp) * probe 1/2
    % grat_info = info for gratings for each trial = n_trials*3*2
    %             dimensions: n_trials*3*2 = trial_n * (grating x position
    %                         in degrees of visual angle from fixation,
    %                         grating y position, tilt in degrees) * grating 1/2
    % validity = validity for each trial (1 = invalid, 2 = valid);
    % delays = delay of probe presentation for each trial (1 to 13);
    % respCue = response cue for each trial (1 = left or 2 = right)
    % congruency = congruency for each trial (1 = one probe on each hemifield,
    %              2 = both probes on response cue hemifield,
    %              3 = both probes on contra-response cue hemifield);

    %% Load the data
    file_loc = [obsdata_loc file];
    load(strrep(file_loc,'\',filesep))

    %% Get Probe data
    %%% Probe identity
    identity = [1 0;2 0;3 0;4 0;5 0;6 0;7 0;8 0;9 0;10 0;11 0;12 0];
    % Probe positions
    positions = [-16 -13.25 -10.5 -7.75 -5 -2.25 0.5 3.25 6 8.75 11.5 14.25];

    exp = getTaskParameters(myscreen,task);
    expProbe = task{1}.probeTask;
    respCorr=task{1}.randVars.responseCorr;
    if size(respCorr,2)~=exp.nTrials
        respCorr(end:exp.nTrials)=false;
    end
    if onlycorrect; correctness=respCorr==1;
    else correctness=logical(ones(1,exp.nTrials));
    end
    nofixbreaktrials=~logical(task{1,1}.randVars.fixBreak);
    if exp.nTrials~=size(nofixbreaktrials,2)
        nofixbreaktrials(end:exp.nTrials)=false;
    end
    
    %% get response cue, target and gratings & probe positions
    
    cue=exp.randVars.cueside;
    a=logical(cue-1);
    a(~logical(exp.randVars.validity<4))=~a(~logical(exp.randVars.validity<4));
    respCue=a+1;

    ntrials=size(task{1,1}.probeTask.probePresented1,2);
    grat_info=zeros(ntrials,3,2);
    % probe_info dimensions: n_trials,(probe identity, nothing, probePosX,
    % probePosY, correctResp), probe1or2
    probe_info=zeros(ntrials,5,2); probe_info(:,5,:)=NaN;
    for i=1:ntrials
        if isfield(task{1,1},'gratingTask')
            grat_info(i,:,1)=task{1,1}.gratingTask.probePresented1{i};
            grat_info(i,:,2)=task{1,1}.gratingTask.probePresented2{i};
        end
        probe_info(i,1:4,1)=task{1,1}.probeTask.probePresented1{i};
        probe_info(i,1:4,2)=task{1,1}.probeTask.probePresented2{i};
    end
    if size(grat_info,1)~=exp.nTrials
        grat_info(end:exp.nTrials,:,:)=0;
        probe_info(end:exp.nTrials,:,:)=0;
    end
    
    %% Get all information of probe - grating relation to only use specific trials
    % depending on probeGratPos
    
    % for each trial check if 1 of the probes' pos = 1 of the grats' pos
    probeGratSamePosOne=grat_info(:,1,1)==probe_info(:,3,1) |...
        grat_info(:,1,2)==probe_info(:,3,1) |...
        grat_info(:,1,1)==probe_info(:,3,2) |...
        grat_info(:,1,2)==probe_info(:,3,2);
    
    % for each trial none of the probes' pos = any of the grats' pos
    probeGratNoOverlap=grat_info(:,1,1)~=probe_info(:,3,1) &...
        grat_info(:,1,2)~=probe_info(:,3,1) &...
        grat_info(:,1,1)~=probe_info(:,3,2) &...
        grat_info(:,1,2)~=probe_info(:,3,2);
    
    % same thing but if both probes were at the gratings' positions
    probeGratSamePosTwo=grat_info(:,1,1)==probe_info(:,3,1) |...
        grat_info(:,1,2)==probe_info(:,3,1) &...
        grat_info(:,1,1)==probe_info(:,3,2) |...
        grat_info(:,1,2)==probe_info(:,3,2);
    
    % target side
    gratTargPos=zeros(1,exp.nTrials);
    for i=1:exp.nTrials; gratTargPos(i) = grat_info(i,1,respCue(i)); end
    probeTargSamePos=gratTargPos==probe_info(:,3,1)' |...
        gratTargPos==probe_info(:,3,2)';
    
    % distractor side
    distrSide=~(respCue-1)+1; gratDistrPos=zeros(1,exp.nTrials);
    for i=1:exp.nTrials; gratDistrPos(i)=grat_info(i,1,distrSide(i)); end
    probeDistrSamePos=gratDistrPos==probe_info(:,3,1)' |...
        gratDistrPos==probe_info(:,3,2)';

    switch probeGratPos
        case 'All'
            probeGratRelation=true(1, exp.nTrials);
        case 'OneSame'
            probeGratRelation=probeGratSamePosOne';
        case 'BothSame'
            probeGratRelation=probeGratSamePosTwo';
        case 'TargetSame'
            probeGratRelation=probeTargSamePos;
        case 'DistrSame'
            probeGratRelation=probeDistrSamePos;
        case 'NoOverlap'
            probeGratRelation=probeGratNoOverlap';
    end
    
    theTrials = find(nofixbreaktrials & correctness & probeGratRelation);
    grat_info=grat_info(theTrials,:,:);
    probe_info=probe_info(theTrials,:,:);
    respCue=respCue(theTrials);
    ntrials=size(theTrials,2);
    
    validity=(exp.randVars.validity(theTrials)<4)+1;
    congruency=exp.randVars.probeside(theTrials);
    delays=exp.randVars.delays(theTrials);
    
    %% Compare presented and reported probes to get number of correct reports
    pboth = NaN(1,ntrials); pnone = NaN(1,ntrials); pone = NaN(1,ntrials);

    % theTrials contains the trial numbers for the conditions set by the parameters
    index = 1;
    for trial_n = theTrials
        if trial_n <= size(expProbe.probeResponse1,2)
            if ~isempty(expProbe.probeResponse1{trial_n})
                idx1 = find(positions == expProbe.probeResponse1{trial_n});
                idx2 = find(positions == expProbe.probeResponse2{trial_n});
                selected_idx1 = expProbe.list{trial_n}==idx1;
                selected_idx2 = expProbe.list{trial_n}==idx2;
                reported1 = identity(selected_idx1,:);
                reported2 = identity(selected_idx2,:);
                presented1 = expProbe.probePresented1{trial_n};
                presented2 = expProbe.probePresented2{trial_n};

                cor1=(reported1(1)==presented1(1))||(reported2(1)==presented1(1));
                cor2=(reported1(1)==presented2(1))||(reported2(1)==presented2(1));
                probe_info(index,5,:)=[cor1,cor2];
                cor12=cor1+cor2;
                pboth(index)=cor12==2; pone(index)=cor12==1; pnone(index)=cor12==0;
            end
        end
        index = index + 1;
    end
end

