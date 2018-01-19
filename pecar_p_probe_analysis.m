function [P1,P2,Delta,probe_info,grat_info,validity,delays,respCue, congruency] =...
    pecar_p_probe_analysis(data_loc, obs,printFg,byvalidity,bycongru,onlycorrect,probeGratPos)
    %
    % Parameters
    %
    % obs = e.g. 'ax'; (observer's initials)
    % printFg = true/false; (figures are printed and saved for each observer)
    % byvalidity = true/false; (separate analysis for valid and invalid trials)
    % bycongru = true/false; (separate analysis by trial congruencies condition)
    % onlycorrect = true/false; (analyze only trials with correct response to grating task)
    % probeGratPos = 'All'/'NoOverlap'/'OneSame'/'BothSame'/'TargetSame'/'DistrSame'
    %                 (only analyses trials in which there was a specific
    %                 probe-grating position relationship)
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
    pboth_all=[]; pone_all=[]; pnone_all=[]; congruency=[];

    obsdata_loc=[data_loc 'subj_' obs '/'];
    files = dir(strrep(obsdata_loc,'/',filesep));
    for file_n = 1:size(files,1)
        filename = files(file_n).name;
        fileL = size(filename,2);
        if fileL == 17 && strcmp(filename(fileL-4+1:fileL),'.mat') && isa(str2double(filename(1:6)),'double')
            [pboth,pone,pnone,probe_info_block,validity_block,delays_block,...
                grat_info_block, respCue_block, congruency_block]=...
                pecar_probe_analysis(obsdata_loc, filename, onlycorrect, probeGratPos);
            probe_info=cat(1,probe_info,probe_info_block);
            validity=cat(2,validity,validity_block);
            delays=cat(2,delays,delays_block);
            grat_info=cat(1,grat_info,grat_info_block);
            respCue=cat(2,respCue,respCue_block);
            congruency=cat(2,congruency,congruency_block);
            pboth_all=cat(2,pboth_all,pboth);
            pone_all=cat(2,pone_all,pone);
            pnone_all=cat(2,pnone_all,pnone);
        end
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


    %% if printFg print individual observers' 
    if printFg
        saveFileLoc = '/figures/probetask/';
        delays_uniq=linspace(40,520,13);
        if bycongru; plot_order=[1 3 5 2 4 6]; else plot_order=[1 2]; end
        
        % Plot PBoth, POne and PNone
        h1=figure('Position',get(groot,'ScreenSize'));plotn=1;        
        for val=val_uniq
            for cong=cong_uniq
                subplot(size(cong_uniq,2),size(val_uniq,2),plot_order(plotn)); hold on;
                errorbar(delays_uniq,pboth_ND(:,val,cong),pbothstd_ND(:,val,cong),...
                    'ro-','LineWidth',2,'MarkerFaceColor',[1 1 1],'MarkerSize',8,'Color',[1 0 0])
                errorbar(delays_uniq,pone_ND(:,val,cong),ponestd_ND(:,val,cong),...
                    'ro-','LineWidth',2,'MarkerFaceColor',[1 1 1],'MarkerSize',8,'Color',[0 1 0])
                errorbar(delays_uniq,pnone_ND(:,val,cong),pnonestd_ND(:,val,cong),...
                    'ro-','LineWidth',2,'MarkerFaceColor',[1 1 1],'MarkerSize',8,'Color',[0 0 1])
                legend('PBoth','POne','PNone')

                set(gca,'YTick',0:.2:1,'FontSize',18,'LineWidth',2','Fontname','Ariel')

                ylabel('Percent correct','FontSize',20,'Fontname','Ariel')
                xlabel('Time from search array onset [ms]','FontSize',20,'Fontname','Ariel')
                if val==2; title('Valid'); elseif val==1; title('Invalid'); end
                ylim([-.05 1.05]); xlim([0 560]); grid on; plotn=plotn+1;
                axis square
            end
        end
        suptitle(sprintf('Raw proba - obs: (%s)\t -',obs));
        namefig=[data_loc saveFileLoc 'probetask_rawProbs_' obs];
        set(gcf,'PaperPositionMode','auto')
        print (namefig, '-djpeg', '-r0');
        close


        %% Plot p1 and p2 for each probe delay
        h2=figure('Position',get(groot,'ScreenSize')); plotn=1;
        for val=val_uniq
            for cong=cong_uniq
                subplot(size(cong_uniq,2),size(val_uniq,2),plot_order(plotn)); hold on;
                plot(delays_uniq,P1(:,val,cong),'ro-','LineWidth',3,...
                    'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',[.96 .37 .15])
                plot(delays_uniq,P2(:,val,cong),'go-','LineWidth',3,...
                    'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',[.13 .7 .15])

                legend('p1','p2','Location','SouthEast')
                set(gca,'YTick',0:.2:1.2,'FontSize',18,'LineWidth',2','Fontname','Ariel')
                set(gca,'XTick',0:100:500,'FontSize',18,'LineWidth',2','Fontname','Ariel')
                if val==2; title('Valid'); elseif val==1; title('Invalid'); end
                ylabel('Probe report probabilities','FontSize',20,'Fontname','Ariel')
                xlabel('Time from search task onset [ms]','FontSize',20,'Fontname','Ariel')
                ylim([-.1 .6]); xlim([0 560]); grid on; plotn=plotn+1;
                axis square
            end
        end
        suptitle(sprintf('P1 & P2 - obs: (%s)\t -',obs));
        namefig=[data_loc saveFileLoc 'p1p2_' obs];
        set(gcf,'PaperPositionMode','auto')
        print (namefig, '-djpeg', '-r0');
        close

        %% Plot p1-p2 for each probe delay
        h3=figure('Position',get(groot,'ScreenSize')); hold on; plotn=1;
        for val=val_uniq
            for cong=cong_uniq
                subplot(size(cong_uniq,2),size(val_uniq,2),plot_order(plotn)); hold on;
                plot(delays_uniq,P1(:,val,cong)-P2(val,cong),'ro-','LineWidth',3,...
                    'MarkerFaceColor',[1 1 1],'MarkerSize',10,'Color',[0.5450 0.5350 0.1500])

                set(gca,'YTick',-.8:.2:.8,'FontSize',14,'LineWidth',2','Fontname','Ariel')
                set(gca,'XTick',0:100:500,'FontSize',14,'LineWidth',2','Fontname','Ariel')
                if val==2; title('Valid'); elseif val==1; title('Invalid'); end
                ylabel('Probe report probabilities','FontSize',12,'Fontname','Ariel')
                xlabel('Time from search task onset [ms]','FontSize',12,'Fontname','Ariel')
                ylim([-.5 .5]); xlim([0 560]); grid on; plotn=plotn+1;
                axis square
            end
        end
        suptitle(sprintf('P1 - P2 - obs: (%s)\t -',obs));
        namefig=[data_loc saveFileLoc 'diff_p1p2_' obs];
        set(gcf,'PaperPositionMode','auto')
        print(namefig, '-djpeg', '-r0');
        close
    end    

