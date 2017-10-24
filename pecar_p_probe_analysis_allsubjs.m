pecar_loc='/Users/mehdisenoussi/Dropbox/postphd/laura/pecar/';
addpath(genpath([pecar_loc 'soft/mrTools-master']))
addpath(genpath([pecar_loc 'soft/mgl-master/']))
data_loc='/Volumes/PNY/PECAR/data/';
delays=40:40:520;

%observers=['cs'; 'sr'; 'ym'; 'ac'; 'al'; 'sa'; 'el'; 'gm'; 'hs'; 'hw'; 'js'; 'ma'; 'nv'];
observers=['ym'; 'ac'; 'al'; 'sa'; 'el'; 'gm'; 'hs'; 'hw'; 'js'; 'ma'; 'nv'];
n_obs=size(observers,1);

probeGratPos='All'; % 'NoOverlap' 'OneSame'  'BothSame'  'TargetSame'  'DistrSame'
byvalidity=true; bycongru=false; onlycorrect=false;

txtval=''; txtcongru=''; txtcorrect='';
if byvalidity; vals=1:2; txtval='_byvalidity'; else vals=[1]; end
if bycongru; congrus=1:3; txtcongru='_bycongrus'; else congrus=[1]; end
if onlycorrect; txtcorrect='_onlycorrect'; end

%% Get data from raw files, compute P1 and P2 for all subjects and save

save_data=false;
for byvalidity=[true]
    for bycongru=[false]
        for onlycorrect=[false]
            txtval=''; txtcongru=''; txtcorrect='';
            if byvalidity; vals=1:2; txtval='_byvalidity'; else vals=[1]; end
            if bycongru; congrus=1:3; txtcongru='_bycongrus'; else congrus=[1]; end
            if onlycorrect; txtcorrect='_onlycorrect';end
            
            P1_all=zeros(13,size(vals,2),size(congrus,2),n_obs);
            P2_all=zeros(13,size(vals,2),size(congrus,2),n_obs);
            probe_info_all={}; grat_info_all={}; delays_all={}; respCue_all={};
            validity_all={}; congruency_all={};

            for obs_i=1:n_obs
                obs=observers(obs_i,:);
                fprintf(sprintf('subject %s\n',obs));
                [P1_all(:,:,:,obs_i), P2_all(:,:,:,obs_i), probe_info_all{obs_i},...
                    grat_info_all{obs_i}, validity_all{obs_i},delays_all{obs_i},...
                    respCue_all{obs_i}, congruency_all{obs_i}] = pecar_p_probe_analysis(...
                    data_loc,obs,true,byvalidity,bycongru,onlycorrect,probeGratPos);
            end
            if save_data
                n_obs_data_filename=[data_loc sprintf('%iobs_P1_P2%s%s%s_probeGratPos_%s',...
                    n_obs,txtval,txtcongru,txtcorrect,probeGratPos)];
                save(n_obs_data_filename, 'observers', 'P1_all','P2_all');
            end
        end
    end
end

%% make bootstrap
tic
repeatnumber=100000;
txtval='_byvalidity';
for bycongru=[false]
    n_obs_data_filename=[data_loc sprintf('%iobs_P1_P2%s%s%s_probeGratPos_%s',...
        n_obs,txtval,txtcongru,txtcorrect,probeGratPos)];
    load(n_obs_data_filename);
    for val=vals
        for cong=congrus
            pecar_p1p2_bootstrap(squeeze(P1_all(:,val,cong,:)),...
                squeeze(P2_all(:,val,cong,:)), repeatnumber, val, probeGratPos);
        end
    end
end
toc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P1 and P2
bycongru=false; byvalidity=true; onlycorrect=false;
if bycongru; plot_order=[1 3 5 2 4 6]; congrus=1:3; txtcongru='_bycongrus';
else plot_order=[1 2 3]; congrus=[1]; txtcongru=''; end
if byvalidity; vals=2:-1:1; txtval='_byvalidity'; else vals=[1]; end
if onlycorrect; txtcorrect='_onlycorrect'; end

% n_obs_data_filename=[data_loc sprintf('%iobs_P1_P2%s%s%s_probeGratPos_%s',...
%         n_obs,txtval,txtcongru,txtcorrect,probeGratPos)];
% load(n_obs_data_filename)
        
figure('Position',get(groot,'ScreenSize'));
plotn=1;
for val=vals
    for cong=congrus
        subplot(size(congrus,2),size(vals,2),plot_order(plotn)); plotn=plotn+1; hold on;
        if val==1; cols=[230, 60, 23; 230, 138, 23]/256.;
        else cols=[51, 78, 198; 51, 180, 255]/256.; end
        errorbar(delays,mean(P1_all(:,val,cong,:),4), std(P1_all(:,val,cong,:),[],4)./sqrt(n_obs),...
            'o-','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',cols(1,:))
        errorbar(delays,mean(P2_all(:,val,cong,:),4), std(P2_all(:,val,cong,:),[],4)./sqrt(n_obs),...
            'o-','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',cols(2,:))
        if val==2; title('Valid'); elseif val==1; title('Invalid'); end
        legend('p1','p2','Location','NorthEast')

        set(gca,'YTick',-.1:.2:.7,'FontSize',13,'LineWidth',2','Fontname','Ariel')
        set(gca,'XTick',0:100:500,'FontSize',13,'LineWidth',2','Fontname','Ariel')

        if plot_order(plotn-1)==5
            ylabel('Probe report probabilities','FontSize',12,'Fontname','Ariel')
            xlabel('Time from search task onset [ms]','FontSize',12,'Fontname','Ariel')
        end
        ylim([0 1]); xlim([0 540])
        axis square
    end
end
suptitle(sprintf('P1 & P2 - %i subj - probeGratPos : %s ',n_obs, probeGratPos));


%% P1 minus P2
Pdiff=P1_all-P2_all;

plotindiv=false;
h2=figure('Position',get(groot,'ScreenSize')); plotn=1;
for val=vals
    for cong=congrus
        subplot(size(congrus,2),size(vals,2),plot_order(plotn)); plotn=plotn+1; hold on;
        if plotindiv;
            for obs_i=1:n_obs
                plot(delays,Pdiff(:,val,cong,obs_i),'ko-','LineWidth',1,...
                    'MarkerFaceColor',[1 1 1],'MarkerSize',6,'Color',[0 0 0])
            end
        end
        errorbar(delays,mean(Pdiff(:,val,cong,:),4), std(Pdiff(:,val,cong,:),[],4)./sqrt(n_obs),...
            'ro-','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',[.2 .2 .2])
        
        set(gca,'YTick',-.1:.2:.7,'FontSize',13,'LineWidth',2','Fontname','Ariel')
        set(gca,'XTick',0:100:500,'FontSize',13,'LineWidth',2','Fontname','Ariel')
        
        if val==2; title('Valid'); elseif val==1; title('Invalid'); end
        if 5==plot_order(plotn-1)
            ylabel('difference P1-P2','FontSize',12,'Fontname','Ariel')
            xlabel('Time from search task onset [ms]','FontSize',12,'Fontname','Ariel')
        end
        ylim([-.3 .5]); xlim([0 540])
        axis square
    end
end
suptitle(sprintf('P1 minus P2 - %i subj - probeGratPos : %s ',n_obs, probeGratPos));


%% Amplitude spectrum

Fs = 1/.04; T = 1/Fs; L = 13-1; t = (0:L-1)*T; f = Fs*(0:(L/2))/L; freqs=f(2:end);
load([data_loc 'fft_p_All_invalid_11subjs.mat'])
load([data_loc 'fft_p_All_valid_11subjs.mat'])
fft_p_all=cat(3,fft_p_invalid, fft_p_valid);
ci=1-(.05/6);

% make custom black to green colormap
C = green_colormap(); C=C*1.9; e=2.7183;
C(:,1)=power(e,linspace(-3,-.6,size(C,1)));
rch=interp1(1:size(C,1),C(:,1),1:.01:size(C,1));
gch=interp1(1:size(C,1),C(:,2),1:.01:size(C,1));
bch=interp1(1:size(C,1),C(:,3),1:.01:size(C,1));
C=[rch;gch;bch]';

Pdiff=P1_all-P2_all;
fft_Pdiff=fft(Pdiff,13,1);
a_fft_Pdiff=abs(fft_Pdiff);
a_fft_Pdiff=a_fft_Pdiff(2:size(freqs,2)+1,:,:,:);

plotindiv=false; % should it plot individual observers' spectra?
plot_surrs_gradient=true; % should it plot the green gradient of surrogate values?
plot_surrs_lvl=false; % should it only plot the significance threshold you chose?
h2=figure('Position',get(groot,'ScreenSize')); plotn=1;
for val=vals
    for cong=congrus
        subplot(size(congrus,2),size(vals,2),plot_order(plotn)); hold on;
        if plotindiv
            for obs_i=1:n_obs
                plot(freqs,a_fft_Pdiff(:,val,cong,obs_i),'ko-','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',6,'Color',[.4 .4 .4])
            end
        end
        
        set(gca,'YTick',.5:.1:1.6,'FontSize',13,'LineWidth',2','Fontname','Ariel')
        set(gca,'XTick',f(2:end),'FontSize',13,'LineWidth',2','Fontname','Ariel')

        if plot_surrs_lvl || plot_surrs_gradient
            allftimefunction=sort(squeeze(fft_p_all(:,:,val)),1);
            percentile = floor(ci*repeatnumber);
            upperlim = allftimefunction(percentile,:);
            expected = mean(fft_p_all(:,:,val),1);
        end
        
        if plot_surrs_lvl
            plot(freqs,upperlim,'k--','LineWidth',1.5);
            plot(freqs,expected,'k-','LineWidth',1.5);
            plot(freqs,mean(a_fft_Pdiff(:,val,cong,:),4),'ko-','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',[.2 .2 .2])
            ylim([.55 1.42]);
            
        elseif plot_surrs_gradient
            surf(freqs,1:100000,allftimefunction);
            hold on; view(0,0); colormap(C); set(gcf,'Renderer','Zbuffer');
            shading interp;

            plot3(freqs,zeros(6),expected,'k-','LineWidth',1.5);
            plot3(freqs,zeros(6),upperlim,'k--','LineWidth',1.5);
            plot3(freqs,zeros(6),mean(a_fft_Pdiff(:,val,cong,:),4),'ko-','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',[.2 .2 .2])
            zlim([.55 1.42]);
            
        else
            plot(freqs,mean(a_fft_Pdiff(:,val,cong,:),4),'ko-','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',[.2 .2 .2])
            ylim([.55 1.42]);
        end

        if val==2; title('Valid'); elseif val==1; title('Invalid'); end
        if plot_order(plotn)==1
            ylabel('Amplitude of oscillation in P1-P2','FontSize',12,'Fontname','Ariel')
            xlabel('Frequency (Hz)','FontSize',12,'Fontname','Ariel')
        end
        xlim([freqs(1)-1 freqs(end)+1]); plotn=plotn+1; axis square
    end
end

suptitle(sprintf('fft diff P1 minus P2 - %i subj - probeGratPos : %s\np<0.05 (bonferroni corrected)',n_obs, probeGratPos));










