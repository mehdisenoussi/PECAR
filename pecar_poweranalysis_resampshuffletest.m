% PECAR_D_resampshuffletest.m
%
% This analysis does not generate new data sets. Instead it uses the
% existing data and runs the analysis on resampled trials (alternative
% hypothesis) compared to shuffling the responses for trials (e.g.
% shuffling the right responses only).
%
% Predictions:
% If the relationship between responses is driving the result then this
% should be lost in the shuffle version and there should be no significant
% difference between P1 and P2.
%
% A power analysis can be carried out in the same way as before. Are the
% distributions of the two versions overlapping?
%
% 9/19/2017 JC Moreland
% 12/28/2017 JC Moreland updated to use D
% 10/07/2018 M Senoussi adapted code to new organization of the data

clear all; close all

seed = 9192017;
rng(seed)
Nsamp = 100000;  % Number of samples to bootstrap
p = .05; % What is our alpha?

datadr = './processed_data/';
datafile = dir([datadr,'datastruct*.mat']);
load(fullfile(datadr, datafile.name))

val_cond = {'valid', 'invalid'};

figure(1)

%% Resample and analyse, then shuffle and reanalyse for each run
for val_ind = 1:2
    % 1 = valid, 2 = invalid
    fprintf(sprintf('Run simulations for %s trials\n', val_cond{val_ind}));
    switch val_ind
        case 1 % valid
            data = valid;
            savename = 'valid';
        case 2 % invalid
            data = invalid;
            savename = 'invalid';
    end
    
    nsubj = length(data);
    
    D_resamp = zeros(nsubj,Nsamp);
    D_shuff = zeros(nsubj,Nsamp);
    
    P1_resamp = zeros(nsubj,Nsamp);
    P2_resamp = zeros(nsubj,Nsamp);
    
    P1_shuff = zeros(nsubj,Nsamp);
    P2_shuff = zeros(nsubj,Nsamp);
    
    for ss = 1:nsubj
        fprintf(sprintf('\tsubject %i\n', ss));
        nt = size(data{ss},1);  % Number of trials
        
        curd = data{ss}; % Pull current subj's data from cell array
        
        r1 = curd(:,5);     % Response 1 (1 = correct, 0 = incorrect)
        r2 = curd(:,10);     % Response 2
        
        for samp_n = 1:Nsamp
            %--- Resampling
            % Sample with replacement a sample of the same size for this
            % subject. Need to sample the index so the pairing of probe 1 and
            % two always remains matched.
            idx = datasample(1:nt,nt);
            r1tmp_resamp = r1(idx);
            r2tmp_resamp = r2(idx);
            
            % Do the normal analysis for this sample and store            
            Pboth = sum(r1tmp_resamp == 1 & r2tmp_resamp == 1)/nt;
            Pnone = sum(r1tmp_resamp == 0 & r2tmp_resamp == 0)/nt;
            [P1_resamp(ss,samp_n), P2_resamp(ss,samp_n), D_resamp(ss,samp_n)] =...
                quadratic_analysis(Pboth, Pnone);
            
            %--- Shuffle
            idx_sh = Shuffle(idx);
            r1tmp_sh = r1(idx_sh); % Only shuffle one set of responses
            r2tmp_sh = r2(idx);
            
            % Sanity check: If this shuffling works then it should keep
            % things paired and using the same index should return the mean
            % of the data under the original method.
            %         r1tmp_sh = r1(idx);
            
            % Do the normal analysis for this sample and store
            Pboth_sh = sum(r1tmp_sh == 1 & r2tmp_sh == 1)/nt;
            Pnone_sh = sum(r1tmp_sh == 0 & r2tmp_sh == 0)/nt;
            [P1_shuff(ss,samp_n), P2_shuff(ss,samp_n), D_shuff(ss,samp_n)] =...
                quadratic_analysis(Pboth_sh, Pnone_sh);
            
        end
    end
    
    %% print power achieved using each measure 
    % Discriminant power achieved
    D_resamp_mean = mean(D_resamp,1);
    D_shuff_mean = mean(D_shuff,1);
    % Where is the critical value between the two Discr. distributions?
    crit_discr = prctile(D_shuff_mean,(1-p)*100);
    power_achieved = sum(D_resamp_mean > crit_discr)/length(D_resamp_mean);
    fprintf('Power achieved by Discr. in %s cond.: %.1f%%\n',...
         val_cond{val_ind}, power_achieved*100)
    
    % P1-P2 power achieved
    P1_resamp_mean = mean(P1_resamp,1);
    P2_resamp_mean = mean(P2_resamp,1);
    P1_shuff_mean = mean(P1_shuff,1);
    P2_shuff_mean = mean(P2_shuff,1);
    P1P2_resamp_diff = P1_resamp_mean - P2_resamp_mean;
    P1P2_shuff_diff = P1_shuff_mean - P2_shuff_mean;
    % Where is the critical value between the two P1-P2 distributions?
    crit_pdiff = prctile(P1P2_shuff_diff,(1-p)*100);
    power_achieved = sum(P1P2_resamp_diff > crit_pdiff)/length(P1P2_resamp_diff);
    fprintf('Power achieved by Pdiff in %s cond.: %.1f%%\n',...
        val_cond{val_ind}, power_achieved*100)
    
    %% Plot distributions of simulated data under the two scenarios above
    
    % Plot results for Discr.
    subplot(2,2,val_ind+2)
    hold on
    histogram(D_shuff_mean,'Binwidth',.002,'Normalization','probability',...
        'FaceColor',.25*[1 1 1],'FaceAlpha',.75,'EdgeAlpha',0)
    histogram(D_resamp_mean,'Binwidth',.002,'Normalization','probability',...
        'FaceColor',.5*[1 1 1],'FaceAlpha',.75,'EdgeAlpha',0)
    
    xlim([-.08,.08]); ylim([0,.12])
    xlabel('Discr'); ylabel('Proportion')
    set(gca,'Linewidth',2,'FontSize',13)
    
    plot([crit_discr, crit_discr],get(gca,'Ylim'),'k-','LineWidth',2)
    plot([0,0],get(gca,'Ylim'),'k--','LineWidth',1.5)
    
    % Plot results for PDiff
    subplot(2,2,val_ind); hold on
    title(val_cond{val_ind})
    
    histogram(P1P2_shuff_diff,...
        'Binwidth',.01,'Normalization','probability',...
        'FaceColor',.25*[1 1 1],'FaceAlpha',.75,'EdgeAlpha',0)
    
    histogram(P1P2_resamp_diff,...
        'Binwidth',.01,'Normalization','probability',...
        'FaceColor',.5*[1 1 1],'FaceAlpha',.75,'EdgeAlpha',0)
    xlim([-.3,.3]); ylim([0,.25])
    xlabel('P1 - P2'); ylabel('Proportion')
    
    plot([crit_pdiff, crit_pdiff],get(gca,'Ylim'),'k-','LineWidth',2)
    plot([0,0],get(gca,'Ylim'),'k--','LineWidth',1.5)
    if val_ind == 1
        legend({'Shuffled','Resampled','p-crit'}, 'Location','Best')
        legend boxoff
    end

    set(gca, 'LineWidth',1.5, 'FontSize',13,'FontName','Ariel')
    
end