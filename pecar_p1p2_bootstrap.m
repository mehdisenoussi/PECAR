function pecar_p1p2_bootstrap(data_loc, subj_pA, subj_pB, repeatnumber,validity, probeGratPos)
    
    n_delays=size(subj_pA,1); n_obs=size(subj_pA,2);
    fft_p = zeros(repeatnumber, n_delays);
    fft_ALL_p = zeros(repeatnumber, n_delays);
    
    for repeat = 1:repeatnumber
        if ~mod(repeat,500); disp(['Repeat number: ' num2str(repeat)]); end
        
        p1=zeros(n_delays,n_obs); p2=zeros(n_delays,n_obs);
        
        for sub = 1:size(subj_pA,2)
            rand_idx_p1 = randsample(1:n_delays,n_delays);
            rand_idx_p2 = randsample(1:n_delays,n_delays);

            p1(:,sub) = subj_pA(rand_idx_p1,sub);
            p2(:,sub) = subj_pB(rand_idx_p2,sub);
        end
        p1p2=p1-p2;
        fft_p(repeat,:) = mean( abs( fft( p1p2,[],2)),2)';
        fft_ALL_p(repeat,:) = abs( fft( mean( p1p2,2)))';
    end
    
    % save average fft amplitude computed by observer
    if validity==2; txtval=''; fft_p_valid=fft_p(:,2:((n_delays-1)/2+1));
    else txtval='in'; fft_p_invalid=fft_p(:,2:((n_delays-1)/2+1)); end
    save([data_loc sprintf('fft_p_%s_%svalid_11subjs.mat',probeGratPos,...
        txtval)], sprintf('fft_p_%svalid',txtval))
    
    % save average fft amplitude computed on average across observers
    if validity==2; txtval=''; fft_p_valid=fft_ALL_p(:,2:((n_delays-1)/2+1));
    else txtval='in'; fft_p_invalid=fft_ALL_p(:,2:((n_delays-1)/2+1)); end
    save([data_loc sprintf( 'fft_ALL_p_%s_%svalid_11subjs.mat',probeGratPos,...
        txtval)], sprintf('fft_p_%svalid',txtval)) 
end