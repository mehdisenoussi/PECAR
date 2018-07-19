% datareformat.m
% Simplify the PECAR data by separating valid and invalid data

clear all; close all

datadr = './processed_data/';
datafile = dir([datadr,'11obs_probe*.mat']);
load(fullfile(datadr, datafile.name))
d = probe_info_all;
validity = validity_all;
delays = delays_all;
nsubj = length(d);

%%
for ss = 1:nsubj
    data(ss).nt = size(d{ss},1);        % Number of trials

    curv = logical(validity{ss} - 1)';  % Pull current subj's validity data (recode to 0-invalid,1-valid)
    curd = d{ss};                       % Pull current subj's data from cell array
    curdelays = delays{ss}';
    
    % Save separate valid and invalid data but also concatenating the delay
    % data onto the last col.
    valid{ss} = [curd(curv,:),curdelays(curv)];
    invalid{ss} = [curd(~curv,:),curdelays(~curv)];
    
end

savefilename = fullfile(datadr,sprintf('datastruct_%iobs_valid_invalid', nsubj));
save(savefilename,'valid','invalid','delays')