<<<<<<< HEAD:Matlab/sfs_mslda.m
function [start_err,start_ss_err,start_t_err,chan_include,chan_selected,error_include,error_ss_include,error_t_include] = sfs_mslda(filename,chanPrev,checkChan,totChan)
% Function: Perform feature extraction for all modalities

% Input: filename (ABXXX_toreclassify.mat or ABXXX_feats_reprocessed.mat), channels already included, candidate channels, total channels to select  
% Output: starting error rate, updated set of channels already included,
% labels of selected channels, updated error rates

% Function dependencies:
% NONE

%%%%%
% Documented by: Blair Hu 08/04/17
%%%%%

warning('off','all')
load(filename);

% This corresponds to features from no delay window
features = feats{1};

% Specify cross-validation
kfold = 10;

usedchan = [];
unusedchan = [];

featsdict = {22,2};
% Create dictionary to map channels with their associated features
% Channel 1: Knee (I)
featsdict{1,1} = 'Ipsi Knee';
% Channel 2: Ankle (I)
featsdict{2,1} = 'Ipsi Ankle';
% Channel 3: Knee (C)
featsdict{3,1} = 'Contra Knee';
% Channel 4: Ankle (C)
featsdict{4,1} = 'Contra Ankle';
% Channel 5: TA (I)
featsdict{5,1} = 'Ipsi TA';
% Channel 6: MG (I)
featsdict{6,1} = 'Ipsi MG';
% Channel 7: SOL (I)
featsdict{7,1} = 'Ipsi SOL';
% Channel 8: BF (I)
featsdict{8,1} = 'Ipsi BF';
% Channel 9: ST (I)
featsdict{9,1} = 'Ipsi ST';
% Channel 10: VL (I)
featsdict{10,1} = 'Ipsi VL';
% Channel 11: RF (I)
featsdict{11,1} = 'Ipsi RF';
% Channel 12: TA (C)
featsdict{12,1} = 'Contra TA';
% Channel 13: MG (C)
featsdict{13,1} = 'Contra MG';
% Channel 14: SOL (C)
featsdict{14,1} = 'Contra SOL';
% Channel 15: BF (C)
featsdict{15,1} = 'Contra BF';
% Channel 16: ST (C)
featsdict{16,1} = 'Contra ST';
% Channel 17: VL (C)
featsdict{17,1} = 'Contra VL';
% Channel 18: RF (C)
featsdict{18,1} = 'Contra RF';
% Channel 19: Shank (I)
featsdict{19,1} = 'Ipsi Shank';
% Channel 20: Thigh (I)
featsdict{20,1} = 'Ipsi Thigh';
% Channel 21: Shank (C)
featsdict{21,1} = 'Contra Shank';
% Channel 22: Thigh (C)
featsdict{22,1} = 'Contra Thigh';

% Find indices for features associated with each channel
for chan = 1:22
    tempindex = strfind(featlabels,featsdict{chan,1});
    tempindex = find(not(cellfun('isempty',tempindex)));
    featsdict{chan,2} = tempindex;
end

% Find the leaving and entering mode for each trigger
for trigind = 1:length(trig)
    leavemode_orig(trigind,1) = str2num(trig{trigind}(1));
    entermode_orig(trigind,1) = str2num(trig{trigind}(3));
end

% Remove steps that start or end with sitting/standing
features = features(find((leavemode_orig > 0) & (leavemode_orig < 6) & (entermode_orig > 0) & (entermode_orig < 6)),:);
legphase = legphase(find((leavemode_orig > 0) & (leavemode_orig < 6) & (entermode_orig > 0) & (entermode_orig < 6)));
leavemode = leavemode_orig(find((leavemode_orig > 0) & (leavemode_orig < 6) & (entermode_orig > 0) & (entermode_orig < 6)));
entermode = entermode_orig(find((leavemode_orig > 0) & (leavemode_orig < 6) & (entermode_orig > 0) & (entermode_orig < 6)));

% Simple error checking
if totChan > length(checkChan)
    error('Not enough input channels.')
end
if ismember(chanPrev,checkChan)
    error('Checking an already selected channel.')
end

% Initialize variables
chan_include = [];
if isempty(chan_include) && isempty(chanPrev)
    chan_prev = [];
    start_err = [];
else
    chan_prev = [];
    for i = 1:length(chanPrev)
        chan_prev = [chan_prev; featsdict{chanPrev(i),2}];
    end
    [start_err,start_ss_err,start_t_err] = evalchan(chan_prev,features,legphase,leavemode,entermode,kfold);
end
chan_exclude = checkChan;
error_include = [];
error_ss_include = [];
error_t_include = [];

% Iterate until the total number of desired channels is reached
for i = 1:totChan
    temperrors = [];
    % Try candidate channels from chan_exclude
    for j = 1:length(chan_exclude)
        disp(['Trying channel ',num2str(chan_exclude(j)),'...'])
        [temperrors(j),sserrors(j),terrors(j)] = evalchan([chan_prev; featsdict{chan_exclude(j),2}],features,legphase,leavemode,entermode,kfold);
    end
    % Find the channel with smallest error rate and add to chan_include
    [min_error, min_chan] = min(temperrors);
    min_error_ss = sserrors(min_chan);
    min_error_t = terrors(min_chan);
    disp(['Added channel ',num2str(chan_exclude(min_chan)),': ',featsdict{chan_exclude(min_chan),1}])
    chan_include = [chan_include; chan_exclude(min_chan)];
    chan_prev = [chan_prev; featsdict{chan_exclude(min_chan),2}];
    error_include = [error_include; min_error];
    error_ss_include = [error_ss_include; min_error_ss];
    error_t_include = [error_t_include; min_error_t];
    % Remove the added channel
    chan_exclude = setdiff(checkChan,chan_include);
end

% Find labels for added channels
chan_selected = {};
for m = 1:length(chan_include)
    chan_selected = [chan_selected; featsdict{chan_include(m),1}];
end    
end

% evalchan replicates mode-specific classification
function [chan_error,chan_error_ss,chan_error_t] = evalchan(usecol,features,legphase,leavemode,entermode,kfold)
% Initialize variable to store all predictions
allphase_pred = [];
% Identify steady-state and transitional steps
ss_or_t = ones(length(entermode),1);
for i = 1:length(entermode)
    if entermode(i) ~= leavemode(i)
        ss_or_t(i) = 0; % Define transitional steps as 0
    end
end
for j = 1:4 % Iterate across legphases
    feats_phase = features(find(legphase == j),:);
    leavemode_phase = leavemode(find(legphase == j));
    entermode_phase = entermode(find(legphase == j));
    ss_or_t_phase = ss_or_t(find(legphase == j));
    
    allmode_pred = [];
    for k = 1:5 % Iterate across possible classes
        featstemp = feats_phase(find(leavemode_phase == k),:);
        entermodetemp = entermode_phase(find(leavemode_phase == k));
        ss_or_t_temp = ss_or_t_phase(find(leavemode_phase == k));
        
        cp = cvpartition(entermodetemp,'KFold',kfold);
        
        mode_pred = [];
        for fold = 1:cp.NumTestSets % Iterate across folds for cross-validation
            train_inds = find(training(cp,fold));
            test_inds = find(test(cp,fold));
            
            feats_train = featstemp(train_inds,usecol);
            mode_train = entermodetemp(train_inds);
            feats_test = featstemp(test_inds,usecol);
            mode_test = entermodetemp(test_inds);
            
            % Perform feature normalization
            [feats_train_norm,mean,sd] = zscore(feats_train);            
            feats_test_norm = (feats_test - repmat(mean,length(test_inds),1))./repmat(sd,length(test_inds),1);
            
            % Perform dimensionality reduction with PCA
            [weight,feats_train_pca,~,~,vaf] = pca(feats_train_norm);
            feats_test_pca = feats_test_norm*weight;
            
            pca_numcomps = min(find(cumsum(vaf) > 95));
            
            % Uses generic Matlab LDA classifier
            test_pred = classify(feats_test_pca(:,1:pca_numcomps),feats_train_pca(:,1:pca_numcomps),mode_train,'linear');
            
            mode_pred(test_inds,1) = test_pred;
            mode_pred(test_inds,2) = mode_test;
            mode_pred(test_inds,3) = ss_or_t_temp(test_inds);
        end
        allmode_pred = [allmode_pred; mode_pred];
    end
    allphase_pred = [allphase_pred; allmode_pred];
end
% Calculate the overall error across both legs for all phases
chan_error = sum(allphase_pred(:,1) ~= allphase_pred(:,2))/length(allphase_pred(:,1));
% Calculate the steady state error
chan_error_ss = sum(allphase_pred(find(allphase_pred(:,3)==1),1) ~= allphase_pred(find(allphase_pred(:,3)==1),2))/length(allphase_pred(find(allphase_pred(:,3)==1),1));
% Calculate the transitional error
chan_error_t = sum(allphase_pred(find(allphase_pred(:,3)==0),1) ~= allphase_pred(find(allphase_pred(:,3)==0),2))/length(allphase_pred(find(allphase_pred(:,3)==0),1));
=======
function [start_err,start_ss_err,start_t_err,chan_include,chan_selected,error_include,error_ss_include,error_t_include] = sfs_mslda(filename,chanPrev,checkChan,totChan)
% Function: Perform feature extraction for all modalities

% Input: filename (ABXXX_toreclassify.mat or ABXXX_feats_reprocessed.mat), channels already included, candidate channels, total channels to select  
% Output: starting error rate, updated set of channels already included,
% labels of selected channels, updated error rates

% Function dependencies:
% NONE

%%%%%
% Documented by: Blair Hu 08/04/17
%%%%%

warning('off','all')
load(filename);

% This corresponds to features from no delay window
features = feats{1};

% Specify cross-validation
kfold = 10;

usedchan = [];
unusedchan = [];

featsdict = {22,2};
% Create dictionary to map channels with their associated features
% Channel 1: Knee (I)
featsdict{1,1} = 'Ipsi Knee';
% Channel 2: Ankle (I)
featsdict{2,1} = 'Ipsi Ankle';
% Channel 3: Knee (C)
featsdict{3,1} = 'Contra Knee';
% Channel 4: Ankle (C)
featsdict{4,1} = 'Contra Ankle';
% Channel 5: TA (I)
featsdict{5,1} = 'Ipsi TA';
% Channel 6: MG (I)
featsdict{6,1} = 'Ipsi MG';
% Channel 7: SOL (I)
featsdict{7,1} = 'Ipsi SOL';
% Channel 8: BF (I)
featsdict{8,1} = 'Ipsi BF';
% Channel 9: ST (I)
featsdict{9,1} = 'Ipsi ST';
% Channel 10: VL (I)
featsdict{10,1} = 'Ipsi VL';
% Channel 11: RF (I)
featsdict{11,1} = 'Ipsi RF';
% Channel 12: TA (C)
featsdict{12,1} = 'Contra TA';
% Channel 13: MG (C)
featsdict{13,1} = 'Contra MG';
% Channel 14: SOL (C)
featsdict{14,1} = 'Contra SOL';
% Channel 15: BF (C)
featsdict{15,1} = 'Contra BF';
% Channel 16: ST (C)
featsdict{16,1} = 'Contra ST';
% Channel 17: VL (C)
featsdict{17,1} = 'Contra VL';
% Channel 18: RF (C)
featsdict{18,1} = 'Contra RF';
% Channel 19: Shank (I)
featsdict{19,1} = 'Ipsi Shank';
% Channel 20: Thigh (I)
featsdict{20,1} = 'Ipsi Thigh';
% Channel 21: Shank (C)
featsdict{21,1} = 'Contra Shank';
% Channel 22: Thigh (C)
featsdict{22,1} = 'Contra Thigh';

% Find indices for features associated with each channel
for chan = 1:22
    tempindex = strfind(featlabels,featsdict{chan,1});
    tempindex = find(not(cellfun('isempty',tempindex)));
    featsdict{chan,2} = tempindex;
end

% Find the leaving and entering mode for each trigger
for trigind = 1:length(trig)
    leavemode_orig(trigind,1) = str2num(trig{trigind}(1));
    entermode_orig(trigind,1) = str2num(trig{trigind}(3));
end

% Remove steps that start or end with sitting/standing
features = features(find((leavemode_orig > 0) & (leavemode_orig < 6) & (entermode_orig > 0) & (entermode_orig < 6)),:);
legphase = legphase(find((leavemode_orig > 0) & (leavemode_orig < 6) & (entermode_orig > 0) & (entermode_orig < 6)));
leavemode = leavemode_orig(find((leavemode_orig > 0) & (leavemode_orig < 6) & (entermode_orig > 0) & (entermode_orig < 6)));
entermode = entermode_orig(find((leavemode_orig > 0) & (leavemode_orig < 6) & (entermode_orig > 0) & (entermode_orig < 6)));

% Simple error checking
if totChan > length(checkChan)
    error('Not enough input channels.')
end
if ismember(chanPrev,checkChan)
    error('Checking an already selected channel.')
end

% Initialize variables
chan_include = [];
if isempty(chan_include) && isempty(chanPrev)
    chan_prev = [];
    start_err = [];
else
    chan_prev = [];
    for i = 1:length(chanPrev)
        chan_prev = [chan_prev; featsdict{chanPrev(i),2}];
    end
    [start_err,start_ss_err,start_t_err] = evalchan(chan_prev,features,legphase,leavemode,entermode,kfold);
end
chan_exclude = checkChan;
error_include = [];
error_ss_include = [];
error_t_include = [];

% Iterate until the total number of desired channels is reached
for i = 1:totChan
    temperrors = [];
    % Try candidate channels from chan_exclude
    for j = 1:length(chan_exclude)
        disp(['Trying channel ',num2str(chan_exclude(j)),'...'])
        [temperrors(j),sserrors(j),terrors(j)] = evalchan([chan_prev; featsdict{chan_exclude(j),2}],features,legphase,leavemode,entermode,kfold);
    end
    % Find the channel with smallest error rate and add to chan_include
    [min_error, min_chan] = min(temperrors);
    min_error_ss = sserrors(min_chan);
    min_error_t = terrors(min_chan);
    disp(['Added channel ',num2str(chan_exclude(min_chan)),': ',featsdict{chan_exclude(min_chan),1}])
    chan_include = [chan_include; chan_exclude(min_chan)];
    chan_prev = [chan_prev; featsdict{chan_exclude(min_chan),2}];
    error_include = [error_include; min_error];
    error_ss_include = [error_ss_include; min_error_ss];
    error_t_include = [error_t_include; min_error_t];
    % Remove the added channel
    chan_exclude = setdiff(checkChan,chan_include);
end

% Find labels for added channels
chan_selected = {};
for m = 1:length(chan_include)
    chan_selected = [chan_selected; featsdict{chan_include(m),1}];
end    
end

% evalchan replicates mode-specific classification
function [chan_error,chan_error_ss,chan_error_t] = evalchan(usecol,features,legphase,leavemode,entermode,kfold)
% Initialize variable to store all predictions
allphase_pred = [];
% Identify steady-state and transitional steps
ss_or_t = ones(length(entermode),1);
for i = 1:length(entermode)
    if entermode(i) ~= leavemode(i)
        ss_or_t(i) = 0; % Define transitional steps as 0
    end
end
for j = 1:4 % Iterate across legphases
    feats_phase = features(find(legphase == j),:);
    leavemode_phase = leavemode(find(legphase == j));
    entermode_phase = entermode(find(legphase == j));
    ss_or_t_phase = ss_or_t(find(legphase == j));
    
    allmode_pred = [];
    for k = 1:5 % Iterate across possible classes
        featstemp = feats_phase(find(leavemode_phase == k),:);
        entermodetemp = entermode_phase(find(leavemode_phase == k));
        ss_or_t_temp = ss_or_t_phase(find(leavemode_phase == k));
        
        cp = cvpartition(entermodetemp,'KFold',kfold);
        
        mode_pred = [];
        for fold = 1:cp.NumTestSets % Iterate across folds for cross-validation
            train_inds = find(training(cp,fold));
            test_inds = find(test(cp,fold));
            
            feats_train = featstemp(train_inds,usecol);
            mode_train = entermodetemp(train_inds);
            feats_test = featstemp(test_inds,usecol);
            mode_test = entermodetemp(test_inds);
            
            % Perform feature normalization
            [feats_train_norm,mean,sd] = zscore(feats_train);            
            feats_test_norm = (feats_test - repmat(mean,length(test_inds),1))./repmat(sd,length(test_inds),1);
            
            % Perform dimensionality reduction with PCA
            [weight,feats_train_pca,~,~,vaf] = pca(feats_train_norm);
            feats_test_pca = feats_test_norm*weight;
            
            pca_numcomps = min(find(cumsum(vaf) > 95));
            
            % Uses generic Matlab LDA classifier
            test_pred = classify(feats_test_pca(:,1:pca_numcomps),feats_train_pca(:,1:pca_numcomps),mode_train,'linear');
            
            mode_pred(test_inds,1) = test_pred;
            mode_pred(test_inds,2) = mode_test;
            mode_pred(test_inds,3) = ss_or_t_temp(test_inds);
        end
        allmode_pred = [allmode_pred; mode_pred];
    end
    allphase_pred = [allphase_pred; allmode_pred];
end
% Calculate the overall error across both legs for all phases
chan_error = sum(allphase_pred(:,1) ~= allphase_pred(:,2))/length(allphase_pred(:,1));
% Calculate the steady state error
chan_error_ss = sum(allphase_pred(find(allphase_pred(:,3)==1),1) ~= allphase_pred(find(allphase_pred(:,3)==1),2))/length(allphase_pred(find(allphase_pred(:,3)==1),1));
% Calculate the transitional error
chan_error_t = sum(allphase_pred(find(allphase_pred(:,3)==0),1) ~= allphase_pred(find(allphase_pred(:,3)==0),2))/length(allphase_pred(find(allphase_pred(:,3)==0),1));
>>>>>>> 9ee07b7770d5b0c88658f9c3fb8701cbd201db21:Matlab/sfs_mslda.m
end