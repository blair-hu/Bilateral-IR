function combinesubs
% Function: Combine all data across subjects for classification in Python

% Input: NONE
% Output: "AllSubs_feats_reprocessed.mat"

% Function dependencies: NONE

%%%%%
% Documented by: Blair Hu 08/04/17
%%%%%

subjnums = {'AB156','AB185','AB186','AB188','AB189','AB190','AB191','AB192','AB193','AB194'};

legphase_combined = [];
subject_combined = {};
trig_combined = {};

for subs = 1:length(subjnums)
    load([subjnums{subs},'_feats_reprocessed_032818TD.mat']);

    if ~exist('feats_combined')
        feats_combined = cell(1,length(feats));
    end
    
	for i = 1:length(feats)
        feats_combined{i} = [feats_combined{i}; feats{i}];
    end
    legphase_combined = [legphase_combined; legphase];
    subject_combined = [subject_combined; subject];
    trig_combined = [trig_combined; trig];
end

save('AllSubs_feats_reprocessed_032818TD.mat','feats_combined','featlabels','legphase_combined','subject_combined','trig_combined');
end