function summarize_AKBK
% Function: Combine results from all subjects to make summary bar plots for
% different classifiers

% Input: NONE
% Output: NONE

% Function dependencies:
% errorbar_groups.m

%%%%%
% Documented by: Blair Hu 08/04/17
%%%%%

subjnums = {'AB156','AB185','AB186','AB188','AB189','AB190','AB191','AB192','AB193','AB194'};
phasekeys = {'RHC','RTO','LHC','LTO'};

withinall = [];
withinall_ss = [];
withinall_t = [];

withinerror = {};
withinsserror = {};
withinterror = {};

% Specify which classification results to summarize
delay = 0;

% Specify whether to save the error rates combined across all subjects for
% the specified delay
savefile = 1;

close all

set(0,'DefaultAxesFontName', 'Palatino Linotype')

for i = 1:length(subjnums)
    disp(['Subject: ',subjnums{i}])
    % This should reflect the output of the Jupyter notebook
    load([subjnums{i},'_WithinSubjectResults_',num2str(delay),'Delay_ModeSpecific_AKBK_Classifiers_10Fold_Redo.mat']);
    eval(['RHC = ',subjnums{i},'.RHC;']);
    eval(['RTO = ',subjnums{i},'.RTO;']);
    eval(['LHC = ',subjnums{i},'.LHC;']);
    eval(['LTO = ',subjnums{i},'.LTO;']);
    
    withinsub{i,1} = RHC;
    withinsub{i,2} = RTO;
    withinsub{i,3} = LHC;
    withinsub{i,4} = LTO;
    
    withintemp = [withinsub{i,1}.Pred; withinsub{i,2}.Pred; withinsub{i,3}.Pred; withinsub{i,4}.Pred];
    trig = [withinsub{i,1}.Trig; withinsub{i,2}.Trig; withinsub{i,3}.Trig; withinsub{i,4}.Trig];
    
    classifiercount = 0;
    
    ss_or_t = nan(length(trig),1);
    
    % Identify steady-state and transitional steps and the gait phase
    for j = 1:length(ss_or_t)
        trigstr = trig{j};
        if ~isnan(trigstr)
            leavemodestr = trigstr(1);
            entermodestr = trigstr(3);
            if strcmp(leavemodestr,entermodestr)
                ss_or_t(j) = 1; % 1 means steady-state
            else
                ss_or_t(j) = 0; % 0 means transitional
            end
            leavemode(j,1) = str2double(trig{j}(1));
            entermode(j,1) = str2double(trig{j}(3));
            enterphase(j,1) = str2double(trig{j}(4));
        end
    end
    
    hcinds = find(enterphase == 1);
    toinds = find(enterphase == 2);
    ssinds = find(ss_or_t == 1);
    tinds = find(ss_or_t == 0);
    
    hc_ss = confusionmat(withintemp(intersect(hcinds,ssinds),9),withintemp(intersect(hcinds,ssinds),1));
    hc_t = confusionmat(withintemp(intersect(hcinds,tinds),9),withintemp(intersect(hcinds,tinds),1));
    to_ss = confusionmat(withintemp(intersect(toinds,ssinds),9),withintemp(intersect(toinds,ssinds),1));
    to_t = confusionmat(withintemp(intersect(toinds,tinds),9),withintemp(intersect(toinds,tinds),1));
    
    % Iterate across classifiers (there are 8 sensor sets)
    for classifier = 1:3:12
        classifiercount = classifiercount + 1;
        disp(['Classifier: ',num2str(classifiercount)]);
        
        % Allows us to add a row (for labeling sensor set) and a column
        % (for labeling classifier)
        headershift = 1;
        
        withinclassifiertemp = withintemp(:,[classifier:classifier+2 13]);
        
        % Gets error rates for each modality(columns)/classifier(rows)
        [withinerror{classifiercount+headershift,1+headershift}(i),withinsserror{classifiercount+headershift,1+headershift}(i),withinterror{classifiercount+headershift,1+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,4),withinclassifiertemp(:,1),ssinds,tinds);
        [withinerror{classifiercount+headershift,2+headershift}(i),withinsserror{classifiercount+headershift,2+headershift}(i),withinterror{classifiercount+headershift,2+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,4),withinclassifiertemp(:,2),ssinds,tinds);
        [withinerror{classifiercount+headershift,3+headershift}(i),withinsserror{classifiercount+headershift,3+headershift}(i),withinterror{classifiercount+headershift,3+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,4),withinclassifiertemp(:,3),ssinds,tinds);
    end
    % Save predictions across subjects
    withinall = [withinall; withintemp];
    withinall_ss = [withinall_ss; withintemp(ssinds,:)];
    withinall_t = [withinall_t; withintemp(tinds,:)];
end

% Add header labels
classifierheader = {'PCA LDA';'NORM SVM';'NORM ANN';'PCA SRC'};
modalityheader = {'AK','BK','IPSI_ALL'};

withinerror(2:5,1) = classifierheader;
withinsserror(2:5,1) = classifierheader;
withinterror(2:5,1) = classifierheader;

withinerror(1,2:4) = modalityheader;
withinsserror(1,2:4) = modalityheader;
withinterror(1,2:4) = modalityheader;

% Save summarized error rates, if specified
if savefile 
    save(['AllSubs_',num2str(delay),'Delay_AKBK_ErrorRates_Redo.mat'],'withinerror','withinsserror','withinterror');
end

% Calculate mean and SEM across subjects for each modality/classifier
for i = 2:5 % Iterate across classifiers
    for j = 2:4 % Iterate across sensor sets
        withinerror_mean{i,j} = mean(100*withinerror{i,j});
        withinerror_std{i,j} = std(100*withinerror{i,j});
        withinsserror_mean{i,j} = mean(100*withinsserror{i,j});
        withinsserror_std{i,j} = std(100*withinsserror{i,j});
        withinterror_mean{i,j} = mean(100*withinterror{i,j});
        withinterror_std{i,j} = std(100*withinterror{i,j});
    end
end

% Add header labels
withinerror_mean(2:5,1) = classifierheader;
withinerror_std(2:5,1) = classifierheader;
withinsserror_mean(2:5,1) = classifierheader;
withinsserror_std(2:5,1) = classifierheader;
withinterror_mean(2:5,1) = classifierheader;
withinterror_std(2:5,1) = classifierheader;

withinerror_mean(1,2:4) = modalityheader;
withinerror_std(1,2:4) = modalityheader;
withinsserror_mean(1,2:4) = modalityheader;
withinsserror_std(1,2:4) = modalityheader;
withinterror_mean(1,2:4) = modalityheader;
withinterror_std(1,2:4) = modalityheader;

chanuse = 2:4;

LDA_withinmean = cell2mat(withinerror_mean(2,chanuse));
LDA_withinstd = cell2mat(withinerror_std(2,chanuse));
LDA_withinssmean = cell2mat(withinsserror_mean(2,chanuse));
LDA_withinssstd = cell2mat(withinsserror_std(2,chanuse));
LDA_withintmean = cell2mat(withinterror_mean(2,chanuse));
LDA_withintstd = cell2mat(withinterror_std(2,chanuse));

SVM_withinmean = cell2mat(withinerror_mean(3,chanuse));
SVM_withinstd = cell2mat(withinerror_std(3,chanuse));
SVM_withinssmean = cell2mat(withinsserror_mean(3,chanuse));
SVM_withinssstd = cell2mat(withinsserror_std(3,chanuse));
SVM_withintmean = cell2mat(withinterror_mean(3,chanuse));
SVM_withintstd = cell2mat(withinterror_std(3,chanuse));

ANN_withinmean = cell2mat(withinerror_mean(4,chanuse));
ANN_withinstd = cell2mat(withinerror_std(4,chanuse));
ANN_withinssmean = cell2mat(withinsserror_mean(4,chanuse));
ANN_withinssstd = cell2mat(withinsserror_std(4,chanuse));
ANN_withintmean = cell2mat(withinterror_mean(4,chanuse));
ANN_withintstd = cell2mat(withinterror_std(4,chanuse));

SRC_withinmean = cell2mat(withinerror_mean(5,chanuse));
SRC_withinstd = cell2mat(withinerror_std(5,chanuse));
SRC_withinssmean = cell2mat(withinsserror_mean(5,chanuse));
SRC_withinssstd = cell2mat(withinsserror_std(5,chanuse));
SRC_withintmean = cell2mat(withinterror_mean(5,chanuse));
SRC_withintstd = cell2mat(withinterror_std(5,chanuse));

withinbar = [LDA_withinmean; SVM_withinmean; ANN_withinmean];
withinbarerr = [LDA_withinstd; SVM_withinstd; ANN_withinstd]/sqrt(10);

withinssbar = [LDA_withinssmean; SVM_withinssmean; ANN_withinssmean];
withinssbarerr = [LDA_withinssstd; SVM_withinssstd; ANN_withinssstd]/sqrt(10); 

withintbar = [LDA_withintmean; SVM_withintmean; ANN_withintmean];
withintbarerr = [LDA_withintstd; SVM_withintstd; ANN_withintstd]/sqrt(10);

% Create results matrices (for use with errorbar_groups) for single
% modalities plot and fused modalities plot
% withinbar_single = reshape(withinbar(1,1:9),3,3);
% withinbarerr_single = reshape(withinbarerr(1,1:9),3,3);
% withinbar_fused = withinbar(1:3,10:12)';
% withinbarerr_fused = withinbarerr(1:3,10:12)';

% Make bar plots
[~,~,singlebar_handles] = errorbar_groups(withinbar(1,:)',withinbarerr(1,:)','bar_colors',[128 128 255; 0 0 255; 0 0 128; 128 255 128]/255,'bar_names',{''},'bar_width',0.75,'errorbar_width',1.5,'optional_bar_arguments',{'LineWidth',0.1},...
    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',2})
ylabel('Overall error (%)');
set(gca,'FontSize',18,'FontWeight','bold');
% legend({'AK','BK','IPSI'},'Location','northoutside','Orientation','horizontal')

[~,~,singlebar_handles] = errorbar_groups(withinssbar(1,:)',withinssbarerr(1,:)','bar_colors',[128 128 255; 0 0 255; 0 0 128; 128 255 128]/255,'bar_names',{''},'bar_width',0.75,'errorbar_width',1.5,'optional_bar_arguments',{'LineWidth',0.1},...
    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',2})
ylabel('Steady-state error (%)');
set(gca,'FontSize',18,'FontWeight','bold');
% legend({'AK','BK','IPSI'},'Location','northoutside','Orientation','horizontal')

[~,~,singlebar_handles] = errorbar_groups(withintbar(1,:)',withintbarerr(1,:)','bar_colors',[128 128 255; 0 0 255; 0 0 128; 128 255 128]/255,'bar_names',{''},'bar_width',0.75,'errorbar_width',1.5,'optional_bar_arguments',{'LineWidth',0.1},...
    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',2})
ylabel('Transitional error (%)');
set(gca,'FontSize',18,'FontWeight','bold');
% legend({'AK','BK','IPSI'},'Location','northoutside','Orientation','horizontal')
end

% Function computes overall, steady-state, and transitional error rates
% from predictions
function [errorrate, sserror, terror] = calcerrorrate(true,pred,ssinds,tinds)
errorrate = sum(true~=pred)/length(true);
sserror = sum(true(ssinds) ~= pred(ssinds))/length(ssinds);
terror = sum(true(tinds) ~= pred(tinds))/length(tinds);
end