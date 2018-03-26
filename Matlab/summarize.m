function summarize
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
savefile = 0;

close all

set(0,'DefaultAxesFontName', 'Palatino Linotype')

% 1-29 LDA
% 30-58 SVM
% 59-87 ANN
% 88-116 SRC
% 117 Truth

for i = 1:length(subjnums)
    disp(['Subject: ',subjnums{i}])
    % This should reflect the output of the Jupyter notebook
    load([subjnums{i},'_WithinSubjectResults_',num2str(delay),'Delay_ModeSpecific_AllSensorsClassifiers_10Fold_Redo.mat']);
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
    
%     hc_ss = confusionmat(withintemp(intersect(hcinds,ssinds),117),withintemp(intersect(hcinds,ssinds),7));
%     hc_t = confusionmat(withintemp(intersect(hcinds,tinds),117),withintemp(intersect(hcinds,tinds),7));
%     to_ss = confusionmat(withintemp(intersect(toinds,ssinds),117),withintemp(intersect(toinds,ssinds),7));
%     to_t = confusionmat(withintemp(intersect(toinds,tinds),117),withintemp(intersect(toinds,tinds),7));
    
    % Iterate across classifiers (there are 29 sensor sets)
    for classifier = 1:33:132
        classifiercount = classifiercount + 1;
        disp(['Classifier: ',num2str(classifiercount)]);
        
        % Allows us to add a row (for labeling sensor set) and a column
        % (for labeling classifier)
        headershift = 1;
        
        withinclassifiertemp = withintemp(:,[classifier:classifier+28 117]);
        
        % Gets error rates for each modality(columns)/classifier(rows)
        for sensorcombo = 1:33
            [within_overall, within_ss, within_t] = calcerrorrate(withinclassifiertemp(:,34),withinclassifiertemp(:,sensorcombo),ssinds,tinds);
            withinerror{classifiercount+headershift,sensorcombo+headershift}(i) = within_overall;
            withinsserror{classifiercount+headershift,sensorcombo+headershift}(i) = within_overall;
            withinterror{classifiercount+headershift,sensorcombo+headershift}(i) = within_overall;
        end
%         [withinerror{classifiercount+headershift,1+headershift}(i),withinsserror{classifiercount+headershift,1+headershift}(i),withinterror{classifiercount+headershift,1+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,1),ssinds,tinds);
%         [withinerror{classifiercount+headershift,2+headershift}(i),withinsserror{classifiercount+headershift,2+headershift}(i),withinterror{classifiercount+headershift,2+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,2),ssinds,tinds);
%         [withinerror{classifiercount+headershift,3+headershift}(i),withinsserror{classifiercount+headershift,3+headershift}(i),withinterror{classifiercount+headershift,3+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,3),ssinds,tinds);
%         [withinerror{classifiercount+headershift,4+headershift}(i),withinsserror{classifiercount+headershift,4+headershift}(i),withinterror{classifiercount+headershift,4+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,4),ssinds,tinds);
%         [withinerror{classifiercount+headershift,5+headershift}(i),withinsserror{classifiercount+headershift,5+headershift}(i),withinterror{classifiercount+headershift,5+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,5),ssinds,tinds);
%         [withinerror{classifiercount+headershift,6+headershift}(i),withinsserror{classifiercount+headershift,6+headershift}(i),withinterror{classifiercount+headershift,6+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,6),ssinds,tinds);
%         [withinerror{classifiercount+headershift,7+headershift}(i),withinsserror{classifiercount+headershift,7+headershift}(i),withinterror{classifiercount+headershift,7+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,7),ssinds,tinds);
%         [withinerror{classifiercount+headershift,8+headershift}(i),withinsserror{classifiercount+headershift,8+headershift}(i),withinterror{classifiercount+headershift,8+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,8),ssinds,tinds);
%         [withinerror{classifiercount+headershift,9+headershift}(i),withinsserror{classifiercount+headershift,9+headershift}(i),withinterror{classifiercount+headershift,9+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,9),ssinds,tinds);
%         [withinerror{classifiercount+headershift,10+headershift}(i),withinsserror{classifiercount+headershift,10+headershift}(i),withinterror{classifiercount+headershift,10+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,10),ssinds,tinds);
%         [withinerror{classifiercount+headershift,11+headershift}(i),withinsserror{classifiercount+headershift,11+headershift}(i),withinterror{classifiercount+headershift,11+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,11),ssinds,tinds);
%         [withinerror{classifiercount+headershift,12+headershift}(i),withinsserror{classifiercount+headershift,12+headershift}(i),withinterror{classifiercount+headershift,12+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,12),ssinds,tinds);
%         [withinerror{classifiercount+headershift,13+headershift}(i),withinsserror{classifiercount+headershift,13+headershift}(i),withinterror{classifiercount+headershift,13+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,13),ssinds,tinds);
%         [withinerror{classifiercount+headershift,14+headershift}(i),withinsserror{classifiercount+headershift,14+headershift}(i),withinterror{classifiercount+headershift,14+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,14),ssinds,tinds);
%         [withinerror{classifiercount+headershift,15+headershift}(i),withinsserror{classifiercount+headershift,15+headershift}(i),withinterror{classifiercount+headershift,15+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,15),ssinds,tinds);
%         [withinerror{classifiercount+headershift,16+headershift}(i),withinsserror{classifiercount+headershift,16+headershift}(i),withinterror{classifiercount+headershift,16+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,16),ssinds,tinds);
%         [withinerror{classifiercount+headershift,17+headershift}(i),withinsserror{classifiercount+headershift,17+headershift}(i),withinterror{classifiercount+headershift,17+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,17),ssinds,tinds);
%         [withinerror{classifiercount+headershift,18+headershift}(i),withinsserror{classifiercount+headershift,18+headershift}(i),withinterror{classifiercount+headershift,18+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,18),ssinds,tinds);
%         [withinerror{classifiercount+headershift,19+headershift}(i),withinsserror{classifiercount+headershift,19+headershift}(i),withinterror{classifiercount+headershift,19+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,19),ssinds,tinds);
%         [withinerror{classifiercount+headershift,20+headershift}(i),withinsserror{classifiercount+headershift,20+headershift}(i),withinterror{classifiercount+headershift,20+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,20),ssinds,tinds);
%         [withinerror{classifiercount+headershift,21+headershift}(i),withinsserror{classifiercount+headershift,21+headershift}(i),withinterror{classifiercount+headershift,21+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,21),ssinds,tinds);
%         [withinerror{classifiercount+headershift,22+headershift}(i),withinsserror{classifiercount+headershift,22+headershift}(i),withinterror{classifiercount+headershift,22+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,22),ssinds,tinds);
%         [withinerror{classifiercount+headershift,23+headershift}(i),withinsserror{classifiercount+headershift,23+headershift}(i),withinterror{classifiercount+headershift,23+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,23),ssinds,tinds);
%         [withinerror{classifiercount+headershift,24+headershift}(i),withinsserror{classifiercount+headershift,24+headershift}(i),withinterror{classifiercount+headershift,24+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,24),ssinds,tinds);
%         [withinerror{classifiercount+headershift,25+headershift}(i),withinsserror{classifiercount+headershift,25+headershift}(i),withinterror{classifiercount+headershift,25+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,25),ssinds,tinds);
%         [withinerror{classifiercount+headershift,26+headershift}(i),withinsserror{classifiercount+headershift,26+headershift}(i),withinterror{classifiercount+headershift,26+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,26),ssinds,tinds);
%         [withinerror{classifiercount+headershift,27+headershift}(i),withinsserror{classifiercount+headershift,27+headershift}(i),withinterror{classifiercount+headershift,27+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,27),ssinds,tinds);
%         [withinerror{classifiercount+headershift,28+headershift}(i),withinsserror{classifiercount+headershift,28+headershift}(i),withinterror{classifiercount+headershift,28+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,28),ssinds,tinds);
%         [withinerror{classifiercount+headershift,29+headershift}(i),withinsserror{classifiercount+headershift,29+headershift}(i),withinterror{classifiercount+headershift,29+headershift}(i)] = calcerrorrate(withinclassifiertemp(:,30),withinclassifiertemp(:,29),ssinds,tinds);
    end
    % Save predictions across subjects
    withinall = [withinall; withintemp];
    withinall_ss = [withinall_ss; withintemp(ssinds,:)];
    withinall_t = [withinall_t; withintemp(tinds,:)];
end

% Add header labels
classifierheader = {'PCA LDA';'NORM SVM';'NORM ANN';'PCA SRC'};
modalityheader = {'EMG(I)','GONIO(I)','IMU(I)','EMG+GONIO(I)','EMG+IMU(I)','GONIO+IMU(I)','IPSI ALL','I/EMG(C)','I/GONIO(C)','I/IMU(C)','I/EMG+GONIO(C)','I/EMG+IMU(C)','I/GONIO+IMU(C)','BILAT EMG','BILAT GONIO','BILAT IMU','BILAT ALL',...
    'EMG LL(I)','EMG UL(I)','EMG(I)/EMG LL(C)','EMG(I)/EMG UL(C)','KNEE(I)','ANKLE(I)','GONIO(I)/KNEE(C)','GONIO(I)/ANKLE(C)','SHANK(I)','THIGH(I)','IMU(I)/SHANK(C)','IMU(I)/THIGH(C)','EMG(C)','GONIO(C)','IMU(C)','ALL(C)'};

withinerror(2:5,1) = classifierheader;
withinsserror(2:5,1) = classifierheader;
withinterror(2:5,1) = classifierheader;

withinerror(1,2:34) = modalityheader;
withinsserror(1,2:34) = modalityheader;
withinterror(1,2:34) = modalityheader;

% Save summarized error rates, if specified
if savefile 
    save(['AllSubs_',num2str(delay),'Delay_ErrorRates_032618.mat'],'withinerror','withinsserror','withinterror');
end

% Calculate mean and SEM across subjects for each modality/classifier
for i = 2:5 % Iterate across classifiers
    for j = 2:34 % Iterate across sensor sets
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

withinerror_mean(1,2:34) = modalityheader;
withinerror_std(1,2:34) = modalityheader;
withinsserror_mean(1,2:34) = modalityheader;
withinsserror_std(1,2:34) = modalityheader;
withinterror_mean(1,2:34) = modalityheader;
withinterror_std(1,2:34) = modalityheader;

% Select which sensor set to plot
info = 5;

% Define plot parameters based on selection above (commented section used
% if sparse representation classifier is included)
if info == 1 % EMG only
    titlestr = 'EMG';
    chanuse = [19 20 2 15];
    xticks = [2:3:12];
    xticklabels = {'LL(I)','UL(I)','LL+UL(I)','LL+UL(B)'};
%     chanuse = [19 20 2 21 22 15];
%     xticks = [2:4:24];
%     xticklabels = {'LL(I)','UL(I)','ALL(I)','ALL(I)+LL(C)','ALL(I)+UL(C)','ALL(B)'};
elseif info == 2 % GONIO only
    titlestr = 'GONIO';
    chanuse = [23 24 3 16];
    xticks = [2:3:12];
    xticklabels = {'K(I)','A(I)','K+A(I)','K+A(B)'};
%     chanuse = [23 24 3 25 26 16];
%     xticks = [2:4:24]; 
%     xticklabels = {'K(I)','A(I)','K-A(I)','K-A(I)+K(C)','K-A(I)+A(C)','K-A(B)'};
elseif info == 3 % IMU only
    titlestr = 'IMU';
    chanuse = [27 28 4 17];
    xticks = [2:3:12];
    xticklabels = {'S(I)','T(I)','S+T(I)','S+T(B)'};
%     chanuse = [27 28 4 29 30 17];
%     xticks = [2:4:24];
%     xticklabels = {'S(I)','T(I)','S-T(I)','S-T(I)+S(C)','S-T(I)+T(C)','S-T(B)'};
elseif info == 4 % FUSED
    titlestr = 'FUSED';
    chanuse = [5 6 7 8 18];
    xticks = [2:3:15];
    xticklabels = {['     EMG','\newline','+GONIO(I)'],['  EMG','\newline','+IMU(I)'],['GONIO','\newline','+IMU(I)'],'ALL(I)','ALL(B)'};
%     chanuse = [5 6 7 8 9 10 11 12 13 14 18];
%     xticks = [2:4:44];
%     xticklabels = {'EMG-GONIO(I)','EMG-IMU(I)','GONIO-IMU(I)','ALL(I)','ALL(I)+EMG(C)','ALL(I)+GONIO(C)','ALL(I)+IMU(C)','ALL(I)+EMG-GONIO(C)','ALL(I)+EMG-IMU(C)','ALL(I)+GONIO-IMU(C)','ALL(B)'};    
else % Paper comparisons
    titlestr = '';
    chanuse = [2 31 15 3 32 16 4 33 17 8 34 18];
    xticks = [2:3:24];
    xticklabels = {'EMG(I)','EMG(C)','EMG(B)','GONIO(I)','GONIO(C)','GONIO(B)','IMU(I)','IMU(C)','IMU(B)','ALL(I)','ALL(C)','ALL(B)'};
%     xticks = [2:4:32];
end

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

withinbar = [LDA_withinmean; SVM_withinmean; ANN_withinmean; SRC_withinmean];
withinbarerr = [LDA_withinstd; SVM_withinstd; ANN_withinstd; SRC_withinstd]/sqrt(10);

withinssbar = [LDA_withinssmean; SVM_withinssmean; ANN_withinssmean; SRC_withinssmean];
withinssbarerr = [LDA_withinssstd; SVM_withinssstd; ANN_withinssstd; SRC_withinssstd]/sqrt(10); 

withintbar = [LDA_withintmean; SVM_withintmean; ANN_withintmean; SRC_withintmean];
withintbarerr = [LDA_withintstd; SVM_withintstd; ANN_withintstd; SRC_withintstd]/sqrt(10);

% Create results matrices (for use with errorbar_groups) for single
% modalities plot and fused modalities plot
withinbar_single = reshape(withinbar(1,1:9),3,3);
withinbarerr_single = reshape(withinbarerr(1,1:9),3,3);
withinbar_fused = withinbar(1:3,10:12)';
withinbarerr_fused = withinbarerr(1:3,10:12)';

% Make bar plots
[~,~,singlebar_handles] = errorbar_groups(withinbar_single,withinbarerr_single,'bar_colors',[128 128 255; 0 0 255; 0 0 128; 128 255 128; 0 255 0; 0 128 0; 255 128 128; 255 0 0; 128 0 0]/255,'bar_names',{'EMG','GONIO','IMU'},'bar_width',0.75,'errorbar_width',1.5,'optional_bar_arguments',{'LineWidth',1.5},...
    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',2})
ylabel('Overall error (%)');
set(gca,'FontSize',18,'FontWeight','bold');
legend({'Ipsilateral','Contralateral','Bilateral'},'Location','northoutside','Orientation','horizontal')

[~,~,fusedbar_handles] = errorbar_groups(withinbar_fused,withinbarerr_fused,'bar_colors',[0.8 0.8 0.8; 0.5 0.5 0.5; 0.2 0.2 0.2; 0.8 0.8 0.8; 0.5 0.5 0.5; 0.2 0.2 0.2; 0.8 0.8 0.8; 0.5 0.5 0.5; 0.2 0.2 0.2;],'bar_names',{'LDA','SVM','ANN'},'bar_width',0.75,'errorbar_width',1.5,'optional_bar_arguments',{'LineWidth',1.5},...
    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',2})
ylabel('Overall error (%)');
set(gca,'FontSize',18,'FontWeight','bold');
legend({'Ipsilateral','Contralateral','Bilateral'},'Location','northoutside','Orientation','horizontal')
end

% Function computes overall, steady-state, and transitional error rates
% from predictions
function [errorrate, sserror, terror] = calcerrorrate(true,pred,ssinds,tinds)
errorrate = sum(true~=pred)/length(true);
sserror = sum(true(ssinds) ~= pred(ssinds))/length(ssinds);
terror = sum(true(tinds) ~= pred(tinds))/length(tinds);
end