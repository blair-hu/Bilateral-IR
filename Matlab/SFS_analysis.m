function SFS_analysis
% Function: Combine results of SFS for each subject and make results graph

% Input: NONE
% Output: Order of channels selected and their corresponding overall error
% rates (for each iteration)
% If plotonly = 1, then make double axis results figure

% Function dependencies:
% SFS_BHH.m

%%%%%
% Documented by: Blair Hu 08/04/17
%%%%%

close all

set(0,'DefaultAxesFontName', 'Palatino Linotype')

% Specify whether you want to just plot or do entire SFS
plotonly = 1;

if plotonly == 0
    % Use uipickfiles to select "ABXXX_toreclassify.mat"
    files = uipickfiles;

    % Iterate across the number of subjects
    for i = 1:length(files)
        disp(['Subject: ',num2str(i)])

        % Perform sequential forward selection for CONTRA sensors
        
        % Start with IPSI ALL, can select ANY CONTRA
        [contra_start_err,contra_start_ss_err,contra_start_t_err,contra_chan_include,contra_chan_selected,contra_error_include,contra_error_ss_include,contra_error_t_include] = SFS_BHH(files{i},[1:2 5:11 19:20],[3:4 12:18 21:22],11);

        % Record each subject's channel selection order
        contra_fs(:,i) = contra_chan_include;

        % Record each subject's overall error rates in a separate column        
        contra_err(:,i) = 100*[contra_start_err; contra_error_include];
        contra_ss_err(:,i) = 100*[contra_start_ss_err; contra_error_ss_include];
        contra_t_err(:,i) = 100*[contra_start_t_err; contra_error_t_include];
    end
    save('SFS_results_redo_contra.mat','contra_err','contra_ss_err','contra_t_err','contra_fs');
else
    load('SFS_results_redo_contra.mat')
    % Compute metrics to plot
    contra_mean = mean(contra_err,2);
    contra_sem = std(contra_err,0,2)/sqrt(10);
    contra_ss_mean = mean(contra_ss_err,2);
    contra_ss_sem = std(contra_ss_err,0,2)/sqrt(10);
    contra_t_mean = mean(contra_t_err,2);
    contra_t_sem = std(contra_t_err,0,2)/sqrt(10);
    
    % Determine the proportion of features of each modality for each iteration
    save_addedfeats = zeros(10,3);
    for i = 1:10 % iterate through each subject
        for j = 1:10 % iterate through each selection cycle
            tempadded = 1:j;
            IMU = 0;
            EMG = 0;
            GONIO = 0;
            for k = 1:length(tempadded)            
                if ismember(contra_fs(k,i),[21 22])
                    IMU = IMU + 36;
                elseif ismember(contra_fs(k,i),[3 4])
                    GONIO = GONIO + 12;
                elseif ismember(contra_fs(k,i),[12:18])
                    EMG = EMG + 10;
                end
            end
            addedfeats(j,:) = [EMG GONIO IMU];
        end
        subj_addedfeats{i} = addedfeats;
        save_addedfeats = save_addedfeats + addedfeats;
    end

    avg_addedfeats = save_addedfeats/10;
    avg_addedfeats = repmat([70 24 72],10,1) + avg_addedfeats;
    avg_addedfeats = [70 24 72; avg_addedfeats; 140 48 144]./repmat([140 48 144],12,1);

    figure;
    featbar = bar([1:12],100*avg_addedfeats)
    featbar(1).FaceColor = [0 0 128]/255;
    featbar(2).FaceColor = [0 128 0]/255;
    featbar(3).FaceColor = [128 0 0]/255;
    set(gca,'YTick',[50:10:100])
    set(gca,'YLim',[50 100])
    set(gca,'XTick',[1:12])
    set(gca,'XTickLabel',{'ALL(I)','+1','+2','+3','+4','+5','+6','+7','+8','+9','+10','ALL(B)'})
    set(gca,'FontSize',18,'FontWeight','bold');
    l{1}='EMG'; 
    l{2}='GONIO'; 
    l{3}='IMU'; 
    legend(featbar,l,'Location','northoutside','Orientation','horizontal');
    ylabel('Features added (%)')
    box(gca,'off')
    
    figure;
    plot([1:12],contra_mean,'k.','MarkerSize',24);
    hold on;
    errorbar([1:12],contra_mean,contra_sem,'LineWidth',2,'Color','k');
    set(gca,'YTick',[1:0.5:2.5])
    set(gca,'YLim',[0.75 2.75])
    set(gca,'XTick',[1:12])
    set(gca,'XTickLabel',{'ALL(I)','+1','+2','+3','+4','+5','+6','+7','+8','+9','+10','ALL(B)'})
    set(gca,'FontSize',18,'FontWeight','bold');
    ylabel('Overall error (%)')
    box(gca,'off')    
end
end