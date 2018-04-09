function makerep
% Function: Make plots of representative EMG, GONIO, and IMU signals and
% raster plot of features across each modality (by concatenating to
% consecutive circuits)

% Input: Use uipickfiles to select CircuitXXX_resegmented.mat (two
% consecutive trials)
% Output: NONE

% Function dependencies: 
% getIMUfeats.m
% getEMGfeats.m
% getGONIOfeats.m
% b2r.m (used for the raster plot)

% Note: AB191 Circuit 11/12 was used for the publication figures
% Axes limits may need to be adjusted for other subjects/trials

%%%%%
% Documented by: Blair Hu 08/04/17
%%%%%

close all

set(0,'DefaultAxesFontName', 'Palatino Linotype')

files = uipickfiles;

IMU_sigs_combined = [];
EMG_sigs_combined = [];
GONIO_sigs_combined = [];
MODE_combined = [];

IMU_feats_combined = [];
EMG_feats_combined = [];
GONIO_feats_combined = [];

for j = 1:length(files)
    load(files{j});
    
    IMU_DAQ = output_struct{14};
    EMG_DAQ = output_struct{15};
    GONIO_DAQ = output_struct{16};
    MODE_DAQ = output_struct{17};
    
    % Keep walking indices only
    walkinds = intersect(find(MODE_DAQ > 0),find(MODE_DAQ < 6));
    
    EMG_DAQ = EMG_DAQ(walkinds,:);
    GONIO_DAQ = GONIO_DAQ(walkinds,:);
    % Ignore waist IMU signals
    IMU_DAQ = IMU_DAQ(walkinds,1:24);
    MODE_DAQ = MODE_DAQ(walkinds,:);
    
    % Define the sampling rate
    F_s = 1000; 
    
    % Pre-processing (filtering)
    EMG_DAQ_filt = HPfilt(1000,3,20,EMG_DAQ);
    EMG_DAQ_filt = NOTCHfilt(1000,3,[59.5 119.5 179.5],[60.5 120.5 180.5],EMG_DAQ_filt);
    
    GONIO_DAQ_filt = LPfilt(1000,3,10,GONIO_DAQ);
    GONIO_VEL = numderiv(GONIO_DAQ_filt,F_s,0);
    GONIO_DAQ_ADDVEL = [GONIO_DAQ_filt(:,1) GONIO_VEL(:,1) GONIO_DAQ_filt(:,2) GONIO_VEL(:,2) GONIO_DAQ_filt(:,3) GONIO_VEL(:,3) GONIO_DAQ_filt(:,4) GONIO_VEL(:,4)];
    
    IMU_DAQ_filt = LPfilt(1000,3,25,IMU_DAQ);
    
    % Define sliding windows for feature extraction
    window_inc = 30;
    window_length = 300;
    window_start = 1:window_inc:(length(MODE_DAQ)-window_length-1);

    IMU_feats_all = zeros(144,length(window_start));
    EMG_feats_all = zeros(140,length(window_start));
    GONIO_feats_all = zeros(48,length(window_start));

    window_count = 0;

    for i = window_start
       window_count = window_count + 1;
       temp_window = i:(i+window_length-1); 
       IMU_window = IMU_DAQ_filt(temp_window,:);
       EMG_window = EMG_DAQ_filt(temp_window,:);
       GONIO_window = GONIO_DAQ_ADDVEL(temp_window,:);
       [IMU_feats,IMU_labels] = getIMUfeats(IMU_window,1:24,1:6);
       [EMG_feats,EMG_labels] = getEMGfeats(EMG_window,1:14,1:10);
       [GONIO_feats,GONIO_labels] = getGONIOfeats(GONIO_window,1:8,1:6);

       IMU_feats_all(:,window_count) = IMU_feats';
       EMG_feats_all(:,window_count) = EMG_feats';
       GONIO_feats_all(:,window_count) = GONIO_feats';
    end
    
    IMU_feats_combined = [IMU_feats_combined IMU_feats_all];
    EMG_feats_combined = [EMG_feats_combined EMG_feats_all];
    GONIO_feats_combined = [GONIO_feats_combined GONIO_feats_all];
    
    IMU_sigs_combined = [IMU_sigs_combined; IMU_DAQ_filt];
    EMG_sigs_combined = [EMG_sigs_combined; EMG_DAQ_filt];
    GONIO_sigs_combined = [GONIO_sigs_combined; GONIO_DAQ_ADDVEL];
    MODE_combined = [MODE_combined; MODE_DAQ];
end

IMU_z_score = zscore(IMU_feats_combined,0,2);
EMG_z_score = zscore(EMG_feats_combined,0,2);
GONIO_z_score = zscore(GONIO_feats_combined,0,2);

% Rearrange by channel, not feature
[IMU_z_score_rearrange,~] = rearrangefeats(IMU_z_score,IMU_labels,24);
[EMG_z_score_rearrange,~] = rearrangefeats(EMG_z_score,EMG_labels,14);
[GONIO_z_score_rearrange,~] = rearrangefeats(GONIO_z_score,GONIO_labels,8);

% Plot raw signals (right = turquoise, left = purple)
% Define transparency to be 70%
turq = [0 153 153]/255;
turqlt = turq; turqlt(4) = 0.25;
turq(4) = 0.7;
purp = [102 0 204]/255;
purplt = purp; purplt(4) = 0.25;
purp(4) = 0.7;

% Specify whether to make representative signals plots
plotsigs = 1;

if plotsigs
    figure('position',[0,0,1365,300])
    plot(EMG_sigs_combined(:,1),'Color',turqlt); hold on; plot(EMG_sigs_combined(:,8),'Color',purplt);
    ylabel(['TA',char(10),'(V)'],'FontWeight','bold');
    axis([0 length(EMG_sigs_combined) -0.55 0.55])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'FontSize',18)
    
    figure('position',[0,0,1365,300])
    plot(EMG_sigs_combined(:,2),'Color',turqlt); hold on; plot(EMG_sigs_combined(:,9),'Color',purplt);
    ylabel(['MG',char(10),'(V)'],'FontWeight','bold');
    axis([0 length(EMG_sigs_combined) -0.55 0.55])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'FontSize',18)
    
    figure('position',[0,0,1365,300])
    plot(EMG_sigs_combined(:,5),'Color',turqlt); hold on; plot(EMG_sigs_combined(:,12),'Color',purplt);
    ylabel(['ST',char(10),'(V)'],'FontWeight','bold');
    axis([0 length(EMG_sigs_combined) -0.35 0.35])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'FontSize',18)
    
    figure('position',[0,0,1365,300])
    plot(EMG_sigs_combined(:,6),'Color',turqlt); hold on; plot(EMG_sigs_combined(:,13),'Color',purplt);
    ylabel(['VL',char(10),'(V)'],'FontWeight','bold');
    axis([0 length(EMG_sigs_combined) -0.2 0.2])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'FontSize',18)
    
    % Ankle
    figure('position',[0,0,1365,300])
    plot(GONIO_sigs_combined(:,1),'Color',turq,'LineWidth',1.5); hold on; plot(GONIO_sigs_combined(:,5),'Color',purp,'LineWidth',1.5);
    ylabel(['Ankle',char(10),'(deg)'],'FontWeight','bold');
    axis([0 length(GONIO_sigs_combined) -25 25])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'FontSize',18)
     
    % Knee
    figure('position',[0,0,1365,300])
    plot(GONIO_sigs_combined(:,3),'Color',turq,'LineWidth',1.5); hold on; plot(GONIO_sigs_combined(:,7),'Color',purp,'LineWidth',1.5);
    ylabel(['Knee',char(10),'(deg)'],'FontWeight','bold');
    axis([0 length(GONIO_sigs_combined) -105 15])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'FontSize',18)
    
	% Ankle Vel
    figure('position',[0,0,1365,300])
    plot(GONIO_sigs_combined(:,2),'Color',turq,'LineWidth',1.5); hold on; plot(GONIO_sigs_combined(:,6),'Color',purp,'LineWidth',1.5);
    ylabel(['Ankle Velocity',char(10),'(deg/s)'],'FontWeight','bold');
    axis([0 length(GONIO_sigs_combined) -250 250])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'FontSize',18)
     
    % Knee Vel
    figure('position',[0,0,1365,300])
    plot(GONIO_sigs_combined(:,4),'Color',turq,'LineWidth',1.5); hold on; plot(GONIO_sigs_combined(:,8),'Color',purp,'LineWidth',1.5);
    ylabel(['Knee Velocity',char(10),'(deg/s)'],'FontWeight','bold');
    axis([0 length(GONIO_sigs_combined) -500 500])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'FontSize',18)
    
    % Shank Acc
    figure('position',[0,0,1365,450])
    subplot(311)
    plot(IMU_sigs_combined(:,1),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,13),'Color',purp,'LineWidth',1.5);
    ylabel(['Shank',char(10),'A_X'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -3 1.5])
    set(gca,'FontSize',18)
    
    subplot(312)
    plot(IMU_sigs_combined(:,2),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,14),'Color',purp,'LineWidth',1.5);
    ylabel(['Shank',char(10),'A_Y'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -1.5 1.5])
    set(gca,'FontSize',18)
    
    subplot(313)
    plot(IMU_sigs_combined(:,3),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,15),'Color',purp,'LineWidth',1.5);
    ylabel(['Shank',char(10),'A_Z'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -2 2])
    set(gca,'FontSize',18)
    
    % Shank Vel
    figure('position',[0,0,1365,450])
    subplot(312)
    plot(IMU_sigs_combined(:,4),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,16),'Color',purp,'LineWidth',1.5);
    ylabel(['Shank',char(10),'G_Y'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -300 500])
    set(gca,'FontSize',18)
    
    subplot(313)
    plot(IMU_sigs_combined(:,5),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,17),'Color',purp,'LineWidth',1.5);
    ylabel(['Shank',char(10),'G_Z'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -250 250])
    set(gca,'FontSize',18)
    
    subplot(311)
    plot(IMU_sigs_combined(:,6),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,18),'Color',purp,'LineWidth',1.5);
    ylabel(['Shank',char(10),'G_X'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -350 350])
    set(gca,'FontSize',18)
    
    % Thigh Acc
    figure('position',[0,0,1365,450])
    subplot(311)
    plot(IMU_sigs_combined(:,7),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,19),'Color',purp,'LineWidth',1.5);
    ylabel(['Thigh',char(10),'A_X'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -3 1])
    set(gca,'FontSize',18)
    
    subplot(312)
    plot(IMU_sigs_combined(:,8),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,20),'Color',purp,'LineWidth',1.5);
    ylabel(['Thigh',char(10),'A_Y'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -1.5 1.5])
    set(gca,'FontSize',18)
    
    subplot(313)
    plot(IMU_sigs_combined(:,9),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,21),'Color',purp,'LineWidth',1.5);
    ylabel(['Thigh',char(10),'A_Z'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -1.75 1.75])
    set(gca,'FontSize',18)
    
    % Thigh Vel
    figure('position',[0,0,1365,450])
    subplot(312)
    plot(IMU_sigs_combined(:,10),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,22),'Color',purp,'LineWidth',1.5);
    ylabel(['Thigh',char(10),'G_Y'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -350 350])
    set(gca,'FontSize',18)
    
    subplot(313)
    plot(IMU_sigs_combined(:,11),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,23),'Color',purp,'LineWidth',1.5);
    ylabel(['Thigh',char(10),'G_Z'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -150 150])
    set(gca,'FontSize',18)
    
    subplot(311)
    plot(IMU_sigs_combined(:,12),'Color',turq,'LineWidth',1.5); hold on; plot(IMU_sigs_combined(:,24),'Color',purp,'LineWidth',1.5);
    ylabel(['Thigh',char(10),'G_X'],'FontWeight','bold');
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis([0 length(IMU_sigs_combined) -350 350])
    set(gca,'FontSize',18)
end

% Mode
figure('position',[0,0,1365,150])
plot(MODE_combined,'k','LineWidth',2);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[1 2 3 4 5])
set(gca,'yticklabel',{'LW','RA','RD','SA','SD'},'FontSize',16,'FontWeight','bold')
axis([0 length(MODE_combined) 0.5 8])

% Specify whether to make raster plot
plotraster = 1;

% Here we are only plotting the mean value feature (which corresponds to
% the first feature for each channel)
if plotraster
    % Plot raster plots
    figure()
    imagesc(IMU_z_score_rearrange(1:6:144,:)); 
    colormap(b2r(-4,4));
    caxis([-4 4]);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[1:6:24])
    set(gca,'yticklabel',{'Shank','Thigh','Shank','Thigh'},'FontSize',16,'FontWeight','bold')
    cb = colorbar;
    cb.Label.String = 'Z-score';  
    cb.Location = 'southoutside';
    
    figure()
    imagesc(GONIO_z_score_rearrange(1:6:48,:));
    colormap(b2r(-4,4));
    caxis([-4 4]);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[1:2:8])
    set(gca,'yticklabel',{'Ankle','Knee','Ankle','Knee'},'FontSize',16,'FontWeight','bold')

    figure()
    imagesc(EMG_z_score_rearrange(1:10:140,:));
    colormap(b2r(-4,4));    
    caxis([-4 4]);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[1:14])
    set(gca,'yticklabel',{'TA','MG','SOL','BF','ST','VL','RF','TA','MG','SOL','BF','ST','VL','RF'},'FontSize',16,'FontWeight','bold')
end
end

% Used to order rows of raster plot by feature or by channel
function [zscores_new,labels_new] = rearrangefeats(zscores,labels,numchans)
zscores_new = [];
labels_new = [];
for i = 1:numchans
    chaninds = i:numchans:size(zscores,1);
    zscores_new = [zscores_new; zscores(chaninds,:)];
    labels_new = [labels_new; labels(chaninds)];
end
end