function plottrajectories
% Function: Plot kinematic trajectories for each mode-specific classifier
% to confirm windowing

% Input: Use uipickfiles to select CircuitXXX_resegmented.mat (can select up to all for one subject)
% Output: NONE

% Function dependencies: 
% uipickfiles.m
% LPfilt.m

%%%%%
% Documented by: Blair Hu 08/04/17
%%%%%

close all
% Select up to all CircuitXXX_resegmented.mat files from a subject
files = uipickfiles;

% Specify the sampling rate and the maximum time after each gait event to
% plot the trajectory
Fs = 1000;
forwarddelay = 200;
saveall = {};
    
% Collect kinematic trajectories from each circuit
for i = 1:length(files)
    temp = load(files{i});
    RHC = temp.output_struct{2};
    RHC_TRIG = temp.output_struct{3};
    RTO = temp.output_struct{4};
    RTO_TRIG = temp.output_struct{5};
    LHC = temp.output_struct{8};
    LHC_TRIG = temp.output_struct{9};
    LTO = temp.output_struct{10};
    LTO_TRIG = temp.output_struct{11};
    GONIO = temp.output_struct{16};
    
    GONIO_LP = LPfilt(Fs,1,10,GONIO);
    
    MODE = temp.output_struct{17};
    
    % Right ankle (1)
    % Right knee (2)
    % Left ankle (3)
    % Left knee (4)
    
    for j = 1:length(RHC)
        windowtemp = RHC(j)-299:RHC(j)+forwarddelay;
        iknee = GONIO_LP(windowtemp,2);
        iankle = GONIO_LP(windowtemp,1);
        cknee = GONIO_LP(windowtemp,4);
        cankle = GONIO_LP(windowtemp,3);
        savetemp = GONIO_LP(windowtemp,[2 1 4 3]);
        trigtemp = RHC_TRIG(j);
        leavemode = str2num(trigtemp{1}(1));
        entermode = str2num(trigtemp{1}(3));
        modetemp = MODE(RHC(j)-299:RHC(j)+300);
        phasetemp = 1;
        saveall = [saveall; {savetemp trigtemp leavemode entermode modetemp phasetemp}];
    end
    
	for j = 1:length(RTO)
        windowtemp = RTO(j)-299:RTO(j)+forwarddelay;
        iknee = GONIO_LP(windowtemp,2);
        iankle = GONIO_LP(windowtemp,1);
        cknee = GONIO_LP(windowtemp,4);
        cankle = GONIO_LP(windowtemp,3);
        savetemp = GONIO_LP(windowtemp,[2 1 4 3]);
        trigtemp = RTO_TRIG(j);
        leavemode = str2num(trigtemp{1}(1));
        entermode = str2num(trigtemp{1}(3));
        modetemp = MODE(RTO(j)-299:RTO(j)+300);
        phasetemp = 2;
        saveall = [saveall; {savetemp trigtemp leavemode entermode modetemp phasetemp}];
    end
    
    for j = 1:length(LHC)
        windowtemp = LHC(j)-299:LHC(j)+forwarddelay;
        iknee = GONIO_LP(windowtemp,4);
        iankle = GONIO_LP(windowtemp,3);
        cknee = GONIO_LP(windowtemp,2);
        cankle = GONIO_LP(windowtemp,1);
        savetemp = GONIO_LP(windowtemp,[4 3 2 1]);
        trigtemp = LHC_TRIG(j);
        leavemode = str2num(trigtemp{1}(1));
        entermode = str2num(trigtemp{1}(3));
        modetemp = MODE(LHC(j)-299:LHC(j)+300);
        phasetemp = 3;
        saveall = [saveall; {savetemp trigtemp leavemode entermode modetemp phasetemp}];
    end
    
    for j = 1:length(LTO)
        windowtemp = LTO(j)-299:LTO(j)+forwarddelay;
        iknee = GONIO_LP(windowtemp,4);
        iankle = GONIO_LP(windowtemp,3);
        cknee = GONIO_LP(windowtemp,2);
        cankle = GONIO_LP(windowtemp,1);
        savetemp = GONIO_LP(windowtemp,[4 3 2 1]);
        trigtemp = LTO_TRIG(j);
        leavemode = str2num(trigtemp{1}(1));
        entermode = str2num(trigtemp{1}(3));
        modetemp = MODE(LTO(j)-299:LTO(j)+300);
        phasetemp = 4;
        saveall = [saveall; {savetemp trigtemp leavemode entermode modetemp phasetemp}];
    end
    
    clear temp
    
    % Temp cell array definitions:
    % Right heel contact (2)
    % Right heel contact trigger (3)
    % Right toe off (4)
    % Right toe off trigger (5)
    % Left heel contact (8)
    % Left heel contact trigger (9)
    % Left toe off (10)
    % Left toe off trigger (11)
    % IMU channels (14)
    % EMG channels (15)
    % GONIO channels (16)
    % MODE (17)
end

% Categorize each trajectory
LWLW_RHC = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 1);
LWLW_RTO = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 2);
LWLW_LHC = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 3);
LWLW_LTO = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 4);

LWRA_RHC = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 2 & cell2mat(saveall(:,6)) == 1);
LWRA_RTO = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 2 & cell2mat(saveall(:,6)) == 2);
LWRA_LHC = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 2 & cell2mat(saveall(:,6)) == 3);
LWRA_LTO = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 2 & cell2mat(saveall(:,6)) == 4);

LWRD_RHC = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 3 & cell2mat(saveall(:,6)) == 1);
LWRD_RTO = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 3 & cell2mat(saveall(:,6)) == 2);
LWRD_LHC = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 3 & cell2mat(saveall(:,6)) == 3);
LWRD_LTO = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 3 & cell2mat(saveall(:,6)) == 4);

LWSA_RHC = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 4 & cell2mat(saveall(:,6)) == 1);
LWSA_RTO = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 4 & cell2mat(saveall(:,6)) == 2);
LWSA_LHC = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 4 & cell2mat(saveall(:,6)) == 3);
LWSA_LTO = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 4 & cell2mat(saveall(:,6)) == 4);

LWSD_RHC = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 5 & cell2mat(saveall(:,6)) == 1);
LWSD_RTO = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 5 & cell2mat(saveall(:,6)) == 2);
LWSD_LHC = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 5 & cell2mat(saveall(:,6)) == 3);
LWSD_LTO = find(cell2mat(saveall(:,3)) == 1 & cell2mat(saveall(:,4)) == 5 & cell2mat(saveall(:,6)) == 4);

RARA_RHC = find(cell2mat(saveall(:,3)) == 2 & cell2mat(saveall(:,4)) == 2 & cell2mat(saveall(:,6)) == 1);
RARA_RTO = find(cell2mat(saveall(:,3)) == 2 & cell2mat(saveall(:,4)) == 2 & cell2mat(saveall(:,6)) == 2);
RARA_LHC = find(cell2mat(saveall(:,3)) == 2 & cell2mat(saveall(:,4)) == 2 & cell2mat(saveall(:,6)) == 3);
RARA_LTO = find(cell2mat(saveall(:,3)) == 2 & cell2mat(saveall(:,4)) == 2 & cell2mat(saveall(:,6)) == 4);

RALW_RHC = find(cell2mat(saveall(:,3)) == 2 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 1);
RALW_RTO = find(cell2mat(saveall(:,3)) == 2 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 2);
RALW_LHC = find(cell2mat(saveall(:,3)) == 2 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 3);
RALW_LTO = find(cell2mat(saveall(:,3)) == 2 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 4);

RDRD_RHC = find(cell2mat(saveall(:,3)) == 3 & cell2mat(saveall(:,4)) == 3 & cell2mat(saveall(:,6)) == 1);
RDRD_RTO = find(cell2mat(saveall(:,3)) == 3 & cell2mat(saveall(:,4)) == 3 & cell2mat(saveall(:,6)) == 2);
RDRD_LHC = find(cell2mat(saveall(:,3)) == 3 & cell2mat(saveall(:,4)) == 3 & cell2mat(saveall(:,6)) == 3);
RDRD_LTO = find(cell2mat(saveall(:,3)) == 3 & cell2mat(saveall(:,4)) == 3 & cell2mat(saveall(:,6)) == 4);

RDLW_RHC = find(cell2mat(saveall(:,3)) == 3 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 1);
RDLW_RTO = find(cell2mat(saveall(:,3)) == 3 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 2);
RDLW_LHC = find(cell2mat(saveall(:,3)) == 3 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 3);
RDLW_LTO = find(cell2mat(saveall(:,3)) == 3 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 4);

SASA_RHC = find(cell2mat(saveall(:,3)) == 4 & cell2mat(saveall(:,4)) == 4 & cell2mat(saveall(:,6)) == 1);
SASA_RTO = find(cell2mat(saveall(:,3)) == 4 & cell2mat(saveall(:,4)) == 4 & cell2mat(saveall(:,6)) == 2);
SASA_LHC = find(cell2mat(saveall(:,3)) == 4 & cell2mat(saveall(:,4)) == 4 & cell2mat(saveall(:,6)) == 3);
SASA_LTO = find(cell2mat(saveall(:,3)) == 4 & cell2mat(saveall(:,4)) == 4 & cell2mat(saveall(:,6)) == 4);

SALW_RHC = find(cell2mat(saveall(:,3)) == 4 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 1);
SALW_RTO = find(cell2mat(saveall(:,3)) == 4 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 2);
SALW_LHC = find(cell2mat(saveall(:,3)) == 4 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 3);
SALW_LTO = find(cell2mat(saveall(:,3)) == 4 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 4);

SDSD_RHC = find(cell2mat(saveall(:,3)) == 5 & cell2mat(saveall(:,4)) == 5 & cell2mat(saveall(:,6)) == 1);
SDSD_RTO = find(cell2mat(saveall(:,3)) == 5 & cell2mat(saveall(:,4)) == 5 & cell2mat(saveall(:,6)) == 2);
SDSD_LHC = find(cell2mat(saveall(:,3)) == 5 & cell2mat(saveall(:,4)) == 5 & cell2mat(saveall(:,6)) == 3);
SDSD_LTO = find(cell2mat(saveall(:,3)) == 5 & cell2mat(saveall(:,4)) == 5 & cell2mat(saveall(:,6)) == 4);

SDLW_RHC = find(cell2mat(saveall(:,3)) == 5 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 1);
SDLW_RTO = find(cell2mat(saveall(:,3)) == 5 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 2);
SDLW_LHC = find(cell2mat(saveall(:,3)) == 5 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 3);
SDLW_LTO = find(cell2mat(saveall(:,3)) == 5 & cell2mat(saveall(:,4)) == 1 & cell2mat(saveall(:,6)) == 4);

black = [0 0 0];
red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
mag = [1 0 1];

% For each leg/gait event, plot the mode-specific kinematic trajectories
%%%%% LW %%%%%
figure(1);
plotmode_phase(saveall,LWLW_RHC,black);
plotmode_phase(saveall,LWRA_RHC,red);
plotmode_phase(saveall,LWRD_RHC,green);
plotmode_phase(saveall,LWSA_RHC,blue);
plotmode_phase(saveall,LWSD_RHC,mag);

figure(2);
plotmode_phase(saveall,LWLW_RTO,black);
plotmode_phase(saveall,LWRA_RTO,red);
plotmode_phase(saveall,LWRD_RTO,green);
plotmode_phase(saveall,LWSA_RTO,blue);
plotmode_phase(saveall,LWSD_RTO,mag);

figure(3);
plotmode_phase(saveall,LWLW_LHC,black);
plotmode_phase(saveall,LWRA_LHC,red);
plotmode_phase(saveall,LWRD_LHC,green);
plotmode_phase(saveall,LWSA_LHC,blue);
plotmode_phase(saveall,LWSD_LHC,mag);

figure(4);
plotmode_phase(saveall,LWLW_LTO,black);
plotmode_phase(saveall,LWRA_LTO,red);
plotmode_phase(saveall,LWRD_LTO,green);
plotmode_phase(saveall,LWSA_LTO,blue);
plotmode_phase(saveall,LWSD_LTO,mag);

%%%%% RA %%%%%
figure(5);
plotmode_phase(saveall,RARA_RHC,red);
plotmode_phase(saveall,RALW_RHC,black);

figure(6);
plotmode_phase(saveall,RARA_RTO,red);
plotmode_phase(saveall,RALW_RTO,black);

figure(7);
plotmode_phase(saveall,RARA_LHC,red);
plotmode_phase(saveall,RALW_LHC,black);

figure(8);
plotmode_phase(saveall,RARA_LTO,red);
plotmode_phase(saveall,RALW_LTO,black);

%%%%% RD %%%%%
figure(9);
plotmode_phase(saveall,RDRD_RHC,green);
plotmode_phase(saveall,RDLW_RHC,black);

figure(10);
plotmode_phase(saveall,RDRD_RTO,green);
plotmode_phase(saveall,RDLW_RTO,black);

figure(11);
plotmode_phase(saveall,RDRD_LHC,green);
plotmode_phase(saveall,RDLW_LHC,black);

figure(12);
plotmode_phase(saveall,RDRD_LTO,green);
plotmode_phase(saveall,RDLW_LTO,black);

%%%%% SA %%%%%
figure(13);
plotmode_phase(saveall,SASA_RHC,blue);
plotmode_phase(saveall,SALW_RHC,black);

figure(14);
plotmode_phase(saveall,SASA_RTO,blue);
plotmode_phase(saveall,SALW_RTO,black);

figure(15);
plotmode_phase(saveall,SASA_LHC,blue);
plotmode_phase(saveall,SALW_LHC,black);

figure(16);
plotmode_phase(saveall,SASA_LTO,blue);
plotmode_phase(saveall,SALW_LTO,black);

%%%%% SD %%%%%
figure(17);
plotmode_phase(saveall,SDSD_RHC,mag);
plotmode_phase(saveall,SDLW_RHC,black);

figure(18);
plotmode_phase(saveall,SDSD_RTO,mag);
plotmode_phase(saveall,SDLW_RTO,black);

figure(19);
plotmode_phase(saveall,SDSD_LHC,mag);
plotmode_phase(saveall,SDLW_LHC,black);

figure(20);
plotmode_phase(saveall,SDSD_LTO,mag);
plotmode_phase(saveall,SDLW_LTO,black);
end

% Plots the kinematic trajectories for a given leg-gait-event-trigger
% combination with the specified color
function plotmode_phase(saveall,indices,color)
color_lt = color; color_lt(4) = 0.2;
toavg = cell(1,4);
for k = 1:length(indices)
    subplot(221)
    plot(saveall{indices(k),1}(:,1),'Color',color_lt); hold on;
    subplot(223)
    plot(saveall{indices(k),1}(:,2),'Color',color_lt); hold on;
    subplot(222)
    plot(saveall{indices(k),1}(:,3),'Color',color_lt); hold on;
    subplot(224)
    plot(saveall{indices(k),1}(:,4),'Color',color_lt); hold on;
    toavg{1} = [toavg{1} saveall{indices(k),1}(:,1)];
    toavg{2} = [toavg{2} saveall{indices(k),1}(:,2)];
    toavg{3} = [toavg{3} saveall{indices(k),1}(:,3)];
    toavg{4} = [toavg{4} saveall{indices(k),1}(:,4)];
    if k == length(indices)
        subplot(221)
        title('Knee (Ipsi)');
        ylabel('Angle(deg)');
        plot(mean(toavg{1},2),'Color',color,'LineWidth',3)
        set(gca,'XTick',[0:50:500])
        set(gca,'XTickLabel',{'-300','-250','-200','-150','-100','-50','0','50','100','150','200'})
        subplot(223)
        title('Ankle (Ipsi)');
        ylabel('Angle(deg)');
        plot(mean(toavg{2},2),'Color',color,'LineWidth',3)
        set(gca,'XTick',[0:50:500])
        set(gca,'XTickLabel',{'-300','-250','-200','-150','-100','-50','0','50','100','150','200'})
        subplot(222)
        title('Knee (Contra)');
        ylabel('Angle(deg)');
        plot(mean(toavg{3},2),'Color',color,'LineWidth',3)
        set(gca,'XTick',[0:50:500])
        set(gca,'XTickLabel',{'-300','-250','-200','-150','-100','-50','0','50','100','150','200'})
        subplot(224)
        title('Ankle (Contra)');
        ylabel('Angle(deg)');
        plot(mean(toavg{4},2),'Color',color,'LineWidth',3)
        set(gca,'XTick',[0:50:500])
        set(gca,'XTickLabel',{'-300','-250','-200','-150','-100','-50','0','50','100','150','200'})
    end
end
end