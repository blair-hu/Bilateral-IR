function plotfeats
% Function: Plot feature space for each mode-specific classifier

% Input: File selected using uipickfiles
% Output: NONE

% Function dependencies: 
% gscatter3.m

%%%%%
% Documented by: Blair Hu 08/04/17
%%%%%

% Load one file (ABXXX_toclassify.mat) 
files = uipickfiles;
load(files{1});

close all

for i = 1:length(trig)
leavemode(i,1) = str2num(trig{i}(1));
entermode(i,1) = str2num(trig{i}(3));
end

% Keep walking indices only
keepinds = find((leavemode > 0) & (leavemode < 6) & (entermode > 0) & (entermode < 6));

legphase = legphase(keepinds);
leavemode = leavemode(keepinds);
entermode = entermode(keepinds);
feats = feats{1}(keepinds,:);

% Specify which phase to plot
% RHC (1)
% RTO (2)
% LHC (3)
% LTO (4)
phase = 1;

lw_lw = find((legphase == phase) & (leavemode == 1) & (entermode == 1));
lw_ra = find((legphase == phase) & (leavemode == 1) & (entermode == 2));
lw_rd = find((legphase == phase) & (leavemode == 1) & (entermode == 3));
lw_sa = find((legphase == phase) & (leavemode == 1) & (entermode == 4));
lw_sd = find((legphase == phase) & (leavemode == 1) & (entermode == 5));

ra_ra = find((legphase == phase) & (leavemode == 2) & (entermode == 2));
ra_lw = find((legphase == phase) & (leavemode == 2) & (entermode == 1));

rd_rd = find((legphase == phase) & (leavemode == 3) & (entermode == 3));
rd_lw = find((legphase == phase) & (leavemode == 3) & (entermode == 1));

sa_sa = find((legphase == phase) & (leavemode == 4) & (entermode == 4));
sa_lw = find((legphase == phase) & (leavemode == 4) & (entermode == 1));

sd_sd = find((legphase == phase) & (leavemode == 5) & (entermode == 5));
sd_lw = find((legphase == phase) & (leavemode == 5) & (entermode == 1));

% Normalize the features
feats = zscore(feats,0,2);
% Perform dimensionality reduction on the features using PCA
[~,feats_pca] = pca(feats);

figure;
gscatter3(feats_pca([lw_lw; lw_ra; lw_rd; lw_sa; lw_sd],1:3),[repmat(1,length(lw_lw),1); repmat(2,length(lw_ra),1); repmat(3,length(lw_rd),1); repmat(4,length(lw_sa),1);repmat(5,length(lw_sd),1)]);

figure;
gscatter3(feats_pca([ra_ra; ra_lw],1:3),[repmat(2,length(ra_ra),1); repmat(1,length(ra_lw),1)],[1 0 0; 0 0 0]);

figure;
gscatter3(feats_pca([rd_rd; rd_lw],1:3),[repmat(3,length(rd_rd),1); repmat(1,length(rd_lw),1)],[0 1 0; 0 0 0]);

figure;
gscatter3(feats_pca([sa_sa; sa_lw],1:3),[repmat(4,length(sa_sa),1); repmat(1,length(sa_lw),1)],[0 0 1; 0 0 0]);

figure;
gscatter3(feats_pca([sd_sd; sd_lw],1:3),[repmat(5,length(sd_sd),1); repmat(1,length(sd_lw),1)],[1 1 0; 0 0 0]);
end