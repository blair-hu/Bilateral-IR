%%%%%
% Load packages needed for ANOVA and multiple comparisons
%%%%%

# install.packages('xlsx')
# install.packages("lme4",
   repos=c("http://lme4.r-forge.r-project.org/repos",
      getOption("repos")[["CRAN"]]))
# install.packages('ggpubr')
# install.packages('multcomp')
# install.packages('multcompView')
# install.packages('car')
# install.packages('contrast')
# install.packages('lmerTest')
# install.packages('lsmeans')
# install.packages('Rmisc')

setwd('C:/Users/bhu/Git/Bilateral_Intent_Recognition/Stats')

library('xlsx')
library('lme4')
library('ggpubr')
% library('multcomp')
% library('multcompView')
% library('car')
% library('contrast')
library('lmerTest')
library('lsmeans')
library('Rmisc')

%%%%%
% Read in the data from Excel
%%%%%
df <- read.xlsx('AllSubsErrors_ForStats.xlsx',sheetIndex = 1)
df_sub <- subset(df,Modality<5 & (Classifier == 'LDA' | Classifier == 'SVM' | Classifier == 'ANN') & Delay == 0)
df_sub <- df_sub[,c(2,3,4,7,11,12,13,14)]
str(df_sub)

% Define single modality logical masks
% Modality 
E_I_L <- subset(df,Modality == 1 & Laterality == 1 & Delay == 0 & Classifier == 'LDA')
E_C_L <- subset(df,Modality == 1 & Laterality == -1 & Delay == 0 & Classifier == 'LDA')
E_B_L <- subset(df,Modality == 1 & Laterality == 2 & Delay == 0 & Classifier == 'LDA')
G_I_L <- subset(df,Modality == 2 & Laterality == 1 & Delay == 0 & Classifier == 'LDA')
G_C_L <- subset(df,Modality == 2 & Laterality == -1 & Delay == 0 & Classifier == 'LDA')
G_B_L <- subset(df,Modality == 2 & Laterality == 2 & Delay == 0 & Classifier == 'LDA')
I_I_L <- subset(df,Modality == 3 & Laterality == 1 & Delay == 0 & Classifier == 'LDA')
I_C_L <- subset(df,Modality == 3 & Laterality == -1 & Delay == 0 & Classifier == 'LDA')
I_B_L <- subset(df,Modality == 3 & Laterality == 2 & Delay == 0 & Classifier == 'LDA')

A_I_L <- subset(df,Modality == 4 & Laterality == 1 & Delay == 0 & Classifier == 'LDA')
A_I_S <- subset(df,Modality == 4 & Laterality == 1 & Delay == 0 & Classifier == 'SVM')
A_I_A <- subset(df,Modality == 4 & Laterality == 1 & Delay == 0 & Classifier == 'ANN')
A_C_L <- subset(df,Modality == 4 & Laterality == -1 & Delay == 0 & Classifier == 'LDA')
A_C_S <- subset(df,Modality == 4 & Laterality == -1 & Delay == 0 & Classifier == 'SVM')
A_C_A <- subset(df,Modality == 4 & Laterality == -1 & Delay == 0 & Classifier == 'ANN')
A_B_L <- subset(df,Modality == 4 & Laterality == 2 & Delay == 0 & Classifier == 'LDA')
A_B_S <- subset(df,Modality == 4 & Laterality == 2 & Delay == 0 & Classifier == 'SVM')
A_B_A <- subset(df,Modality == 4 & Laterality == 2 & Delay == 0 & Classifier == 'ANN')

% Form groups to compare all combinations of sensors and classifiers
df_sub_test = rbind(E_I_L,E_C_L,E_B_L,G_I_L,G_C_L,G_B_L,I_I_L,I_C_L,I_B_L,A_I_L,A_I_S,A_I_A,A_C_L,A_C_S,A_C_A,A_B_L,A_B_S,A_B_A)
df_sub_test <- df_sub_test[,c(2,3,4,7,11,12,13,14)]
str(df_sub_test)

% Form groups to compare LDA classification only across EMG(I/B), GONIO(I/B), IMU(I/B), and ALL(I/B)
df_sub_lda = rbind(E_I_L,E_B_L,G_I_L,G_B_L,I_I_L,I_B_L,A_I_L,A_B_L)
df_sub_lda <- df_sub_lda[,c(2,3,4,7,11,12,13,14)]
str(df_sub_lda)

% Form groups to compare classifiers for ALL(I) and ALL(B) sensor sets only
df_sub_class = rbind(A_I_L,A_I_S,A_I_A,A_B_L,A_B_S,A_B_A)
df_sub_class <- df_sub_class[,c(2,3,4,7,11,12,13,14)]
str(df_sub_class)

% Form groups to compare classifiers for ALL(I), ALL(C), ALL(B) sensor sets only
df_sub_all = rbind(A_I_L,A_I_S,A_I_A,A_C_L,A_C_S,A_C_A,A_B_L,A_B_S,A_B_A)
df_sub_all <- df_sub_all[,c(2,3,4,7,11,12,13,14)]
str(df_sub_all)

%%%%%
% Convert modality and laterality keys into categorical factors
%%%%%
df_sub$Modality <- as.factor(df_sub$Modality)
df_sub$Laterality <- as.factor(df_sub$Laterality)

df_sub_test$Modality <- as.factor(df_sub_test$Modality)
df_sub_test$Laterality <- as.factor(df_sub_test$Laterality)

df_sub_lda$Modality <- as.factor(df_sub_lda$Modality)
df_sub_lda$Laterality <- as.factor(df_sub_lda$Laterality)

df_sub_class$Classifier <- as.factor(df_sub_class$Classifier)
df_sub_class$Laterality <- as.factor(df_sub_class$Laterality)

df_sub_all$Modality <- as.factor(df_sub_all$Modality)
df_sub_all$Laterality <- as.factor(df_sub_all$Laterality)

%%%%%
% Print summary statistics for each combination of fixed factors
%%%%%
sum = summarySE(df_sub,measurevar='Overall',groupvars=c('Classifier','Modality','Laterality'))
print(sum)

sum_SS = summarySE(df_sub,measurevar='SS',groupvars=c('Classifier','Modality','Laterality'))
print(sum_SS)

sum_T = summarySE(df_sub,measurevar='T',groupvars=c('Classifier','Modality','Laterality'))
print(sum_T)

%%%%%
% Boxplots of main effect and interaction
%%%%%
boxplot(Overall ~ Modality*Laterality*Classifier, data = df_sub, xlab = 'Modality x Laterality x Classifier',ylab = 'Overall error')

par(mfrow=c(1,3)) # put the next two plots side by side
plot(Overall ~ Modality, data=df_sub_lda, main="Modality")
plot(Overall ~ Laterality, data=df_sub_lda, main="Laterality")
plot(Overall ~ Classifier, data=df_sub_lda, main="Classifier")

ggboxplot(df_sub_lda,x='Modality',y='Overall',color='Classifier')
ggboxplot(df_sub_lda,x='Laterality',y='Overall',color='Classifier')

%%%%%
% Fit the linear models with random effect (subject) and conduct ANOVA
%%%%%
% (1) Determine effects of modality/laterality for LDA classification only
lme_model = lmer(Overall ~ Modality*Laterality + (1|Subject), data = df_sub_lda, REML=TRUE)
anova(lme_model)
summary(lme_model)
% Determine significance of random effect
rand(lme_model)

% (2) Determine effect of classifier on ALL(I) and ALL(B) sensor sets
lme_model_class = lmer(Overall ~ Classifier*Laterality + (1|Subject), data = df_sub_class, REML=TRUE)
anova(lme_model_class)
summary(lme_model_class)
% Determine significance of random effect
rand(lme_model_class)

%%%%%
% Check model assumptions
%%%%%
hist(residuals(lme_model),col='darkgray')
% The distribution of residuals should be approximately normal
plot(fitted(lme_model),residuals(lme_model))
% The residuals should be unbiased and evenly distributed around 0

%%%%%
% The interaction is significant so look at the simple effects (paired t-tests)
%%%%%
% Compare EMG Ipsi vs. Contra and Ipsi vs. Bilat
with(df_sub, t.test(Overall[Modality==1 & Classifier == 'LDA' & Laterality==1], Overall[Modality==1 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub, t.test(Overall[Modality==1 & Classifier == 'LDA' & Laterality==1], Overall[Modality==1 & Classifier == 'LDA' & Laterality==2],paired=TRUE))
with(df_sub, t.test(SS[Modality==1 & Classifier == 'LDA' & Laterality==1], SS[Modality==1 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub, t.test(SS[Modality==1 & Classifier == 'LDA' & Laterality==1], SS[Modality==1 & Classifier == 'LDA' & Laterality==2],paired=TRUE))
with(df_sub, t.test(T[Modality==1 & Classifier == 'LDA' & Laterality==1], T[Modality==1 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub, t.test(T[Modality==1 & Classifier == 'LDA' & Laterality==1], T[Modality==1 & Classifier == 'LDA' & Laterality==2],paired=TRUE))

% Compare GONIO Ipsi vs. Contra and Ipsi vs. Bilat
with(df_sub, t.test(Overall[Modality==2 & Classifier == 'LDA' & Laterality==1], Overall[Modality==2 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub, t.test(Overall[Modality==2 & Classifier == 'LDA' & Laterality==1], Overall[Modality==2 & Classifier == 'LDA' & Laterality==2],paired=TRUE))
with(df_sub, t.test(SS[Modality==2 & Classifier == 'LDA' & Laterality==1], SS[Modality==2 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub, t.test(SS[Modality==2 & Classifier == 'LDA' & Laterality==1], SS[Modality==2 & Classifier == 'LDA' & Laterality==2],paired=TRUE))
with(df_sub, t.test(T[Modality==2 & Classifier == 'LDA' & Laterality==1], T[Modality==2 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub, t.test(T[Modality==2 & Classifier == 'LDA' & Laterality==1], T[Modality==2 & Classifier == 'LDA' & Laterality==2],paired=TRUE))

% Compare IMU Ipsi vs. Contra and Ipsi vs. Bilat
with(df_sub, t.test(Overall[Modality==3 & Classifier == 'LDA' & Laterality==1], Overall[Modality==3 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub, t.test(Overall[Modality==3 & Classifier == 'LDA' & Laterality==1], Overall[Modality==3 & Classifier == 'LDA' & Laterality==2],paired=TRUE))
with(df_sub, t.test(SS[Modality==3 & Classifier == 'LDA' & Laterality==1], SS[Modality==3 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub, t.test(SS[Modality==3 & Classifier == 'LDA' & Laterality==1], SS[Modality==3 & Classifier == 'LDA' & Laterality==2],paired=TRUE))
with(df_sub, t.test(T[Modality==3 & Classifier == 'LDA' & Laterality==1], T[Modality==3 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub, t.test(T[Modality==3 & Classifier == 'LDA' & Laterality==1], T[Modality==3 & Classifier == 'LDA' & Laterality==2],paired=TRUE))

% Compare ALL Ipsi vs. Contra and Ipsi vs. Bilat
with(df_sub_all, t.test(Overall[Modality==4 & Classifier == 'LDA' & Laterality==1], Overall[Modality==4 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub_all, t.test(Overall[Modality==4 & Classifier == 'LDA' & Laterality==1], Overall[Modality==4 & Classifier == 'LDA' & Laterality==2],paired=TRUE))
with(df_sub_all, t.test(SS[Modality==4 & Classifier == 'LDA' & Laterality==1], SS[Modality==4 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub_all, t.test(SS[Modality==4 & Classifier == 'LDA' & Laterality==1], SS[Modality==4 & Classifier == 'LDA' & Laterality==2],paired=TRUE))
with(df_sub_all, t.test(T[Modality==4 & Classifier == 'LDA' & Laterality==1], T[Modality==4 & Classifier == 'LDA' & Laterality==-1],paired=TRUE))
with(df_sub_all, t.test(T[Modality==4 & Classifier == 'LDA' & Laterality==1], T[Modality==4 & Classifier == 'LDA' & Laterality==2],paired=TRUE))

with(df_sub_all, t.test(Overall[Modality==4 & Classifier == 'SVM' & Laterality==1], Overall[Modality==4 & Classifier == 'SVM' & Laterality==-1],paired=TRUE))
with(df_sub_all, t.test(Overall[Modality==4 & Classifier == 'SVM' & Laterality==1], Overall[Modality==4 & Classifier == 'SVM' & Laterality==2],paired=TRUE))
with(df_sub_all, t.test(SS[Modality==4 & Classifier == 'SVM' & Laterality==1], SS[Modality==4 & Classifier == 'SVM' & Laterality==-1],paired=TRUE))
with(df_sub_all, t.test(SS[Modality==4 & Classifier == 'SVM' & Laterality==1], SS[Modality==4 & Classifier == 'SVM' & Laterality==2],paired=TRUE))
with(df_sub_all, t.test(T[Modality==4 & Classifier == 'SVM' & Laterality==1], T[Modality==4 & Classifier == 'SVM' & Laterality==-1],paired=TRUE))
with(df_sub_all, t.test(T[Modality==4 & Classifier == 'SVM' & Laterality==1], T[Modality==4 & Classifier == 'SVM' & Laterality==2],paired=TRUE))

with(df_sub_all, t.test(Overall[Modality==4 & Classifier == 'ANN' & Laterality==1], Overall[Modality==4 & Classifier == 'ANN' & Laterality==-1],paired=TRUE))
with(df_sub_all, t.test(Overall[Modality==4 & Classifier == 'ANN' & Laterality==1], Overall[Modality==4 & Classifier == 'ANN' & Laterality==2],paired=TRUE))
with(df_sub_all, t.test(SS[Modality==4 & Classifier == 'ANN' & Laterality==1], SS[Modality==4 & Classifier == 'ANN' & Laterality==-1],paired=TRUE))
with(df_sub_all, t.test(SS[Modality==4 & Classifier == 'ANN' & Laterality==1], SS[Modality==4 & Classifier == 'ANN' & Laterality==2],paired=TRUE))
with(df_sub_all, t.test(T[Modality==4 & Classifier == 'ANN' & Laterality==1], T[Modality==4 & Classifier == 'ANN' & Laterality==-1],paired=TRUE))
with(df_sub_all, t.test(T[Modality==4 & Classifier == 'ANN' & Laterality==1], T[Modality==4 & Classifier == 'ANN' & Laterality==2],paired=TRUE))

% Compare laterality and classifiers for ALL(I) and ALL(B) sensor sets
dataALL <- subset(df_sub, Modality == 4)
all_lme_model = lmer(Overall ~ Laterality*Classifier + (1|Subject), data = dataALL, REML=TRUE)
anova(all_lme_model)
summary(all_lme_model)
% Determine significance of random effect
rand(all_lme_model)

with(df_sub_all, t.test(Overall[Classifier == 'LDA' & Laterality==1], Overall[Classifier == 'SVM' & Laterality==1],paired=TRUE))
with(df_sub_all, t.test(Overall[Classifier == 'LDA' & Laterality==1], Overall[Classifier == 'ANN' & Laterality==1],paired=TRUE))
with(df_sub_all, t.test(Overall[Classifier == 'LDA' & Laterality==2], Overall[Classifier == 'SVM' & Laterality==2],paired=TRUE))
with(df_sub_all, t.test(Overall[Classifier == 'LDA' & Laterality==2], Overall[Classifier == 'ANN' & Laterality==2],paired=TRUE))

%%%%%
% Stats for SFS
%%%%%
sfs <- read.xlsx('AllSubsSFS_ForStats.xlsx',sheetIndex = 1)
str(sfs)

sfs_lme_model = lmer(Overall ~ Iteration + (1|Subject), data = sfs, REML=TRUE)
anova(sfs_lme_model)
summary(sfs_lme_model)
% Determine significance of random effect
rand(sfs_lme_model)

% Check model assumptions
hist(residuals(sfs_lme_model),col='darkgray')
% The distribution of residuals should be approximately normal
plot(fitted(sfs_lme_model),residuals(sfs_lme_model))

% Compare baseline to adding one sensor
with(sfs, t.test(Overall[Iteration == 0], Overall[Iteration == 1],paired = TRUE))
with(sfs, t.test(SS[Iteration == 0], SS[Iteration == 1],paired = TRUE))
with(sfs, t.test(T[Iteration == 0], T[Iteration == 1],paired = TRUE))

%%%%%
% Coding examples:
%%%%%
https://rcompanion.org/rcompanion/d_08.html