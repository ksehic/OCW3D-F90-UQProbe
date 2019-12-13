%% Low-dimensional representation for the U-learning of a Gaussian process based on OceanWave3D
%{
---------------------------------------------------------------------------
Created by:
Kenan Šehić (kense@dtu.dk; kenosehic@gmail.com)
Department of Applied Mathematics and Computer Science
Technical University of Denmark
Licence: Copyright (C) 2019 Kenan Šehić DTU Compute, Technical University of Denmark

Cite: Šehić K., Bredmose H., Sørensen J.D., Karamehmedović M.: Low-dimensional representation of wave generation to quantify extreme events, TBD
---------------------------------------------------------------------------
Version December 2019
---------------------------------------------------------------------------
Description:
* Discrete not-optimal learning of GP to quantify exceedance probability
(failure probability) based on the U-function (the probability of
misclassification). Low-dimensional representation substitutes
high-dimensional initial uncertainties to defined initial waves for
low-dimensional classification factors, which are used then as design
parameters for learning.
---------------------------------------------------------------------------
%}
%% Define seed and filename
rng('shuffle')

sss = rng;

fn = sss.Seed;

filename = sprintf('%s_%d','OW3D_UQrun',fn);

%% Wave generation and model parameters

TDur=256; % time duration
dt=0.25; % time step

Hs=4; %significant wave height
Tp=9; % time period

fHighCut=0.3;
f0 = 0.05;

df=1/TDur; %freq step

f=f0:df:fHighCut; % freq

omega=2*pi*f; % ang freq

dd=length(omega); 

tSpan=0:dt:TDur; % time span

%% Generate random Fourier coefficients from standard normal distribution
input = randn(1,2*dd); %

on = 0; %% not making wave file

%% generate initial waves
N = 50000; % How many data points

Xbc = zeros(length(tSpan),N);

fprintf("Generating inital data...\n")

for i=1:N
    
    input = randn(1,2*dd);

    etabc = wavegen(input,TDur,dt,Hs,Tp,f0,fHighCut,on); % generate initial surface elevation based on Jonswap

    Xbc(:,i) = etabc;
    
    fprintf("Run %g of %g \n",i,N);

end

fprintf("Finished.\n")

%% Generate classification factors based on the initial data at the boundary

Xbc = Xbc.';

BCmax = max(Xbc.').'; % max
BCvar = zeros(N,1); % variance
BCsk = zeros(N,1); % skewness
BCku = zeros(N,1); % kurtosis
BCrms = zeros(N,1); % root mean squared
BCent = zeros(N,1); % approximateEntropy(X)
BCp50 = zeros(N,1); % theshold for 50%
BCp75 = zeros(N,1); % theshold for 75%
BCp90 = zeros(N,1); % theshold for 90%
BCmode = zeros(N,1); % most frequent value

for i=1:N % Possibility for parallel computing, some factor does not require for loop but can be used directly on matrix N x d
   
    BCvar(i,1) = var(Xbc(i,:));
    BCsk(i,1) = skewness(Xbc(i,:));
    BCku(i,1) = kurtosis(Xbc(i,:));
    BCrms(i,1) = rms(Xbc(i,:));
    
    BCent(i,1) = approximateEntropy(Xbc(i,:));
    
    BCp50(i,1) = prctile(Xbc(i,:),50);
    BCp75(i,1) = prctile(Xbc(i,:),75);
    BCp90(i,1) = prctile(Xbc(i,:),90);
    
    BCmode(i,1) = mode(Xbc(i,:));
    
    fprintf("Estimating Classification Factors -> Iteration: %g of %g\n",i,N);
    
end

BCdata = [BCmax,BCvar,BCsk,BCku,BCrms,BCent,BCp50,BCp75,BCp90,BCmode]; % Collect classification estimations within initial data

%% Prepare initial design

[BCdatas,id]=sort(BCdata(:,1),'descend');

BCdata=BCdata(id,:);
Xbc = Xbc(id,:);

[rr1,~]=find(BCdata(:,1)<=2.5);
[rr2,~]=find(BCdata(:,1)<=3 & BCdata(:,1)>2.5);
[rr3,~]=find(BCdata(:,1)<=3.5 & BCdata(:,1)>3);
[rr4,~]=find(BCdata(:,1)<=4 & BCdata(:,1)>3.5);
[rr5,~]=find(BCdata(:,1)>4);

BCdata1 = BCdata(rr1,:);
BCdata2 = BCdata(rr2,:);
BCdata3 = BCdata(rr3,:);
BCdata4 = BCdata(rr4,:);
BCdata5 = BCdata(rr5,:);

Xbc1 = Xbc(rr1,:);
Xbc2 = Xbc(rr2,:);
Xbc3 = Xbc(rr3,:);
Xbc4 = Xbc(rr4,:);
Xbc5 = Xbc(rr5,:);

ws1 = randperm(length(rr1),10);
ws2 = randperm(length(rr2),10);
ws3 = randperm(length(rr3),15);
ws4 = randperm(length(rr4),20);
ws5 = randperm(length(rr5),5);

BCdataX = [BCdata1(ws1,:);BCdata2(ws2,:);BCdata3(ws3,:);BCdata4(ws4,:);BCdata5(ws5,:)];
XbcX = [Xbc1(ws1,:);Xbc2(ws2,:);Xbc3(ws3,:);Xbc4(ws4,:);Xbc5(ws5,:)];

clear rr1 rr2 rr3 rr4 rr5
clear BCdata1 BCdata2 BCdata3 BCdata4 BCdata5
clear Xbc1 Xbc2 Xbc3 Xbc4 Xbc5

%% Data visualization

figure(1)
subplot(1,2,1)
scatter(BCdata(:,1)./max(BCdata(:,1)),sum((BCdata(:,2:end)./max(BCdata(:,2:end))).^2,2),20,'filled');
grid on
hold on
scatter(BCdataXu(1:60,1)./max(BCdata(:,1)),sum((BCdataXu(1:60,2:end)./max(BCdata(:,2:end))).^2,2),100,'filled');
xlabel({'K_1';'(a)'});
ylabel({'$\sum_{i=2}^{10} K^2_i$'},'Interpreter','latex')
set(gca,'fontsize', 28);
%legend('Generated samples','Initial Design points','U-Design points','Location','northwest','FontSize',22);
legend('Generated samples N','Initial Design points N_0','Location','northwest','FontSize',22);
xlim([0 1])
axis square
hold off

%% Run OceanWave3D

maxdata=zeros(length(BCdataX(:,1)),137);

for i=1:length(BCdataX(:,1)) % Possibility for parallel computing
%% create folder
mkdir run_case

%% Generate initial wave files

ini_wave(XbcX(i,:),dt);

%% copy OceanWave3D file in run_case folder
copyfile ReadKinematics.m run_case
copyfile BuildStencilEven.m run_case
copyfile initial_waves.iwf run_case
copyfile beach_137 run_case
copyfile OceanWave3D.inp run_case

%% Change directory
cd run_case

%% Run OceanWave3D
command = 'OceanWave3D_devel OceanWave3D.inp';

[status,s2] = system(command);

%% Collect data at wind turbine
eta_ref = ReadKinematics;

%% Pick max wave height at node 100th
eta_max = max(eta_ref);

maxdata(i,:) = eta_max;

%% Return to original folder
cd ..

%% Delete run_case
rmdir run_case s

%% Print
fprintf("%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf("%%%%%%%%%% Done sample %g of %g %%%%%%%%%%%%%% \n",i,length(BCdataX(:,1)));
fprintf("%%%%%%%%%%%%%%%%%%%%%%%% \n");
end

%% Failure threshold for the quantity of interest
fup = 6; %[m]
flow = 0; %[m]
fstep = 100; % exceedance discretization
tr = linspace(flow,fup,fstep); % failure

%% Define a location within OceanWave3D domain where is wind turbine

xloc = 73;

%% Fit GP
BCdataY = maxdata(:,xloc);

gppls = fitrgp(BCdataX,BCdataY,'KernelFunction','ardsquaredexponential','BasisFunction','pureQuadratic','FitMethod','exact','PredictMethod','exact');

%% Predict the initial data

[z,zs] = predict(gppls,BCdata);

%% Estimate the 95% CI

z1=z+1.96*zs;
z2=z-1.96*zs;

%% Calculate failure probabilities for tr

Pfgp = zeros(fstep,1);
Pfgp1 = zeros(fstep,1);
Pfgp2 = zeros(fstep,1);

for j=1:length(tr)

    Pfgp(j) = sum((z>=tr(j)))/length(z);
    Pfgp1(j) = sum((z1>=tr(j)))/length(z);
    Pfgp2(j) = sum((z2>=tr(j)))/length(z);  
    
end

%% Plot exceedance probabilities

figure(2)
subplot(1,2,1)
hold on
plot(tr,Pfgp,'-og','LineWidth',3);
plot(tr,Pfgp1,'--r','LineWidth',3);
plot(tr,Pfgp2,'--r','LineWidth',3);
grid on
xlabel({'\eta_{max} [m]','(a)'});
ylabel('P_F');
set(gca,'fontsize', 28);
set(gca, 'YScale', 'log');
xlim([0 tr(end)])
ylim([10^(-5) 1])
legend('\mu_{GP} with N_{GP} = 60','\mu_{GP} \pm \alpha \cdot \sigma_{GP}','Location','southwest','FontSize',22)
axis square
title('Initial Learning Form','FontWeight','normal','FontSize',28)

%% Prepaper the data for U-learning based on the U-function

Udata1 = BCdata1;
Udata2 = BCdata2;
Udata3 = BCdata3;
Udata4 = BCdata4;
Udata5 = BCdata5;

Udata1(ws1,:)=[]; % clear used samples
Udata2(ws2,:)=[];
Udata3(ws3,:)=[]; 
Udata4(ws4,:)=[];
Udata5(ws5,:)=[];

Udata = [Udata1;Udata2;Udata3;Udata4;Udata5]; % Learning data

clear Udata1 Udata2 Udata3 Udata4 Udata5

Uxbc1 = Xbc1;
Uxbc2 = Xbc2;
Uxbc3 = Xbc3;
Uxbc4 = Xbc4;
Uxbc5 = Xbc5;

Uxbc1(ws1,:)=[];
Uxbc2(ws2,:)=[];
Uxbc3(ws3,:)=[];
Uxbc4(ws4,:)=[];
Uxbc5(ws5,:)=[];

Uxbc = [Uxbc1;Uxbc2;Uxbc3;Uxbc4;Uxbc5]; % Learning initial surface elevations

clear Uxbc1 Uxbc2 Uxbc3 Uxbc4 Uxbc5

BCdataXu= BCdataX; % GP design set
BCdataYu = BCdataY; % GP correspoding respose
XbcXu = XbcX; % Correspoding initial surface elevations

%% Learning function U-function based on the probability of misclassification 
% Ref: R Schöbi, Bruno Sudret. Imprecise structural reliability analysis
% using PC-Kriging. 25th European Safety and Reliability Conference (ESREL2015), Sep 2015, Zurich, Switzerland. ff10.1201/b19094-549ff.ffhal-01247145f

UU=@(mm,ss,aa) abs(mm-aa)./ss;

%% Define failure/exceedance of interest

tru = 4.6667;
floc = 78; % location within the exceedance array in tr

%% Iteration
i=1;

check = abs(Pfgp1(floc)-Pfgp2(floc))/Pfgp(floc); % uncertainty measure based on the 95%CI

fprintf("Iteration %g with total num of samples %g and uncertainty measure %g\n",i,length(BCdataYu(:,1)),check);

i=2;

rj=1;

%% Define how many samples to include for learning
ss = 5; 

while true % Learning process
%% Predict and estimate U-values   
    
[zu,zsu] = predict(gppls,Udata); % Predict based on the learning data

UU1 = zeros(length(Udata(:,1)),length(tru));

for i1=1:length(tru)
    
    UU1(:,i1)=UU(zu,zsu,tru(i1)); % Estimate U-values based on the learning data

end

%% Defines samples with lowest U-values

xc=1:ss;

for i0 =1:length(tru)
   
    [us,id] = sort(UU1(:,i0));
    
    id1 = id(1:ss);
    
    BCdataXu0(xc,:) = Udata(id1,:);
    
    Uxbc0(xc,:) = Uxbc(id1,:);
    
    UU1(id1,:) = [];
    
    Udata(id1,:) = [];
    
    Uxbc(id1,:) = [];
    
    xc = ss + xc;
    
end

%% Run OceanWave3D for samples with lower U-values

for i1 = 1:length(BCdataXu0(:,1)) % Possibility for parallel computing
    
    outputOC = runOC(Uxbc0(i1,:),dt); % Run OceanWave3D
    
    BCdataYu0(i1,1) = outputOC(:,xloc); % Estimate the quantity of interest
    
    runsOC(rj,:) = outputOC; % Collect OceanWave3D Kinematics
    
    rj=rj+1;
    
end

%% Learn GP with new evaluations

BCdataXu = [BCdataXu;BCdataXu0];

BCdataYu = [BCdataYu;BCdataYu0];

XbcXu = [XbcXu;Uxbc0];

gppls = fitrgp(BCdataXu,BCdataYu,'KernelFunction','ardsquaredexponential','BasisFunction','pureQuadratic','FitMethod','exact','PredictMethod','exact'); % Fit GP

%% Estimate the initial data 
[z,zs] = predict(gppls,BCdata);

%% Estimate the 95%CI
z1=z+1.96*zs;
z2=z-1.96*zs;

%% Estimate failure/exceedance probabilities
Pfgp = zeros(fstep,1);
Pfgp1 = zeros(fstep,1);
Pfgp2 = zeros(fstep,1);

for j=1:length(tr)
    
    Pfgp(j) = sum((z>=tr(j)))/length(z);
    Pfgp1(j) = sum((z1>=tr(j)))/length(z);
    Pfgp2(j) = sum((z2>=tr(j)))/length(z);  
    
end

clear UU1 BCdataXu0 BCdataYu0 Uxbc0

%% Check uncertainty measure

check(i) = abs(Pfgp1(floc)-Pfgp2(floc))/Pfgp(floc);

fprintf("It %g with num of samples %g and dis %g\n",i,length(BCdataYu(:,1)),check(i));

if check(i) <= 2 && check(i)>=0 % Stopping Criterion for uncertainty measure
   
    break
    
end

i=i+1;

end

%% Plot when GP stops learning

figure(2)
subplot(1,2,2)
plot(tr,Pfgp,'-og','LineWidth',3);
hold on
plot(tr,Pfgp1,'--r','LineWidth',3);
plot(tr,Pfgp2,'--r','LineWidth',3);
grid on
xlabel({'\eta_{max} [m]','(a)'});
ylabel('P_F');
set(gca,'fontsize', 28);
set(gca, 'YScale', 'log');
xlim([0 tr(end)])
ylim([10^(-5) 1])
axis square
legend('\mu_{GP} with N_{GP} = 60','\mu_{GP} \pm \alpha \cdot \sigma_{GP}','Location','southwest','FontSize',28)
title('Final Learning Form','FontWeight','normal','FontSize',28)

%% Collect OceanWave3D Kinematic

OCdata = [maxdata;runsOC];

%% Global sensitivity Pearson correlation coefficient for spatial locations based on the OceanWave3D evalutions

Rq = zeros(10,137);

for j=1:137
    
    for i=1:10

        R = corrcoef(BCdataXu(:,i)./max(BCdataXu(:,i)),OCdata(:,j));

        Rq(i,j) = R(1,2);


    end
    
end

%% Sensitivity plots

subplot(2,2,1)
bar(1:10,Rq(:,18))
grid on
xlabel({'K_i';'(a)'});
ylabel('\rho');
xlim([1 10])
ylim([-1 1])
set(gca,'fontsize', 28);
axis square
xticks(1:10);
yticks(-1:0.5:1);
title("x=100m",'FontSize',24,'FontWeight','normal')

subplot(2,2,2)
bar(1:10,Rq(:,38))
grid on
xlabel({'K_i';'(b)'});
ylabel('\rho');
xlim([1 10])
ylim([-1 1])
set(gca,'fontsize', 28);
axis square
xticks(1:10);
yticks(-1:0.5:1);
title("x=217.64m",'FontSize',24,'FontWeight','normal')

subplot(2,2,3)
bar(1:10,Rq(:,73))
grid on
xlabel({'K_i';'(c)'});
ylabel('\rho');
xlim([1 10])
ylim([-1 1])
set(gca,'fontsize', 28);
axis square
xticks(1:10);
yticks(-1:0.5:1);
title("x=423.53m",'FontSize',24,'FontWeight','normal')

subplot(2,2,4)
bar(1:10,Rq(:,100))
grid on
xlabel({'K_i';'(d)'});
ylabel('\rho');
xlim([1 10])
ylim([-1 1])
set(gca,'fontsize', 28);
axis square
xticks(1:10);
yticks(-1:0.5:1);
title("x=582.35m",'FontSize',24,'FontWeight','normal')

%% Histograms

subplot(2,2,1)
histogram(OCdata(:,18),'Normalization','pdf')
grid on
xlabel({'\eta_{max} [m]';'(a)'});
ylabel('Pr');
xlim([2 8]);
ylim([0 1]);
set(gca,'fontsize', 28);
axis square
title("x=100m",'FontSize',24,'FontWeight','normal')

subplot(2,2,2)
histogram(OCdata(:,38),'Normalization','pdf')
grid on
xlabel({'\eta_{max} [m]';'(b)'});
ylabel('Pr');
xlim([2 8]);
ylim([0 1]);
set(gca,'fontsize', 28);
axis square
title("x=217.64m",'FontSize',24,'FontWeight','normal')

subplot(2,2,3)
histogram(OCdata(:,73),'Normalization','pdf')
grid on
xlabel({'\eta_{max} [m]';'(c)'});
ylabel('Pr');
xlim([2 8]);
ylim([0 1]);
set(gca,'fontsize', 28);
axis square
title("x=423.53m",'FontSize',24,'FontWeight','normal')

subplot(2,2,4)
histogram(OCdata(:,100),'Normalization','pdf')
grid on
xlabel({'\eta_{max} [m]';'(d)'});
ylabel('Pr');
xlim([2 8]);
ylim([0 1]);
set(gca,'fontsize', 28);
axis square
title("x=582.35m",'FontSize',24,'FontWeight','normal')

%% This is THE END...





















