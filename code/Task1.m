%% initialize program
clear all
close all
clc
%% data inpute
% number of iterations
Kn=1;
% number of steps
Jn=101;
% covariance of noise
Q=10;
R=1;
%% the Monte Carlo simulations
%% 1000 times experiments
xState(:,1)=randn(1000,1); 
 for i=2:Jn
   % generate new state
   trueState(:,i)=f_Func(xState(:,i-1),i);
   xState(:,i)=trueState(:,i)+sqrt(Q)*randn(1000,1);  
 end
%% plot figures of states x0,x1,x50,x100

figure()
histogram(xState(:,1))
mea=mean(xState(:,1))
stand=std(xState(:,1))
title('distribution of states at K=0','FontSize',14)
xlabel('values of states','FontSize',14)
ylabel('number of instances','FontSize',14)
txt1 = ['mean: ',num2str(mea)];
txt2 = ['standard deviation: ',num2str(stand)];
%text(1,110,txt1,'FontSize',14)
%text(1,100,txt2,'FontSize',14)

% state at k=1
figure()
histogram(xState(:,2))
mea=mean(xState(:,2))
stand=std(xState(:,2))
title('distribution of states at K=1','FontSize',14)
xlabel('values of states','FontSize',14)
ylabel('number of instances','FontSize',14)
txt1 = ['mean: ',num2str(mea)];
txt2 = ['standard deviation: ',num2str(stand)];
%text(1,110,txt1,'FontSize',14)
%text(1,100,txt2,'FontSize',14)
GMModel = fitgmdist(xState(:,2),2)

% state at k=50
figure()
histogram(xState(:,51))
mea=mean(xState(:,51))
stand=std(xState(:,51))
title('distribution of states at K=50','FontSize',14)
xlabel('values of states','FontSize',14)
ylabel('number of instances','FontSize',14)
txt1 = ['mean: ',num2str(mea)];
txt2 = ['standard deviation: ',num2str(stand)];
%text(5,110,txt1,'FontSize',14)
%text(5,120,txt2,'FontSize',14)
GMModel = fitgmdist(xState(:,51),2)

% state at k=100
figure()
histogram(xState(:,101))
mea=mean(xState(:,101))
stand=std(xState(:,101))
title('distribution of states at K=100','FontSize',14)
xlabel('values of states','FontSize',14)
ylabel('number of instances','FontSize',14)
txt1 = ['mean: ',num2str(mea)];
txt2 = ['standard deviation: ',num2str(stand)];
%text(5,200,txt1,'FontSize',14)
%text(5,210,txt2,'FontSize',14)
GMModel = fitgmdist(xState(:,101),2)