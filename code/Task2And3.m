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
% filter settings
maxError=10;
stateChange=20;
nParticle=1000;
gridxl=100;
gridyl=100;
%% the extended kalman filter
% generate initial state
xState(1)=randn();   
%xState(1)=0;
% generate initial estimates and weights
xEstimate(1)=xState(1);
xEstimate2(1,:)=xState(1)*ones(1,nParticle)+randn(1,nParticle);
wEstimate2(1,:)=ones(1,nParticle)/nParticle;
xEstimate3(1,:)=linspace(-20,20,gridxl);
wEstimate3(1,:)=normpdf(xEstimate3(1,:),xState,sqrt(Q));
wEstimate3(1,:)=wEstimate3(1,:)./sum(wEstimate3(1,:));
% generate initial covaiance
PEstimate=1;  
zMeasure(1)=0;
%% states simulation and estimation
for i=2:Jn
   %% states simulation
   % generate new state
   trueState(i)=f_Func(xState(i-1),i);
   xState(i)=trueState(i)+sqrt(Q)*randn;  
   % perform measurement
   zMeasure(i)=g_Func(xState(i))+sqrt(R)*randn;
   %% states estimation
   % extended kalman filtering 
   [xEstimate(i), PEstimate(i)] = ExtendKalman(xEstimate(i-1), PEstimate(i-1),zMeasure(i),R,Q,i); 
   % grid-based filtering 
   [xEstimate3(i,:),wEstimate3(i,:),estimate3(i)] = MaxGridBasedFilter(xEstimate3(i-1,:),wEstimate3(i-1,:),zMeasure(i),Q,R,i,gridxl,gridyl);     
   % particle filtering 
   [xEstimate2(i,:),wEstimate2(i,:),estimate(i)] = GenericParticleFilter(xEstimate2(i-1,:),zMeasure(i),Q,R,i);     

end
%% figures of true states and estimations from different filters
% Comparison of real state and predicted state by extended Kalman filter
 figure()
 plot(trueState(1,:))
 hold on
 plot(xEstimate(1,:))
 title('Comparison of real state and predicted state by extended Kalman filter','FontSize',14)
 legend('true state' ,'extended kalman filter')
 xlabel('steps','FontSize',14)
 ylabel('states','FontSize',14)
 MSE=mean((xEstimate-trueState).^2);
 txt1 = ['RMSE: ',num2str(sqrt(MSE))];
 text(71,-25,txt1,'FontSize',14)
 hold off
% find max absolute error
max(xEstimate-trueState)
 
 % Comparison of real state and predicted state by grid-based filter
 figure()
 plot(trueState(1,:))
 hold on
 plot(estimate3(1,:))
 title('Comparison of real state and predicted state by grid-based filter','FontSize',14)
 legend('true state' ,'grid-based filter')
 xlabel('steps','FontSize',14)
 ylabel('states','FontSize',14)
 MSE4=mean((estimate3-trueState).^2);
 txt1 = ['RMSE: ',num2str(sqrt(MSE4))];
 text(83,-10,txt1,'FontSize',14)
 hold off
% find max absolute error
max(estimate3-trueState)


% Comparison of real state and predicted state by particle filter
 figure()
 plot(trueState(1,:))
 hold on
 plot(estimate(1,:))
 legend('true state','particle filter')
 title('Comparison of real state and predicted state by particle filter','FontSize',14)
 xlabel('steps','FontSize',14)
 ylabel('states','FontSize',14)
 MSE2=mean((estimate-trueState).^2);
 txt1 = ['RMSE: ',num2str(sqrt(MSE2))];
 text(83,-10,txt1,'FontSize',14)
 hold off
% find max absolute error
max(estimate-trueState)

%% comparison of different filters
 figure()
 boxplot([abs(xEstimate(1,:)-trueState(1,:))',abs(estimate3(1,:)-trueState(1,:))',abs(estimate(1,:)-trueState(1,:))'],'Notch','on','Labels',{'extended kalman','grid-based method','particle filter'})
 %boxplot([((xEstimate-trueState).^2)',((estimate3-trueState).^2)',((estimate-trueState).^2)'],'Notch','on','Labels',{'extended kalman','grid-based method','particle filter'})
 title('Comparison of different filters','FontSize',14)
 xlabel('Different Strategies','FontSize',14)
 ylabel('Absolute Error','FontSize',14)