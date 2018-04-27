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
% initilize RMSE 
 MSEsEKF =0;
 MSEsEKF2=0;
 MSEsGB =0;
 MSEsGB2=0;
 MSEsPF =0;
 MSEsPF2=0;
% filter settings
maxError=10;
stateChange=20;
nParticle=100;
gridxl=100;
gridyl=100;
Test=20;
%% the extended kalman filter
% generate initial state
xState(1)=randn(); 
zMeasure(1)=randn(); 
%xState(1)=0;
% initial state extended kalman filtering 
xEstimate(1)=xState(1);
% initial state iterated extended kalman filtering 
xEstimate2(1)=xState(1);
% initial state grid-based filtering 
xEstimate3(1,:)=linspace(-20,20,gridxl);
xEstimate32(1,:)=xEstimate3(1,:);
% initial state particle filtering 
xEstimate4(1,:)=xState(1)*ones(1,nParticle)+randn(1,nParticle);
% initial state generic particle filtering 
xEstimate5(1,:)=xState(1)*ones(1,nParticle)+randn(1,nParticle);

% initial weights grid-based filtering  
wEstimate3(1,:)=normpdf(xEstimate3(1,:),xState,sqrt(Q));
wEstimate32(1,:)=wEstimate3(1,:);
wEstimate3(1,:)=wEstimate3(1,:)./sum(wEstimate3(1,:));
% initial weights particle filtering  
wEstimate4(1,:)=ones(1,nParticle)/nParticle;
wEstimate5(1,:)=ones(1,nParticle)/nParticle;

% generate initial covaiance
PEstimate=1;  
PEstimate2=1;  
for j=1:Test
%% j times experiments
for i=2:Jn
   % generate new state
   trueState(i)=f_Func(xState(i-1),i);
   xState(i)=trueState(i)+sqrt(Q)*randn;  
   % perform measurement
   zMeasure(i)=g_Func(xState(i))+sqrt(R)*randn;
   %% kalman
   % extended kalman filtering 
   [xEstimate(i), PEstimate(i)] = ExtendKalman(xEstimate(i-1), PEstimate(i-1),zMeasure(i),R,Q,i); 
   % iterated extended kalman filtering 
   [xEstimate2(i), PEstimate2(i)] = IteratedExtendKalman(xEstimate2(i-1), PEstimate2(i-1),zMeasure(i),R,Q,i); 
   
   %% grid-based filtering 
   % grid-based filtering 
    [xEstimate3(i,:),wEstimate3(i,:),estimate3(i)] = GridBasedFilter(xEstimate3(i-1,:),wEstimate3(i-1,:),zMeasure(i),Q,R,i,gridxl,gridyl);
    [xEstimate32(i,:),wEstimate32(i,:),estimate32(i)] = MaxGridBasedFilter(xEstimate32(i-1,:),wEstimate32(i-1,:),zMeasure(i),Q,R,i,gridxl,gridyl); 
   %% particles 
   % SIR particle filtering 
    [xEstimate4(i,:),wEstimate4(i,:),estimate(i)] = SIRParticleFilter(xEstimate4(i-1,:),zMeasure(i),Q,R,i);     
   % generic particle filtering 
    [xEstimate5(i,:),wEstimate5(i,:),estimate2(i)] = GenericParticleFilter(xEstimate5(i-1,:),zMeasure(i),Q,R,i);    
end
 MSEsEKF =sqrt(mean((xEstimate-trueState).^2))+MSEsEKF;
 MSEsEKF2=sqrt(mean((xEstimate2-trueState).^2))+MSEsEKF2;
 MSEsGB =sqrt(mean((estimate3-trueState).^2))+MSEsGB;
 MSEsGB2=sqrt(mean((estimate32-trueState).^2))+MSEsGB2;
 MSEsPF =sqrt(mean((estimate-trueState).^2))+MSEsPF;
 MSEsPF2=sqrt(mean((estimate2-trueState).^2))+MSEsPF2;
end
%% calculate RMSE
 MSEsEKF =MSEsEKF/Test
 MSEsEKF2=MSEsEKF2/Test
 MSEsGB =MSEsGB/Test
 MSEsGB2=MSEsGB2/Test
 MSEsPF =MSEsPF/Test
 MSEsPF2=MSEsPF2/Test

%% plot comparison of different filters

 figure()
 boxplot([abs(xEstimate(1,:)-trueState(1,:))',abs(xEstimate2(1,:)-trueState(1,:))'],'Notch','on','Labels',{'extended kalman','iterated extended kalman'})
 title('Absolute error of different kalman filters','FontSize',14)
 xlabel('Different Strategies','FontSize',14)
 ylabel('Absolute Error','FontSize',14)
 txt1 = ['RMSE: ',num2str( MSEsEKF)];
 text(0.7,50,txt1,'FontSize',14)
 txt1 = ['RMSE: ',num2str( MSEsEKF2)];
 text(1.7,50,txt1,'FontSize',14)
 %}

 figure()
 boxplot([abs(estimate32(1,:)-trueState(1,:))',abs(estimate3(1,:)-trueState(1,:))'],'Notch','on','Labels',{'grid-based max estimate','grid-based average estimate'})
 title('Absolute error of different grid-based filters','FontSize',14)
 xlabel('Different Strategies','FontSize',14)
 ylabel('Absolute Error','FontSize',14)
 txt1 = ['RMSE: ',num2str( MSEsGB)];
 text(1.7,20,txt1,'FontSize',14)
 txt1 = ['RMSE: ',num2str( MSEsGB2)];
 text(0.7,20,txt1,'FontSize',14)
 %}

 figure()
 boxplot([abs(estimate(1,:)-trueState(1,:))',abs(estimate2(1,:)-trueState(1,:))'],'Notch','on','Labels',{'Particle Filter','SIR Particle Filter'})
 title('Absolute error of different particle filters','FontSize',14)
 xlabel('Different Strategies','FontSize',14)
 ylabel('Absolute Error','FontSize',14)
 txt1 = ['RMSE: ',num2str( MSEsPF)];
 text(1.7,10,txt1,'FontSize',14)
 txt1 = ['RMSE: ',num2str( MSEsPF2)];
 text(0.7,10,txt1,'FontSize',14)
 %}