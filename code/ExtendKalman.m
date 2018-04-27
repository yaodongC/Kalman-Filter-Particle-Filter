function [state,errorCov,kl]=ExtendKalman(prevState,prevCov,measure,measNoiseCov,proNoiseCov,kn)
%% initial variables
x=prevState;
Pcov=prevCov;
Q=proNoiseCov;
R=measNoiseCov;
z=measure;
%% kalman begain:
%% prediction part:
% step 1 project the state
xPrej=f_Func(x,kn);
zPrej=g_Func(xPrej);
% step 2 project the error covariance
% step 2.1 claculate jaccobian matrix A and H
% use matlab jaccobian matrix (it is slow, but stable)
%{ 
syms variable1 
AJacc=jacobian([variable1/2+25*variable1/(1+variable1 ^2)+8*cos(1.2*kn)], [variable1]);
A=vpa(subs(AJacc,{variable1},{prevState}));
syms variable2 
HJacc=jacobian([variable2^2/20], [variable2]);
H=vpa(subs(HJacc,{variable2},{xPrej}));
%}
 % use third party jaccobian, it is faster than matlab's.
 % cite : Youngmok Yun, 2013.05.04
 A=NumJacob(@f_Func,x,kn);
 H=NumJacob(@g_Func,xPrej);
 %}
 %% update correction part
 %% use if input are matrix
% step 2.2 calculate error covariance
Pcov=A*Pcov*A'+Q;
%% update correction part
% step 3 compute Kalman filter gain
KG=Pcov*H'*inv(H*Pcov*H'+R);
% step 4 Update the new estimate
state=xPrej+KG*(z-zPrej);          
% step 5 Update the error covariance
errorCov=Pcov-KG*Pcov'*H; 
%}
%% if not use this (this assignment only has one state) 
%{
% step 2.2 calculate error covariance
Pcov=A*Pcov*A+Q;
% step 3 compute Kalman filter gain
KG=Pcov*H*(H*Pcov*H+R)^(-1);
% step 4 Update the new estimate
state=xPrej+KG*(z-zPrej);          
% step 5 Update the error covariance
%errorCov=(1-KG*H)*Pcov*(1-KG*H)+KG*R*KG; % error covariance on the PPT
errorCov=(1-KG*H)*Pcov; 
%}
%kl=Pcov*log(Pcov);
end