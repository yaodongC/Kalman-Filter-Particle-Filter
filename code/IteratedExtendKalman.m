function [state,errorCov]=IteratedExtendKalman(prevState,prevCov,measure,measNoiseCov,proNoiseCov,kn)
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
 % use matlab jaccobian matrix
%{
syms variable1 
AJacc=jacobian([variable1/2+25*variable1/(1+variable1 ^2)+8*cos(1.2*kn)], [variable1]);
A=vpa(subs(AJacc,{variable1},{prevState}));
syms variable2 
HJacc=jacobian([variable2^2/20], [variable2]);
H=vpa(subs(HJacc,{variable2},{xPrej}));
%}
 % use third party
 A=NumJacob(@f_Func,x,kn);
% step 2.2 calculate error covariance Pcov=A*Pcov*A+Q;
 Pcov=A*Pcov*A+Q;
 xTemp= xPrej;
 pTemp= Pcov;
 
  for i=1:100
   H=NumJacob(@g_Func, xTemp);
 %}
 % use if input are matrix
 %{
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

%% update correction part
% step 3 compute Kalman filter gain
KG=Pcov*H*(H*Pcov*H+R)^(-1);
% errorCov
pTemp=(1-KG*H)*Pcov;
% step 4 Update the new estimate
xTemp=xPrej+KG*(z-g_Func(xTemp)-H*(xPrej-xTemp));          
% step 5 Update the error covariance
 end
 errorCov= pTemp;
 state= xTemp;
end