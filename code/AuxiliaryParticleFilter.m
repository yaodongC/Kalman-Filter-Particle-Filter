function [state2,particleWeight,estimate]=AuxiliaryParticleFilter(prevState2,preWeight2,measure2,processNoiseCov2,measureNoiseCov2,kn)
%% initial variables
x2=prevState2;
z2=measure2;
Q2=processNoiseCov2;
R2=measureNoiseCov2;
lengthN=length(prevState2);
%pindex=preIndex;
%k=kn
%% partical filter begain:
% previous weight are all set equal
wPre=preWeight2;
%% prediction part:
% step 2 drew a sample from particle index
mu=randi([1 lengthN],1,lengthN);
%mu=pindex(mu);
wPre2=wPre(mu);
tempState=x2(mu);
%tempState=x2;
%xPrejTemp=f_Func(tempState,kn)+sqrt(Q2)*randn(size(tempState));
xPrejTemp=f_Func(tempState,kn)+sqrt(Q2)*randn(size(tempState));
zPrejTemp=g_Func(xPrejTemp);
particleWeightTemp = wPre2.*normpdf((z2-zPrejTemp),0,sqrt(R2));
%particleWeightTemp = normpdf((z2-zPrejTemp),0,sqrt(R2));
%particleWeightTemp = wPre.*normpdf((z2-zPrejTemp),0,sqrt(R2));
%ji=round(particleWeightTemp*lengthN*lengthN);
%ji=randi([1 lengthN],1,lengthN);
% use optimal importance density
particleWeightTemp = particleWeightTemp./sum(particleWeightTemp);
%% Resampling part:
% step 3 resample from the cummulative distribution of the probability distribution
indexT = randsample(1:lengthN, lengthN, true, particleWeightTemp);
x3=x2(indexT);
%x3=x2(ji);
%x3=wPre.*zPrejTemp.*xPrejTemp;
% step 1 project the state and measurements
xPrej2=f_Func(x3,kn)+sqrt(Q2)*randn(size(x3));
%xPrej2=f_Func(x2,kn)+sqrt(Q2)*randn(size(x2));
zPrej2=g_Func(xPrej2);
% step 2 Generate the weights for each of these particles.
% use prior as importance density
 particleWeight = normpdf((z2-zPrej2),0,sqrt(R2))./particleWeightTemp;
% use optimal importance density
 particleWeight = particleWeight./sum(particleWeight);
%% Resampling part:
% step 3 resample from the cummulative distribution of the probability distribution
Neff=1/sum(particleWeight.^2);
if Neff < 0.75*lengthN
index = randsample(1:lengthN, lengthN, true, particleWeight);
% extract new particles
state2=xPrej2(:,index); 
else
state2=xPrej2;
end
%% update the state estimate
estimate=mean(state2);
%% update the state estimate
estimate=mean(state2);
end