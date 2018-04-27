function [state2,particleWeight,estimate]=SIRParticleFilter(prevState2,measure2,processNoiseCov2,measureNoiseCov2,kn)
%% initial variables
x2=prevState2;
z2=measure2;
Q2=processNoiseCov2;
R2=measureNoiseCov2;
lengthN=length(prevState2);
%% partical filter begain:
% previous weight are all set equal
wPre=ones(1,lengthN)/lengthN;
%% prediction part:
% step 1 project the state and measurements
xPrej2=f_Func(x2,kn)+sqrt(Q2)*randn(size(x2));
zPrej2=g_Func(xPrej2);
% step 2 Generate the weights for each of these particles.
% use prior as importance density
 particleWeight = wPre.*normpdf((z2-zPrej2),0,sqrt(R2));
% normalize weights
 particleWeight = particleWeight./sum(particleWeight);
%% Resampling part:
% step 3 resample every time
index = randsample(1:lengthN, lengthN, true, particleWeight);
state2=xPrej2(:,index); 
%% update the state estimate
estimate=mean(state2);
end