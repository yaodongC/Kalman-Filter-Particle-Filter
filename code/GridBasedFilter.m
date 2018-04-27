function [state2,newGridWeight,estimate]=GridBasedFilter(prevState2,prevWeight2,measure2,processNoiseCov2,measureNoiseCov2,kn,Xlength,Ylength)
%% initial variables
x2=prevState2;
Q2=processNoiseCov2;
R2=measureNoiseCov2;
z2=measure2;
wPre=prevWeight2;
lengthN=length(prevState2);
gridX=Xlength;
gridY=Ylength;
%% grid-base filter begain:
%% prediction grid:
% go through grid
% step 1 project the state 
Prej2X=f_Func(prevState2,kn)+sqrt(Q2)*randn(size(x2));
for i=1:1:gridY
 % step 2 Generate the weights for each of these states.   
   Prob=normpdf((Prej2X(i)-Prej2X),0,sqrt(Q2));
  % sum(Prob)
   gridWeight(i)=sum(wPre(i)*Prob);
end
% normalize
gridWeight = gridWeight./sum(gridWeight);
zPrej2=g_Func(Prej2X);
%% update grid
postProb=normpdf((z2-zPrej2),0,sqrt(R2));
% new weight
newGridWeight= gridWeight.*postProb;
% normalize
newGridWeight = newGridWeight./sum(newGridWeight);
% calculate estimate
estimate = sum(Prej2X.*newGridWeight);
%estimate = Prej2X(find(max(newGridWeight)));
state2=Prej2X;
end