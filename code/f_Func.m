function state=f_Func(preState,k)

state=0.5*preState+(25*preState)./(1+preState.^2)+8*cos(1.2*(k-1)); 
