function result=Ornstein_Uhlenbeck(tau,sd,t)
%Input definition:
%tau: Noise correlation time. Scalar.
%sd: Noise standard deviation. Scalar.
%t: Time array. 1D array.

%Output definition:
%result: 1D array of noise figures in the Ornstein-Uhlenbeck process
%corresponding to time in the input t.

result=nan(size(t));
result(1)=normrnd(0,sd);
rndarray=normrnd(0,sqrt(t(2:end)-t(1:end-1)));
for h=2:numel(result)
    result(h)=result(h-1)-(1/tau)*result(h-1)*(t(h)-t(h-1))+sqrt(2*(1/tau)*(sd^2))*rndarray(h-1);
end