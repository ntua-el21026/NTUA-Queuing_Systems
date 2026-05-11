close all;
clear all;
%load dataset with 1248 elements
data = load('dataset.mat');


y = double(data.values');


% 95% prediction interval
alpha = 0.95;


%number of observations for training set
nob = 1200;
t = [1:nob]';

%apply linear regression with (d+5) parameters, d+1 for polynomial terms
%and 4 for harmonics (sin and cos with T=T_A and T=T_B respectively)
d = 2;
T_A = 1;
T_B = 2;

X=zeros(nob, d+5);

for k = 1:d+1
  X(:,k)=t.^(k-1);
end

X(:,d+2)=cos(2*pi*t/T_A);
X(:,d+3)=sin(2*pi*t/T_A);
X(:,d+4)=cos(2*pi*t/T_B);
X(:,d+5)=sin(2*pi*t/T_B);

%find coefficient values of linear regression
[b,bint,r,rint,stats] = regress(y(1:nob), X, alpha);

t= [1:length(y)]';
%apply the linear regression
yreg = zeros(length(y));
for k = 1:d+1
  yreg = yreg + b(k) * t.^(k-1);
end

yreg = yreg + b(d+2) * cos(2*pi*t/T_A) + b(d+3) * sin(2*pi*t/T_A);
yreg = yreg + b(d+4) * cos(2*pi*t/T_B) + b(d+5) * sin(2*pi*t/T_B);

%plot original data and linear regression at the same figure
figure
plot(t(1:nob), y(1:nob), 'k')
hold
plot(t(1:nob), yreg(1:nob), 'r')
hold off
legend('data', 'regression')
title('Data with linear regression application')
npast=48;

% compute the mean squared error of the predictions obtained from the linear regression model

mse = sum(r.^2)/nob;

display(mse);

%find prediction intervals
%Cf. slide 22, equation 5.4 and slide 17,
cint = norminv((1+alpha)/2)*sqrt(sum(r.^2)/(nob - d - 5));

%apply forecasting
figure
plot(t(nob-npast:nob), y(nob-npast:nob), 'k')

% plot the original test dataset, the point predictions and the prediction intervals at the same figure
hold on
plot(t(nob+1:length(y)), y(nob+1:length(y)), 'o-b')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg)), 'r')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg))+cint, 'r-.')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg))-cint, 'r-.')
title('Forecast data')
hold off


% compute the CRPS metric for each prediction
for k=nob+1:length(y)
    mu(k-nob)=yreg(k);
    sigma(k-nob)=cint/norminv((1+alpha)/2);
end


for k=nob+1:length(y)
 ans=(y(k)-mu(k-nob))/sigma(k-nob);
 crps(k-nob)=sigma(k-nob)*(ans*(2*normcdf(ans)-1)+2*normpdf(ans)-1/sqrt(pi));
end

medianCRPS=median(crps);
display(medianCRPS);
