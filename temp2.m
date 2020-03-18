pkg load optim
clear ;
%%close all;
clc;

ird = dlmread("data/dpc-covid19-ita-andamento-nazionale.csv", ',');

rplusd = ird(:,3)
rplusd_delta=[0; rplusd(2:size(rplusd)) - rplusd(1:size(rplusd)-1)];
infected= ird(:,1);

m = length(infected);
% Add a column of all ones (intercept term) to x
X = [ones(m, 1) infected];
% Calculate theta
theta = (pinv(X'*X))*X'*rplusd

plot(infected,rplusd,'+b');