pkg load optim
clear ;
%%close all;
clc;

recovered = dlmread("data/covid-19-data-it-recovered.csv", ',');
recovered_delta=[[1:size(recovered)]',[0; recovered(2:size(recovered),2) - recovered(1:size(recovered)-1,2)]];

deaths = dlmread("data/covid-19-data-it-deaths.csv", ',');
deaths_delta=[[1:size(deaths)]',[0; deaths(2:size(deaths),2) - deaths(1:size(deaths)-1,2)]];

start=1;
x_drdt=recovered_delta(start:size(recovered_delta,1),1);
%y_drdt=recovered_delta(start:size(recovered_delta,1),2);
y_drdt=recovered_delta(start:size(recovered_delta,1),2)+deaths_delta(start:size(recovered_delta,1),2);

%%DR = @(x,p) p(1) * x ;
%%DR = @(x,p) p(1) * exp(p(2) * x) ;
DR = @(x,p) p(1) ./ (1+exp(-p(2)*(x-p(3)))); init_DR=[0,0,0];
[f_drdt, p_drdt, cvg, iter]=leasqr(x_drdt,y_drdt,init_DR,DR);

plot(x_drdt,y_drdt,'o',x_drdt,DR(x_drdt,p_drdt));
disp(p_drdt);


