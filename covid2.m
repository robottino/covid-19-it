pkg load optim
clear ;
%%close all;
clc;
more off;

%data fase 2: 4/5/2020
% gilet arancioni: 30/5/2020

global beta gamma start_date;

function xdot = f (x,t)
  r = 0.25;
  k = 1.4;
  a = 1.5;
  b = 0.16;
  c = 0.9;
  d = 0.8;
  xdot(1) = r * x(1) * (1 - x(1)/k) - a*x(1)*x(2)/(1 + b*x(1));
  xdot(2) = c*a*x(1)*x(2)/(1 + b*x(1)) - d*x(2);     
endfunction

function xdot = sir(x,t)
  global beta gamma start_date;
  N = x(1)+x(2)+x(3);
  b_beta = 1 * beta/N;
  %%if (t - start_date>30) b_beta = b_beta * 0.7; endif
  
  %% S
  xdot(1) = - b_beta* x(2) * x(1);

  %% I
  xdot(2) = (b_beta * x(2) * x(1)) - gamma * x(2);

  %% R
  xdot(3) = gamma * x(2);
endfunction

function y = ddt(x,x0)
  y=[x0;x(2:length(x)) - x(1:length(x)-1)];
endfunction

function y = theta(X,yy)
  y= (pinv(X'*X))*X'*yy;
endfunction

function y = mm(x,w)
  y=filter(ones(w,1)/w, 1, x);
endfunction

data_orig = dlmread("data/dpc-covid19-ita-andamento-nazionale.csv", ',');
#{

1:ricoverati_con_sintomi
2:terapia_intensiva
3:totale_ospedalizzati
4:isolamento_domiciliare
5:totale_positivi
6:variazione_totale_positivi
7:nuovi_positivi
8:dimessi_guariti
9:deceduti
10:casi_da_sospetto_diagnostico
11:casi_da_screening
12:totale_casi
13:tamponi
14:casi_testati

#}

days_back=0;
%%data = data_orig(1:size(data_orig)(1)-days_back,1:size(data_orig)(2));
data = data_orig(1:end-days_back,1:end);

tamponi=data(:,13);
dtdt=ddt(tamponi,0);

%% guestimate beta (not rigorous)
x_it = [1:length(data)]';
y_it = data(:,12);
%%I = @(x,p) p(1) ./ (1+exp(-p(2)*(x-p(3)))); init_I=[0,0,0];
I = @(x,p) p(4) + p(1) ./ (1+exp(-p(2)*(x-p(3)))); init_I=[0,0.1,10,0];
[f_it, p_it, cvg_it, iter_it] = leasqr (x_it, y_it, init_I, I);
beta = p_it(2);
S0=p_it(1);
beta=1/9;
beta=1.3848e-01;
S0=2.*1.2010e+05;
%S0=1.4968e+05*2;

%% estimate gamma = (dR/dt) / I
%% in SIR model, deaths are considered in the "R" state.
recovered=data(:,8);
deaths=data(:,9);
removed=recovered+deaths;
%%drdt=[0;removed(2:length(removed)) - removed(1:length(removed)-1)];
drdt = ddt(removed,0);
infected=data(:,5);
gamma=theta(infected,drdt);
gamma=1/3.1824e+01;
%gamma=3.5638e-02;

didt=data(:,7);

%% estimate intense therapy / infected ratio
tins=data(:,2);
ti_ratio=mean(tins ./ infected);

%%https://www.truenumbers.it/coronavirus-terapia-intensiva/
intense_care_spots = 6864;

%% lockdown date: 2020, 3, 9
start_date = datenum (2020, 2, 24);
num_days = 100;

t = linspace (0, 220, 1000)'; t=t+start_date;
x = lsode ("sir", [S0; 4000; 0], t);

current_timestamp=datenum(datevec(date()));

II=[t,[-1000;diff(x(:,2))]];
[m1,row_peak] = min(min(abs(II),[],2));
[m2,row_today]=min(min([abs(t-current_timestamp)],[],2));

figure(1);
plot(
%%semilogy(
  t,x,'linewidth',2
  ,t,x(:,2) * ti_ratio,'linewidth',2
  ,t,ones(length(t),1)*intense_care_spots,'--','linewidth',1
  ,t(row_today),x(row_today:row_today,2),'o','linewidth',2
  
  ,x_it+start_date,infected,'.'
  ,x_it+start_date,recovered,'.'
  %%,x_it,y_it,'o'
  %%,t,I(t,p_it)
  %%,t, x(:,2) + x(:,3)
  );
  legend ({"Susceptible","Infected","Removed (recovered+deaths)","Intensive care","Intensive care max capacity","Today"}, "location", "northeast");
 set (gca, "xgrid", "on");
 set (gca, "ygrid", "on");
%% strtitle=["Peak: " str(10) "infected, @ " datestr(t(row),'dd-mmm-yyyy')];
 %%xlabel();
 title(sprintf("Peak: %d infected @ %s",round(x(row_peak:row_peak,2)),datestr(t(row_peak),'dd-mmm-yyyy')));
 datetick ("x", "dd mmm");

figure(2);
  didtperc=didt./dtdt;
  start_date = datenum (2020, 2, 24);
  t=[1:length(didtperc)];
  v3=mm(didtperc,3);
  v7=mm(didtperc,7);
  v14=mm(didtperc,14);
  plot(t+start_date,v3,'.',t+start_date,v7,'.',t+start_date,v14,'.');
  set (gca, "xgrid", "on");
  set (gca, "ygrid", "on");
  set (gca, "xminorgrid", "on");
  set (gca, "yminorgrid", "on");
  datetick ("x", "dd mmm","keeplimits");

  
figure(3);
plot(
%%semilogy(
  %%t,ddt(x(:,2),0)
  %x_it,filter(ones(3,1)/3, 1, didt)
  x_it+start_date, mm(didt,7)
  );
  legend ({"Susceptible","Infected","Removed (recovered+deaths)","Intensive care","Intensive care max capacity","Today"}, "location", "northeast");
 set (gca, "xgrid", "on");
 set (gca, "ygrid", "on");
 datetick ("x", "dd mmm");
 
figure(4);
 
 format short g;
 status = [didt,-ddt(recovered,0),-ddt(deaths,0),-drdt,didt-drdt,dtdt, didt./dtdt.*100]
 plot([1:length(didt)],status);
 legend ({"new infected","new recovered","new deaths","new removed","delta","new tests","% tests"}, "location", "northeast");

%------------ Figura stima infetti

INF = @(x,p) (p(1) ./ (1+exp(-p(2).*(x-p(3))))) - (p(1) ./ (1+exp(-p(4).*(x-p(5)))));
dINF = @(x, f, p, dp, func) p(1) .* ((p(2) .* exp(p(2) .* (p(3) + x)))./(exp(p(2) .* p(3)) + exp(p(2) .* x)).^2 - (p(4) .* exp(p(4) .* (p(5) + x)))./(exp(p(4) .* p(5)) + exp(p(4) .* x)).^2);

init_inf= [180000,1/12,35,1/20,83];

x_inf=[1:length(infected)]';
y_inf=infected;

weights=ones (size(y_inf));
%weights=[1:length(y_inf)]';

minstep = 0.00000001*ones(length(init_inf),1);
maxstep = 0.0001*ones(length(init_inf),1);
options = [minstep, maxstep];

#{

%[f,p,c,i]=leasqr (x_inf, y_inf, init_inf, INF, 0.00000001,10000, ones (size (y_inf)) ,0.001 * ones (size (init_inf)), 'dfdp',options);
[f,p,c,i]=leasqr (x_inf, y_inf, init_inf, INF, 0.00000001,10000, weights ,0.001 * ones (size (init_inf)), 'dfdp' ,options);

start_date = datenum (2020, 2, 24);
t_inf=linspace(0,250,200);
figure(5);
plot(x_inf+start_date,y_inf,'.',t_inf+start_date,INF(t_inf,p));
set (gca, "xgrid", "on");
set (gca, "ygrid", "on");
set (gca, "xminorgrid", "on");
set (gca, "yminorgrid", "on");
datetick ("x", "dd mmm","keeplimits");

#} 
 
 
 
 
 