pkg load optim
clear ;
%%close all;
clc;
more off;

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

data_orig = dlmread("data/dpc-covid19-ita-andamento-nazionale.csv", ',');
#{
1:  ricoverati_con_sintomi
2:  terapia_intensiva
3:  totale_ospedalizzati
4:  isolamento_domiciliare
5:  totale_attualmente_positivi
6:  nuovi_attualmente_positivi
7:  dimessi_guariti
8:  deceduti
9:  totale_casi
10:  tamponi
#}

%%data(:,9) = data(:,5)+data(:,7)+data(:,8);

days_back=0;
%%data = data_orig(1:size(data_orig)(1)-days_back,1:size(data_orig)(2));
data = data_orig(1:end-days_back,1:end);

%% guestimate beta (not rigorous)
x_it = [1:length(data)]';
y_it = data(:,9);
%%I = @(x,p) p(1) ./ (1+exp(-p(2)*(x-p(3)))); init_I=[0,0,0];
I = @(x,p) p(4) + p(1) ./ (1+exp(-p(2)*(x-p(3)))); init_I=[0,0.1,10,0];
[f_it, p_it, cvg_it, iter_it] = leasqr (x_it, y_it, init_I, I);
beta = p_it(2);
S0=p_it(1);

%% estimate gamma = (dR/dt) / I
%% in SIR model, deaths are considered in the "R" state.
recovered=data(:,7);
deaths=data(:,8);
removed=recovered+deaths;
%%drdt=[0;removed(2:length(removed)) - removed(1:length(removed)-1)];
drdt = ddt(removed,0);
infected=data(:,5);
gamma=theta(infected,drdt);

%% estimate bets*S = (dI/dt + dR/dt) / I
%%didt=data(:,6);
didt=ddt(infected,0);
didrdt=didt+drdt;
betaS=theta(infected,didrdt);

%% estimate dS/dt = -gamma*I -dI/dt 
dsdt=-gamma * infected -didt;

%% estimate intense therapy / infected ratio
tins=data(:,2);
ti_ratio=mean(tins ./ infected);

%% http://www.salute.gov.it/imgs/C_17_pubblicazioni_1203_ulterioriallegati_ulterioreallegato_10_alleg.pdf
intense_care_spots = 7981;

start_date = datenum (2020, 2, 24);
num_days = 100;

t = linspace (0, 120, 1000)'; t=t+start_date;
x = lsode ("sir", [S0; 450; 0], t);

current_timestamp=datenum(datevec(date()));

II=[t,[-1000;diff(x(:,2))]];
[m1,row_peak] = min(min(abs(II),[],2));
[m2,row_today]=min(min([abs(t-current_timestamp)],[],2));

plot(
%%semilogy(
  t,x,'linewidth',2
  ,t,x(:,2) * ti_ratio,'linewidth',2
  ,t,ones(length(t),1)*intense_care_spots,'--','linewidth',1
  ,t(row_today),x(row_today:row_today,2),'o','linewidth',2
  
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
