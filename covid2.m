pkg load optim
clear ;
%%close all;
clc;

global beta gamma;

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
  global beta gamma;
  N = x(1)+x(2)+x(3);
  b_beta = beta/N;
  
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

function y = theta(xx,yy,delay)
  y= ( xx(1:length(xx)-delay)' * yy(1+delay:length(yy)) ) * pinv(xx' * xx);
endfunction

data = dlmread("data/dpc-covid19-ita-andamento-nazionale.csv", ',');
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

%% guestimate beta (not rigorous)
x_it = [1:length(data)]';
y_it = data(:,9);
I = @(x,p) p(1) ./ (1+exp(-p(2)*(x-p(3)))); init_I=[0,0,0];
[f_it, p_it, cvg_it, iter_it] = leasqr (x_it, y_it, init_I, I);
beta = p_it(2);
S0=p_it(1);

%% estimate gamma = (dR/dt) / I
recovered=data(:,7);
%%drdt=[0;recovered(2:length(recovered)) - recovered(1:length(recovered)-1)];
drdt = ddt(recovered,0);
infected=data(:,5);
gamma=(infected' * drdt) * pinv(infected' * infected);

%% estimate bets*S = (dI/dt + dR/dt) / I
didt=data(:,6);
didrdt=didt+drdt;
betaS=(infected' * didrdt) * pinv(infected' * infected);

%% estimate dS/dt = -gamma*I -dI/dt 
dsdt=-gamma * infected -didt;

%% estimate intense therapy / infected ratio
tins=data(:,2);
ti_ratio=mean(tins ./ infected);

%% http://www.salute.gov.it/imgs/C_17_pubblicazioni_1203_ulterioriallegati_ulterioreallegato_10_alleg.pdf
intense_care_spots = 7981;

start_date = datenum (2020, 2, 24);
num_days = 100;

t = linspace (0, 120, 200)'; t=t+start_date;
x = lsode ("sir", [S0; 450; 0], t);

current_timestamp=datenum(datevec(date()));
t2 = linspace (start_date, current_timestamp, 2)';
today_values = lsode ("sir", [S0; 450; 0], t2); 
today_values=today_values(2:2,1:3);

plot(
  t,x,'linewidth',2
  ,t,x(:,2) * ti_ratio,'linewidth',2
  ,t,ones(length(t),1)*intense_care_spots,'--','linewidth',1
  ,current_timestamp,today_values(2),'o','linewidth',2
  
  %%x_it,y_it,'o'
  %%,t,I(t,p_it)
  %%,t, x(:,2) + x(:,3)
  );
  legend ({"Susceptible","Infected","Recovered","Intensive care","Intensive care max capacity","Today"}, "location", "northeast");
 set (gca, "xgrid", "on");
 set (gca, "ygrid", "on");
 %%xlabel('Days since february 24, 2020');
 %%title('ooo');
 datetick ("x", "dd mmm");
