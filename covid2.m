pkg load optim
clear ;
%%close all;
clc;

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
  N = x(1)+x(2)+x(3);
  beta = 0.343/N;
  %%gamma = 5000/N;
  %%gamma = 0.22683 / 2.5;
  %%gamma = 0.094262;
  %%gamma = 0.086305;
  gamma = 0.11734;
  
  %% S
  xdot(1) = - beta* x(2) * x(1);

  %% I
  xdot(2) = (beta * x(2) * x(1)) - gamma * x(2);

  %% R
  xdot(3) = gamma * x(2);
endfunction

function y = delta(x,x0)
  [x0;x(2:length(x)) - x(1:length(x)-1)]
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

%% estimate gamma
recovered=data(:,7);
infected=data(:,5);
gamma=(infected' * recovered) * pinv(infected' * infected)



t = linspace (0, 200, 200)';
x = lsode ("sir", [66318.20533; 0.22683; 0], t); 
I = @(x,p) p(1) ./ (1+exp(-p(2)*(x-p(3)))); init_I=[0,0,0];

data_it = load('data/covid-19-data-it.txt');
x_it = [1:length(data)]';
y_it = data(:,9);
%%p_it = estimate(I,x_it,y_it,init_I);
[f_it, p_it, cvg_it, iter_it] = leasqr (x_it, y_it, init_I, I);


plot(
  x_it,y_it,'linewidth',2,'o'
  ,t, x(:,2) + x(:,3)
  ,t,x
  );
  legend ({"Total cases Italy","x","Susceptible","Infected","Recovered"}, "location", "east");
 set (gca, "xgrid", "on");
 set (gca, "ygrid", "on");
