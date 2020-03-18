pkg load optim
clear ;
%%close all;
clc;

data = dlmread("data/dpc-covid19-ita-andamento-nazionale.csv", ',');
#{
1:  data
2:  stato
3:  ricoverati_con_sintomi
4:  terapia_intensiva
5:  totale_ospedalizzati
6:  isolamento_domiciliare
7:  totale_attualmente_positivi
8:  nuovi_attualmente_positivi
9:  dimessi_guariti
10:  deceduti
11:  totale_casi
12:  tamponi
#}

recovered=data(:,9);
infected=data(:,7);

m = length(infected);
% Add a column of all ones (intercept term) to x
X = [ones(m, 1) infected];
% Calculate theta
theta = (pinv(X'*X))*X'*recovered

plot(infected,recovered,'+b');