pkg load optim

  function params = estimate (F, x, y, init)
    global verbose; verbose = false;
    [f, p, cvg, iter] = ...
      leasqr (x, y, init, F);
    params=p;      
  endfunction

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
  beta = 0.00001582;
  gamma = 0.00000001;
  
  %% S
  xdot(1) = - (beta * x(2) * x(1)) / N;

  %% I
  xdot(2) = (beta * x(2) * x(1)) / N - gamma * x(2);

  %% R
  xdot(3) = gamma * x(2);
endfunction

function xdot = sir2(x,t)
  beta = 6.8958e+04;
  gamma = 0.5;
  
  %% S
  xdot(1) = - beta * x(2) * x(1);

  %% I
  xdot(2) = (beta * x(2) * x(1)) - gamma * x(2);

  %% R
  xdot(3) = gamma * x(2);
endfunction

%%x = lsode ("f", [1; 2], (t = linspace (0, 50, 50)')); plot(t,x);
 x = lsode ("sir2", [15000; 3; 0], (t = linspace (0, 100, 200)')); plot(t,x,'linewidth',2);

%%  logistica
I = @(x,p) p(1)./(exp(p(2) * x + p(1)) - 1) + p(1)
 
data_china = load('data/covid-19-data-china.txt');
x_china = data_china(:,1);
y_china = data_china(:,2);
p_china = estimate(I,x_china,y_china,[0.1,0.1]);

data_it = load('data/covid-19-data-it.txt');
x_it = data_it(:,1);
y_it = data_it(:,2);
p_it = estimate(I,x_it,y_it,[1,1]);

disp(p_china);
disp(p_it);
