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
  beta = 0.1;
  gamma = 0;
  
  %% S
  xdot(1) = - (beta * x(2) * x(1)) / N;

  %% I
  xdot(2) = (beta * x(2) * x(1)) / N - gamma * x(2);

  %% R
  xdot(3) = gamma * x(2);
endfunction

function xdot = sir2(x,t)
  beta = 0.1;
  gamma = 0.8;
  
  %% S
  xdot(1) = - beta * x(2) * x(1);

  %% I
  xdot(2) = (beta * x(2) * x(1)) - gamma * x(2);

  %% R
  xdot(3) = gamma * x(2);
endfunction

%%x = lsode ("f", [1; 2], (t = linspace (0, 50, 50)')); plot(t,x);
 x = lsode ("sir2", [100; 3; 0], (t = linspace (0, 3, 200)')); plot(t,x,'linewidth',2);

%%  logistica
 %% x(t) = (a e^(a c_1 + b t))/(e^(a c_1 + b t) - 1)
L = @(x,p) (p(1) * exp(p(1) + p(2) * x))./(exp(p(1) + p(2) * x) - 1);
init_l=[1.0,1.0];
data_china = load('data/covid-19-data-china.txt');
x_china = data_china(:,1);
y_china = data_china(:,2);

F = @(x,p) p(1) + p(2) * exp(p(3) * (x - p(4)))./(1 + exp(p(3) * (x - p(4))));
F3 = @(x,p) (p(1) * 1) ./ (1 + (p(1)-1) * exp(p(2) * x));

p_china = estimate(F,x_china,y_china,[0,0,0,0]);
p_china_2 = estimate(F3,x_china,y_china,[0.1,0.1]);

disp(p_china_2);
