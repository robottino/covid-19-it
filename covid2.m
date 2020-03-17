pkg load optim

  function params = estimate (F, x, y, init)
    [f, p, cvg, iter] = leasqr (x, y, init, F);
    params = p;
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
  beta = 0.22683/N;
  gamma = 1000/N;
  
  %% S
  xdot(1) = - beta* x(2) * x(1);

  %% I
  xdot(2) = (beta * x(2) * x(1)) - gamma * x(2);

  %% R
  xdot(3) = gamma * x(2);
endfunction

%%x = lsode ("f", [1; 2], (t = linspace (0, 50, 50)')); plot(t,x);
t = linspace (0, 100, 200)';
 x = lsode ("sir", [66318.20533; 0.22683; 0], t); 

%%  logistica
%%I = @(x,p) p(1)./(exp(p(2) * x + p(1)) - 1) + p(1)

%% L / (1+ exp(-k(x-x0))); L=p(1), k=p(2), x0=p(3)
I = @(x,p) p(1) ./ (1+exp(-p(2)*(x-p(3)))); init_I=[0,0,0];

#{ 
data_china = load('data/covid-19-data-china.txt');
x_china = data_china(:,1);
y_china = data_china(:,2);
p_china = estimate(I,x_china,y_china,init);
#}

data_it = load('data/covid-19-data-it.txt');
x_it = data_it(:,1);
y_it = data_it(:,2);
%%p_it = estimate(I,x_it,y_it,init_I);
[f_it, p_it, cvg_it, iter_it] = leasqr (x_it, y_it, init_I, I);

%%disp(p_china);
disp(p_it);

plot(x_it,y_it,'linewidth',2,'o',
  %%t,I(t,p_it),'linewidth',2,
  t,x(:,2) + x(:,3),
  t,x);
 set (gca, "xgrid", "on");
 set (gca, "ygrid", "on");
