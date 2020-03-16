pkg load optim

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

%%x = lsode ("f", [1; 2], (t = linspace (0, 50, 50)')); plot(t,x);
 x = lsode ("sir", [100; 3; 0], (t = linspace (0, 1000, 200)')); plot(t,x,'linewidth',2);

%%  logistica
 %% x(t) = (a e^(a c_1 + b t))/(e^(a c_1 + b t) - 1)

