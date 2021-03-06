pkg load optim
%% Example for linear inequality constraints.
  %% model function:
  %F = @ (x, p) (1 + exp(-p(2) * (x - p(3)))).^(-1)*p(1);
  
  F = @(x,p) p(1) + p(2) * exp(p(3) * (x - p(4)))./(1 + exp(p(3) * (x - p(4))));
 
  data = load('data/covid-19-data-it.txt');
  %%data = load('data/covid19-it.txt');
  x = data(:,1);
  y = data(:,2);

  p1=min(y);
  p2=max(y)-min(y);
  p3=4*max(diff(y)./diff(x))/p2;
  p4=mean(x);

  %% initial values:
  %init = [50000.0; 30.0; 0.2];
  init=[p1,p2,p3,p4];
  %% other configuration (default values):
  tolerance = .00000001;
  max_iterations = 200;
  weights = ones (size (y));
  dp = [.00001; .00001; .00001; .00001]; % bidirectional numeric gradient stepsize
  dFdp = "dfdp"; % function for gradient (numerical)
  
 minstep = [0.001; 0.001;0.001; 0.001];
 maxstep = [0.08; 0.08; 0.08; 0.08];
 options = [minstep, maxstep];

  %% start leasqr, be sure that 'verbose' is not set
  global verbose; verbose = false;
  [f, p, cvg, iter] = ...
      leasqr (x, y, init, F, tolerance, max_iterations,weights,dp,dFdp)
      
      xx=[0:80]';
      
 
 plot(xx,F(xx,p),'linewidth',4,x,y,'or','linewidth',2);
 
 xlabel('Days since January 22, 2020');
 ylabel('Total cases');
 set (gca, "xgrid", "on");
 set (gca, "ygrid", "on");
 
 axis([20 80 0 70000]);
 
