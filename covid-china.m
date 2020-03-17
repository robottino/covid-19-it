pkg load optim
  F = @(x,p) p(1) + p(2) * exp(p(3) * (x - p(4)))./(1 + exp(p(3) * (x - p(4))));
  dF = @(x,p) (p(2) * p(3) * exp(p(3) * (p(4) + x)))./(exp(p(3) * p(4)) + exp(p(3) * x)).^2
  
  function params = estimate (F, x, y)

    p1=min(y);
    p2=max(y)-min(y);
    p3=4*max(diff(y)./diff(x))/p2;
    p4=mean(x);

    init=[p1,p2,p3,p4];
    
    %% other configuration (default values):
    tolerance = .00000001;
    max_iterations = 200;
    weights = ones (size (y));
    dp = [.00001; .00001; .00001; .00001]; % bidirectional numeric gradient stepsize
    dFdp = "dfdp"; % function for gradient (numerical)
    %% dFdp = dF;
  
    minstep = [0.001; 0.001;0.001; 0.001];
    maxstep = [0.08; 0.08; 0.08; 0.08];
    options = [minstep, maxstep];
    
    global verbose; verbose = false;
    [f, p, cvg, iter] = ...
      leasqr (x, y, init, F, tolerance, max_iterations,weights,dp,dFdp);
  
    params=p;      
  endfunction

  data_china = load('data/covid-19-data-china.txt');
  x_china = data_china(:,1);
  y_china = data_china(:,2);
 
  p_china = estimate(F,x_china,y_china);

%%  data_it = load('data/covid-19-data-it.txt');
  data_it = load('data/ita2test.csv');

  x_it = data_it(:,1);
  y_it = data_it(:,2);
 
  p_it = estimate(F,x_it,y_it);
  p_it_1d = estimate(F,x_it(1:size(x_it)(1)-1,1),y_it(1:size(y_it)(1)-1,1));
  p_it_2d = estimate(F,x_it(1:size(x_it)(1)-2,1),y_it(1:size(y_it)(1)-2,1));
  p_it_3d = estimate(F,x_it(1:size(x_it)(1)-3,1),y_it(1:size(y_it)(1)-3,1));
  
  %%x=[0:100]';
  x=linspace(1,100,200);
  
  plot(x,F(x,p_china),'linewidth',1,'--m'
    ,x_china,y_china,'.b','linewidth',1
    ,x_it,y_it,'or','linewidth',2
    ,x,F(x,p_it),'linewidth',2
%%    ,x,F(x,p_it_1d),'linewidth',2
%%    ,x,F(x,p_it_2d),'linewidth',2
    );

 xlabel('Days since January 22, 2020');
 ylabel('Total cases');
 set (gca, "xgrid", "on");
 set (gca, "ygrid", "on");


