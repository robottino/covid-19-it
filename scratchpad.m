
N=@(x,p) p(1) * normpdf(x,p(2),p(3));
x=[1:length(didt)]';
mcamp = (x' * didt) ./ sum(didt);
gain = max(didt) * dscamp * sqrt(2*pi);
dscamp = sqrt(((x-mcamp)' * (x-mcamp)) ./ length(x));
init=[gain,mcamp,dscamp];
[f,p,c,i] = leasqr (x,didt,init,N);

plot(x,didt,x,N(x,p));

%%---------------------------------

F=@(x,p) p(4) + p(3)*exp(p(1)*(x-p(2)))
init=[1,0,100,0];
[f,p,c,i] = leasqr (x_it,y_it,init,F);

plot(x_it,y_it,'o');

%%---------------------------------

drrdt=ddt(recovered,0);
dddt=ddt(deaths,0);

function f = gsqm(v,w)
  l=length(v);
  res=[];
  for i= 0:l-1
    vs = v(1:end-i);
    ws = w(i+1:end);
    res = [res,mean((ws-vs).^2)];
  endfor
  f = res;
endfunction
mminterval=5;
%%y=mm(drrdt./recovered,mminterval);
x=mm(didt,mminterval);
%%y=mm(drrdt,mminterval);
y=mm(dddt,mminterval);
bar (gsqm(x,y))


sqm = @(x,y,i) mean((x(1:end-i)-y(i+1:end)).^2)
limit=length(didt);
x=[1:limit];
vcorr = []; for i = 0:limit-2
  vcorr = [vcorr,sqm(dddt,didt,i)];
endfor

bar(vcorr);
%%plot(x,vcorr,'o');


drrdt=ddt(recovered,0);
bar([drdt,dddt])

