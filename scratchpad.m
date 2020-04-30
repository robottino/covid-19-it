
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
%    vs = v(1:end-i);
%    ws = w(i+1:end);
    vs = [v(1:end-i);zeros(i,1)];
    ws = [zeros(i,1);w(i+1:end)];
    res = [res,mean((ws-vs).^2)];
  endfor
  f = res;
endfunction
mminterval=1;
%%y=mm(drrdt./recovered,mminterval);
x=mm(didt,mminterval);
%%y=mm(drrdt,mminterval);
y=mm(dddt,mminterval);
%bar (gsqm(x,y))
bar(xcorr(x,y))


%%---------------------------------

pkg load signal;
drrdt=ddt(recovered,0);
dddt=ddt(deaths,0);

mminterval=5;
x=mm(drdt,mminterval);
y=mm(didt,mminterval);

[r,l]= xcorr(x,y,'unbiased');
plot(l,r);

