%function y=ConfLims(deg,prob)
%Return lower and upper confidence limits (relative to 1) given the desired
%degrees of freedom and percent confidence.
function y=ConfLims(deg,prob)


p1=(1-(prob/100))/2;
p2=1-p1;

if deg < 20 
	xmax=40;
elseif deg < 40
	xmax=100;
else 
	xmax=deg*2;
end
	
dx=xmax/1000;
x=0:dx:xmax;
m=size(x);
m=m(2);

if deg < 40
	chi=1./(2.^(deg/2)*gamma(deg/2))* x.^(deg/2-1).*exp(-x./2);
	chid=cumsum(chi)*dx;
else
	chid=1/2*(erf((x-deg)/2/sqrt(deg))+1);
end

j=1;
while chid(j) < p1
	j=j+1;
end
cl=x(j)/deg;
while  chid(j) < p2
	j=j+1;
end
ch=x(j)/deg;
y=[cl,ch];
