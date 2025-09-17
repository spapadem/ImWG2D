function g = green_inf(x1,y1,x2,y2,kn,depth,nmod,R)
% y1 = column array
% y2 = row array
m = size(x1,1);
n = size(x1,2);
g = zeros(m,n);
g1 = zeros(m,n);
k = 1:nmod;
lambda = (k*pi) / (depth);
mu = sqrt(kn*kn - lambda.*lambda);

for n = k

    Xn = @(x) sqrt(2/depth)*sin(lambda(n)*x);
   g = g + 0.5*1i/mu(n)*(exp(1i*mu(n)*abs(x1-x2))...
      ).*Xn(y1).*Xn(y2);

end
