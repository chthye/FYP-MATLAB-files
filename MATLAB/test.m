p1 = 0.5;
A = 1;
mu = 0;
sigma = 0.1;
digits(40);
N = makedist('Normal','mu',mu,'sigma',sigma);
thresh = linspace(-A,A,2001);
%figure(1);
%plot(thresh, pdf(N,thresh))
%figure(2);

e = (0.5)*(cdf(N,thresh+A,'upper'))+0.5*(cdf(N,thresh-A));
%note: upper returns complement of cdf, much more accurate than 1-cdf
e = -vpa(log10(e));
plot(thresh,e)
hold on

thresh2 = linspace(-1,1,101);
numBits = 10000;
e2 = zeros(1,101);

for i = 0:100
    x = randi([0 1],1,numBits);
    conditionX = (x==0);
    x(conditionX) = 1;
    x(~conditionX) = -1;
    n = normrnd(mu,sigma,1,numBits);%mu,sigma,numrows,numcolumns
    y = x+n;
    threshold = -1+0.02*i;
    y2 = y;
    conditionY2 = (y2 >= threshold);
    y2(conditionY2) = 1;
    y2(~conditionY2) = -1;
    z = x-y2;
    z = abs(z);
    e2(i+1) = 0.5*sum(z)/numBits;
end
e2 = -log10(e2)
plot(thresh2,e2)
grid
hold off    
