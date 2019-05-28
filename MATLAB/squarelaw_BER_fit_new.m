p1 = 0.5;
A = 1;
A2 = 1;
B = 0;
DC = 0;
numThresh = 100;
threshMin = B;
threshMax = A; %will be overwritten
threshStep = (threshMax-threshMin)/numThresh;%will be overwritten
mu = 0;
digits(40);

SNR = 19;
sigma = sqrt(0.5*A^2/(2*(10^((SNR)/10))));
N = makedist('Normal','mu',mu,'sigma',sigma);
N1 = makedist('Normal','mu',A,'sigma',sigma);
thresh = linspace(B,A,numThresh+1);

%for Y = (X+N)^2
% E = (0.5)*(cdf(N,sqrt(thresh)-A))-0.5*(cdf(N,-sqrt(thresh)-A))+0.5*cdf(N,-sqrt(thresh)-B)+0.5*(cdf(N,sqrt(thresh)-B,'upper'));
sumcdf = 0;
for x = 0:0.001:0.3
    sumcdf = sumcdf + 0.001*(cdf(N1,real(sqrt(thresh-x)))-cdf(N1,real(sqrt(-thresh-x))))*chi2pdf(x/sigma^2,3)/sigma^2;
end
e = (0.5)*sumcdf+0.5*(chi2cdf(thresh/sigma^2,4,'upper'));

%for Y = (X+N+DC)^2/4
%e = (0.5)*(cdf(N,sqrt(A2*thresh)-A-DC))-0.5*(cdf(N,-sqrt(A2*thresh)-A-DC))+0.5*cdf(N,-sqrt(A2*thresh)-B-DC)+0.5*(cdf(N,sqrt(A2*thresh)-B-DC,'upper'));
%note: upper returns complement of cdf, much more accurate than 1-cdf
e = -log10(e);
% E = -vpa(log10(E));
%find max -log10BER and set to centre
% [M,I] = max(e);
% threshMax = threshMin+2*threshStep*(I-1);
% threshStep = (threshMax-threshMin)/numThresh;
% thresh = linspace(threshMin,threshMax,numThresh+1);


%for Y = (X+N)^2
% e = (0.5)*(cdf(N,sqrt(thresh)-A))-0.5*(cdf(N,-sqrt(thresh)-A))+0.5*cdf(N,-sqrt(thresh)-B)+0.5*(cdf(N,sqrt(thresh)-B,'upper'));
%for Y = (X+N+1)^2/4
%e = (0.5)*(cdf(N,sqrt(A2*thresh)-A-DC))-0.5*(cdf(N,-sqrt(A2*thresh)-A-DC))+0.5*cdf(N,-sqrt(A2*thresh)-B-DC)+0.5*(cdf(N,sqrt(A2*thresh)-B-DC,'upper'));
% e = -vpa(log10(e));


plot(thresh,e)
hold on





thresh2 = linspace(threshMin,threshMax,numThresh+1);
numBits = 1000000;
e2 = zeros(1,numThresh+1);

for i = 0:numThresh
    dataIn = randi([0 1],1,numBits);
    conditionX = (dataIn==0);
    dataIn(conditionX) = A;
    dataIn(~conditionX) = B;
    n = normrnd(mu,sigma,1,numBits);%mu,sigma,numrows,numcolumns
    ns = normrnd(mu,sigma,3,numBits);
    y = dataIn+n;
    % added SJS
    %yh=(y+DC).^2/A2;
    yh = y.^2+sum(ns.^2);
    y=yh;
    
    threshold = threshMin+threshStep*i;
    y2 = y;
    conditionY2 = (y2 >= threshold);
    y2(conditionY2) = A;
    y2(~conditionY2) = B;
    z = dataIn-y2;
    z = abs(z);
    e2(i+1) = sum(z)/numBits;
end

e2 = -log10(e2);
plot(thresh2,e2,'x')

firstinf = find(e2==inf);
if(isempty(firstinf))
    [M,I] = max(e2);
    firstinf = I;
end
e3 = e2(1:firstinf-1);
e3 = reshape(e3,length(e3),1);

thresh3 = reshape(linspace(threshMin,threshMin-threshStep+threshStep*length(e3),length(e3)),length(e3),1);
% ft3 = fittype('a*exp(b*x)+c');
% option3 = fitoptions(ft3);
% option3.Lower = [1 3 -inf];
% option3.Upper = [inf inf  0];
% f3 = fit(thresh3, e3, ft3,'StartPoint',[2 1 5])
% plot(f3)
ft3 = fittype('a*x^2+b*x+c');
% options3 = fitoptions(ft3);
% options3.Lower = [-Inf, -Inf];

f3 = fit(thresh3,e3,ft3,'StartPoint',[0 0 0],'lower',[0 -Inf -Inf]);
plot(f3)

lastinf = find(e2==inf,1,'last');
if(isempty(lastinf))
    [M,I] = max(e2);
    lastinf = I;
end
lastaccept = find(e2>1,1,'last');
if(isempty(lastaccept))
    lastaccept = length(e2);
end
%e4 = e2(lastinf+1:length(e2));
e4 = e2(lastinf+1:lastaccept);
e4 = reshape(e4,length(e4),1);
%thresh4 = reshape(linspace(threshMax+threshStep-threshStep*length(e4),threshMax,length(e4)),length(e4),1);
thresh4 = reshape(linspace(threshMax+threshStep-threshStep*(length(e4)+numThresh-lastaccept),threshMax-threshStep*(numThresh-lastaccept),length(e4)),length(e4),1);
ft4 = fittype('c*x^2+d*x+e');
f4 = fit(thresh4,e4,ft4,'StartPoint',[0 0 0]);
plot(f4)

grid
xlabel('threshold')
ylabel('-log10BER')
hold off    

fun = @(x) f3(x)-f4(x);
[M,I] = max(e);
theoxx = thresh(I);
theoyy = e(I);
estxx = fzero(fun,0);
estyy = -log10(10^(-f3(estxx))+10^(-f4(estxx)));
disp(theoxx)
disp(theoyy)
disp(estxx)
disp(estyy)