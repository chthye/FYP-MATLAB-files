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

SNR = 15;
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
%get optimum threshold for EVT use
[a,optThresh] = min(e);
optThresh = (optThresh-1)/numThresh;
%for Y = (X+N+DC)^2/4
%e = (0.5)*(cdf(N,sqrt(A2*thresh)-A-DC))-0.5*(cdf(N,-sqrt(A2*thresh)-A-DC))+0.5*cdf(N,-sqrt(A2*thresh)-B-DC)+0.5*(cdf(N,sqrt(A2*thresh)-B-DC,'upper'));
%note: upper returns complement of cdf, much more accurate than 1-cdf
e = -log10(e);

plot(thresh,e)
hold on


thresh2 = linspace(threshMin,threshMax,numThresh+1);
numBits = 1000;
numExtremes = 1000;
yLower = zeros(numExtremes,1);
yUpper = zeros(numExtremes,1);
e2 = zeros(1,numThresh+1);

for i = 1:numExtremes
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
    y2 = sort(y);
    indexLower = find(y2<optThresh,1,'last');
    yLower(i) = y2(indexLower);
    indexUpper = find(y2>optThresh,1);
    yUpper(i) = y2(indexUpper);
    
    
    
    
end

lowerFx = linspace(1/(numExtremes+1),numExtremes/(numExtremes+1),numExtremes);
lowerSpace = reshape(sort(yLower),1,numExtremes);
upperFx = flip(linspace(1/(numExtremes+1),numExtremes/(numExtremes+1),numExtremes));
upperSpace = reshape(sort(yUpper),1,numExtremes);
upperSpace = 1 - upperSpace;

[parmhatLower,parmciLower] = gevfit(lowerSpace);
[parmhatUpper,parmciUpper] = gevfit(upperSpace);


thresh2 = linspace(B,A,101);

lowerVal = vpa(1-gevcdf(thresh2,parmhatLower(1),parmhatLower(2),parmhatLower(3)).^(1/numBits));
upperVal = vpa(1-gevcdf(1-thresh2,parmhatUpper(1),parmhatUpper(2),parmhatUpper(3)).^(1/numBits));
BER = (lowerVal + upperVal);

e2 = -log10(BER);
plot(thresh2,e2)
display(max(e2))


hold off
grid
