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

e22 = erfcinv(4*e2);
e222 = -log10(e2);
plot(thresh2,e222,'x')


firstzero = find(e2==0,1);
if(isempty(firstzero))
    firstzero = 50;
end
% e33 = e22(2:firstzero-1);
e33 = e2(2:firstzero-1);
e33 = chi2inv(1-e33,4);

e33 = reshape(e33,length(e33),1);
thresh3 = reshape(linspace(threshMin,threshMin-threshStep+threshStep*length(e33),length(e33)),length(e33),1);


poly33coeff = polyfit(thresh3,e33,1);
% sigma33 = abs(1/(sqrt(2)*poly33coeff(1)));
% mu33 = -poly33coeff(2)*sqrt(2)*sigma33;
sigma33 = abs(1/sqrt(poly33coeff(1)));


lastzero = find(e2==0,1,'last');
if(isempty(lastzero))
    lastzero = 50;
end
e44 = e22(lastzero+1:length(e22));
e44 = reshape(e44,length(e44),1);

thresh4 = reshape(linspace(threshMax+threshStep-threshStep*(length(e44)),threshMax,length(e44)),length(e44),1);




poly44coeff = polyfit(thresh4,e44,1);
sigma44 = abs(1/(sqrt(2)*poly44coeff(1)));
mu44 = poly44coeff(2)*sqrt(2)*sigma44;

% e5 = -log10(0.25*(erfc(abs(mu33-thresh2)/(sigma33*sqrt(2)))+erfc(abs(mu44-thresh2)/(sigma44*sqrt(2)))));
e5 = -log10(0.5*(chi2cdf(thresh2/sigma33^2,4,'upper'))+0.25*(erfc(abs(mu44-thresh2)/(sigma44*sqrt(2)))));

plot(thresh2,e5)

grid
xlabel('threshold')
ylabel('-log10BER')
hold off

[theoyy,theoxx] = max(e);
theoxx = thresh(theoxx);

[estyy2,estxx2] = max(e5);
estxx2 = thresh2(estxx2);

disp(theoxx)
disp(theoyy)

disp(estxx2)
disp(estyy2)