p1 = 0.5;
A = 1;
B = -1;
mu = 0;




SNR = 12;
sigma = sqrt(1/(2*(10^((SNR)/10))));
%using definition SNR = 10*log10(Eb/No), and variance = No/2
digits(40);
N = makedist('Normal','mu',mu,'sigma',sigma);
thresh = linspace(B,A,2001);

e = (0.5)*(cdf(N,thresh-B,'upper'))+0.5*(cdf(N,thresh-A));
%note: upper returns complement of cdf, much more accurate than 1-cdf
e1 = -vpa(log10(e));
plot(thresh,e1)
hold on
display(e1(1001))
numBins = 500;
binStep = 0.5*(A-B)/numBins;
numBits = 10000;
numExtremes = 1000;
lowerBinEdges = linspace(-1,0,numBins+1);
upperBinEdges = linspace(0,1,numBins+1);
yLower = zeros(numExtremes,1);
yUpper = zeros(numExtremes,1);
for i = 1:numExtremes
            dataIn = randi([0 1],1,numBits);
            conditionX = (dataIn==0);
            dataIn(conditionX) = A;
            dataIn(~conditionX) = B;
            n = normrnd(mu,sigma,1,numBits);%mu,sigma,numrows,numcolumns
            y = dataIn+n;
            y2 = sort(y);
%             indexLower = find(y2<0,1,'last');
            indexLower = 1;
            yLower(i) = y2(indexLower);
%             indexUpper = find(y2>0,1);
            indexUpper = numBits;
            yUpper(i) = y2(indexUpper);


end
yLower = -2 - yLower;
yUpper = 2 - yUpper;
% lowerHist = histcounts(yLower,lowerBinEdges);
% lowerHist = lowerHist/sum(lowerHist);
% upperHist = histcounts(yUpper,upperBinEdges);
% upperHist = upperHist/sum(upperHist);
% fullHist = [lowerHist upperHist];
lowerFx = linspace(1/(numExtremes+1),numExtremes/(numExtremes+1),numExtremes);
lowerSpace = reshape(sort(yLower),1,numExtremes);
% lowerFx = cumsum(lowerHist);
% lowerFirstIndex = find(lowerFx>0,1);
% lowerLastIndex = find(lowerFx<0.9999,1,'last');
% lowerSpace = linspace(B-0.5*binStep+binStep*lowerFirstIndex,B-0.5*binStep+binStep*lowerLastIndex,lowerLastIndex-lowerFirstIndex+1);

upperFx = flip(linspace(1/(numExtremes+1),numExtremes/(numExtremes+1),numExtremes));
upperSpace = reshape(sort(yUpper),1,numExtremes);
% upperHistFlipped = flip(upperHist);
% upperFxFlipped = cumsum(upperHistFlipped);
% upperFx = flip(upperFxFlipped);
% upperFirstIndex = find(upperFx<0.9999,1);
% upperLastIndex = find(upperFx>0,1,'last');
% upperSpace = linspace(-0.5*binStep+binStep*upperFirstIndex,-0.5*binStep+binStep*upperLastIndex,upperLastIndex-upperFirstIndex+1);

% lowerLogFx = -log(-log(lowerFx));
% upperLogFx = -log(-log(upperFx));

upperSpace = 1 - upperSpace;
[parmhatLower,parmciLower] = gevfit(lowerSpace);
[parmhatUpper,parmciUpper] = gevfit(upperSpace);

% lowerCoeff = polyfit(lowerSpace,lowerLogFx(lowerFirstIndex:lowerLastIndex),1);
% upperCoeff = polyfit(upperSpace,upperLogFx(upperFirstIndex:upperLastIndex),1);

% lowerCoeff = polyfit(lowerSpace,lowerLogFx,1);
% upperCoeff = polyfit(upperSpace,upperLogFx,1);


thresh2 = linspace(B,A,101);

% lowerVal = polyval(lowerCoeff,thresh2);
% upperVal = polyval(upperCoeff,thresh2);

lowerVal = vpa(1-gevcdf(thresh2,parmhatLower(1),parmhatLower(2),parmhatLower(3)).^(1/numBits));
upperVal = vpa(1-gevcdf(1-thresh2,parmhatUpper(1),parmhatUpper(2),parmhatUpper(3)).^(1/numBits));
BER = (lowerVal + upperVal);

% BER = (exp(-lowerVal)/numBits + exp(-upperVal)/numBits);
e2 = -log10(BER);
plot(thresh2,e2)
display(e2(51))


hold off
grid


% pd2 = makedist('HalfNormal','mu',mu,'sigma',sigma);
% pdf2 = pdf(pd2,linspace(0,1,1001));
% pdf2 = (pdf2/sum(pdf2));
% cdf2 = (cumsum(pdf2));
% cdf3 = cdf2.^numBits;
% pdf3 = diff(cdf3);
% plot(linspace(-1+binStep/2,-binStep/2,numBins),lowerHist)
% hold on
% plot(linspace(-0.9995,-0.0005,1000),pdf3)
% 
% hold off



