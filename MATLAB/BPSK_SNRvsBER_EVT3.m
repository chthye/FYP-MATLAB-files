p1 = 0.5;
A = 1;
B = -1;
mu = 0;

SNRstep = 0.5;
SNRpoints = 7;
numTimes = 9;
estBER = zeros(SNRpoints+1,1);
upperBER = zeros(SNRpoints+1,1);
lowerBER = zeros(SNRpoints+1,1);
theoBER = zeros(SNRpoints+1,1);
SNRspace = linspace(10,10+SNRstep*SNRpoints,SNRpoints+1);


for j = 0:SNRpoints
    SNR = 10+SNRstep*j

    sigma = sqrt(1/(2*(10^((SNR)/10))));
    %using definition SNR = 10*log10(Eb/No), and variance = No/2
    digits(40);
    N = makedist('Normal','mu',mu,'sigma',sigma);
    thresh = linspace(B,A,2001);

    e = (0.5)*(cdf(N,thresh-B,'upper'))+0.5*(cdf(N,thresh-A));
    %note: upper returns complement of cdf, much more accurate than 1-cdf
    e1 = -log10(e);
    % plot(thresh,e1)
    % hold on
    % display(e1(1001))


    numBins = 500;
    binStep = 0.5*(A-B)/numBins;
    numBits = 10000;
    numExtremes = 10000;
    lowerBinEdges = linspace(-1,0,numBins+1);
    upperBinEdges = linspace(0,1,numBins+1);
    estyy = zeros(1,numTimes);
    for k = 1:numTimes
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
                    indexLower = find(y2<0,1,'last');
                    yLower(i) = y2(indexLower);
                    indexUpper = find(y2>0,1);
                    yUpper(i) = y2(indexUpper);


        end

%         lowerHist = histcounts(yLower,lowerBinEdges);
%         lowerHist = lowerHist/sum(lowerHist);
%         upperHist = histcounts(yUpper,upperBinEdges);
%         upperHist = upperHist/sum(upperHist);
% 
%         lowerFx = cumsum(lowerHist);
%         lowerFirstIndex = find(lowerFx>0,1);
%         lowerLastIndex = find(lowerFx<0.9999,1,'last');
%         lowerSpace = linspace(B-0.5*binStep+binStep*lowerFirstIndex,B-0.5*binStep+binStep*lowerLastIndex,lowerLastIndex-lowerFirstIndex+1);
% 
% 
%         upperHist = flip(upperHist);
%         upperFx = cumsum(upperHist);
%         upperFx = flip(upperFx);
%         upperFirstIndex = find(upperFx<0.9999,1);
%         upperLastIndex = find(upperFx>0,1,'last');
%         upperSpace = linspace(-0.5*binStep+binStep*upperFirstIndex,-0.5*binStep+binStep*upperLastIndex,upperLastIndex-upperFirstIndex+1);
% 
%         lowerLogFx = -log(-log(lowerFx));
%         upperLogFx = -log(-log(upperFx));
% 
%         lowerCoeff = polyfit(lowerSpace,lowerLogFx(lowerFirstIndex:lowerLastIndex),1);
%         upperCoeff = polyfit(upperSpace,upperLogFx(upperFirstIndex:upperLastIndex),1);

        % thresh2 = linspace(B,A,101);


        
        lowerFx = linspace(1/(numExtremes+1),numExtremes/(numExtremes+1),numExtremes);
        lowerSpace = reshape(sort(yLower),1,numExtremes);
        upperFx = flip(linspace(1/(numExtremes+1),numExtremes/(numExtremes+1),numExtremes));
        upperSpace = reshape(sort(yUpper),1,numExtremes);
        
        upperSpace = 1 - upperSpace;
        [parmhatLower,parmciLower] = gevfit(lowerSpace);
        [parmhatUpper,parmciUpper] = gevfit(upperSpace);
        
        lowerVal = vpa(1-gevcdf(0,parmhatLower(1),parmhatLower(2),parmhatLower(3)).^(1/numBits));
        upperVal = vpa(1-gevcdf(1,parmhatUpper(1),parmhatUpper(2),parmhatUpper(3)).^(1/numBits));
        BER = vpa(lowerVal + upperVal);
        
        
        
        
%         BER = exp(-lowerVal)/numBits + exp(-upperVal)/numBits;
        e2 = vpa(-log10(BER));
        % plot(thresh2,e2)
        % display(e2(51))

        theoyy = e1(1001);
        estyy(k) = e2



    end

    % hold off
    theoBER(j+1) = theoyy;
    estBER(j+1) = mean(estyy(isfinite(estyy)));
    lowerBER(j+1) = estBER(j+1) - 1.96*std2(estyy(isfinite(estyy)))/sqrt(length(estyy(isfinite(estyy))));
    upperBER(j+1) = estBER(j+1) + 1.96*std2(estyy(isfinite(estyy)))/sqrt(length(estyy(isfinite(estyy))));
end
lowerBER2 = smooth(lowerBER,0.5,'loess');
upperBER2 = smooth(upperBER,0.5,'loess');
plot(SNRspace,theoBER,SNRspace,estBER,'x',SNRspace,lowerBER,'.',SNRspace,upperBER,'.',SNRspace,lowerBER2,SNRspace,upperBER2)
grid
xlabel('SNR (dB)')
ylabel('-log10BER')


