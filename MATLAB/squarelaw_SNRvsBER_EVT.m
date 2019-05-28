threshMin = B;
threshMax = A; %will be overwritten
threshStep = (threshMax-threshMin)/numThresh;%will be overwritten
mu = 0;
digits(40);

SNRstep = 0.5;
SNRpoints = 4;
numTimes = 50;
estyyTable = zeros(numTimes,SNRpoints+1);
estBER = zeros(SNRpoints+1,1);
upperBER = zeros(SNRpoints+1,1);
lowerBER = zeros(SNRpoints+1,1);
theoBER = zeros(SNRpoints+1,1);
SNRspace = linspace(14,14+SNRstep*SNRpoints,SNRpoints+1);

for j = 0:SNRpoints
    SNR = 14+SNRstep*j
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
    theoyy = max(e);
    % plot(thresh,e)
    % hold on


    thresh2 = linspace(threshMin,threshMax,numThresh+1);
    numBits = 100000;
    numExtremes = 1000;
    yLower = zeros(numExtremes,1);
    yUpper = zeros(numExtremes,1);
    e2 = zeros(1,numThresh+1);

    estyy = zeros(1,numTimes);
        for k = 1:numTimes
            disp(k)
            for i = 1:numExtremes
%                 dataIn = randi([0 1],1,numBits);
                dataIn = [zeros(1,numBits/2) ones(1,numBits/2)];
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
                
                
%                 y2 = sort(y);
%                 indexLower = find(y2<optThresh,1,'last');
%                 yLower(i) = y2(indexLower);
%                 indexUpper = find(y2>optThresh,1);
%                 yUpper(i) = y2(indexUpper);
                yLower(i) = max(y(numBits/2+1:numBits));
                yUpper(i) = min(y(1:numBits/2));



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
            estyy(k) = min(BER);
            estyyTable(k,j+1) = estyy(k);
    %         e2 = -log10(BER);
    end
% plot(thresh2,e2)
% display(max(e2))
    theoBER(j+1) = theoyy;
    estBER(j+1) = mean(estyy);
    estBER(j+1) = -log10(estBER(j+1));
    lowerBER(j+1) = -log10(min(estyy));
    upperBER(j+1) = -log10(max(estyy));

end
plot(SNRspace,theoBER,'k',SNRspace,estBER,'x',SNRspace,lowerBER,'.',SNRspace,upperBER,'.')
grid
xlabel('SNR (dB)')
ylabel('-log10BER')
% hold off
% grid
