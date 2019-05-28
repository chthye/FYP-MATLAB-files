p1 = 0.5;
A = 1;
A2 = 1;
B = 0;
%DC = 0; %unused
numThresh = 100;
threshMin = B;
threshMax = A; %will be overwritten
threshStep = (threshMax-threshMin)/numThresh;%will be overwritten
mu = 0;
digits(40);

SNRstep = 0.5;
SNRpoints = 6;
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
    thresh = linspace(B,A,numThresh+1);

    %for Y = (X+N)^2
    % E = (0.5)*(cdf(N,sqrt(thresh)-A))-0.5*(cdf(N,-sqrt(thresh)-A))+0.5*cdf(N,-sqrt(thresh)-B)+0.5*(cdf(N,sqrt(thresh)-B,'upper'));
    e = (0.5)*(cdf(N,sqrt(thresh)-A))-0.5*(cdf(N,-sqrt(thresh)-A))+0.5*(chi2cdf(thresh/sigma^2,1,'upper'));

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
    [theoyy,optThresh] = max(e);

    %for Y = (X+N)^2
    % e = (0.5)*(cdf(N,sqrt(thresh)-A))-0.5*(cdf(N,-sqrt(thresh)-A))+0.5*cdf(N,-sqrt(thresh)-B)+0.5*(cdf(N,sqrt(thresh)-B,'upper'));
    %for Y = (X+N+1)^2/4
    %e = (0.5)*(cdf(N,sqrt(A2*thresh)-A-DC))-0.5*(cdf(N,-sqrt(A2*thresh)-A-DC))+0.5*cdf(N,-sqrt(A2*thresh)-B-DC)+0.5*(cdf(N,sqrt(A2*thresh)-B-DC,'upper'));
    % e = -vpa(log10(e));


    %plot(thresh,e)
    %hold on





    thresh2 = linspace(threshMin,threshMax,numThresh+1);
    numBits = 1000;
    e2 = zeros(1,numThresh+1);
    estyy = zeros(1,numTimes);
    for k = 1:numTimes
%         disp(k)
        for i = 0:numThresh
            dataIn = randi([0 1],1,numBits);
            conditionX = (dataIn==0);
            dataIn(conditionX) = A;
            dataIn(~conditionX) = B;
            n = normrnd(mu,sigma,1,numBits);%mu,sigma,numrows,numcolumns
            y = dataIn+n;
            % added SJS
            %yh=(y+DC).^2/A2;
            yh = y.^2;
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

        %plot(thresh2,e2,'o')

        firstzero = find(e2==0,1);
        if(isempty(firstzero))
            firstzero = optThresh;
        end
        % e33 = e22(2:firstzero-1);
        e33 = e22(2:firstzero-1);
        e33 = reshape(e33,length(e33),1);
        thresh3 = reshape(linspace(threshMin,threshMin-threshStep+threshStep*length(e33),length(e33)),length(e33),1);


        poly33coeff = polyfit(thresh3,e33,1);
        sigma33 = abs(1/(sqrt(2)*poly33coeff(1)));
        mu33 = -poly33coeff(2)*sqrt(2)*sigma33;



        lastzero = find(e2==0,1,'last');
        if(isempty(lastzero))
            lastzero = optThresh;
        end
        e44 = e22(lastzero+1:length(e22));
        e44 = reshape(e44,length(e44),1);

        thresh4 = reshape(linspace(threshMax+threshStep-threshStep*(length(e44)),threshMax,length(e44)),length(e44),1);




        poly44coeff = polyfit(thresh4,e44,1);
        sigma44 = abs(1/(sqrt(2)*poly44coeff(1)));
        mu44 = poly44coeff(2)*sqrt(2)*sigma44;

        e5 = 0.25*(erfc(abs(mu33-thresh2)/(sigma33*sqrt(2)))+erfc(abs(mu44-thresh2)/(sigma44*sqrt(2))));
        estyy(k) = min(e5);
        estyyTable(k,j+1) = estyy(k);
    end
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
