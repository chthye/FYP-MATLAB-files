p1 = 0.5;
A = 1;
B = -1;
mu = 0;
SNRstep = 1;
SNRpoints = 6;
numTimes = 100;
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
    e = -log10(e);
    %plot(thresh,e)
    %hold on

    thresh2 = linspace(B,A,101);
    numBits = 100000;
    e2 = zeros(1,101);
    estyy = zeros(1,numTimes);
    for k = 1:numTimes
        for i = 0:100
            dataIn = randi([0 1],1,numBits);
            conditionX = (dataIn==0);
            dataIn(conditionX) = A;
            dataIn(~conditionX) = B;
            n = normrnd(mu,sigma,1,numBits);%mu,sigma,numrows,numcolumns
            y = dataIn+n;
            % added SJS
            %yh=(y+1).^2-1;
            %y=yh;

            threshold = B+0.01*(A-B)*i;
            y2 = y;
            conditionY2 = (y2 >= threshold);
            y2(conditionY2) = A;
            y2(~conditionY2) = B;
            z = dataIn-y2;
            z = abs(z);
            e2(i+1) = 0.5*sum(z)/numBits;
        end

        e2 = -log10(e2);
        %plot(thresh2,e2,'o')

        firstinf = find(e2==inf);
        if(isempty(firstinf))
            firstinf = 50;
        end
        e3 = e2(1:firstinf-1);
        e3 = reshape(e3,length(e3),1);

        thresh3 = reshape(linspace(-1,-1.02+0.02*length(e3),length(e3)),length(e3),1);
        % ft3 = fittype('a*exp(b*x)+c');
        % option3 = fitoptions(ft3);
        % option3.Lower = [1 3 -inf];
        % option3.Upper = [inf inf  0];
        % f3 = fit(thresh3, e3, ft3,'StartPoint',[2 1 5])
        % plot(f3)
        
%         ft3 = fittype('a*x^2+b*x+c');
%         f3 = fit(thresh3,e3,ft3,'StartPoint',[0 0 0]);
%         p_seb= polyfit(thresh3,e3,2);
%         p2_seb= polyval(p_seb, thresh3);
        poly3coeff = polyfit(thresh3,e3,2);
        poly3y = polyval(poly3coeff, thresh2);
        %plot(f3)

        lastinf = find(e2==inf,1,'last');
        if(isempty(lastinf))
            lastinf = 50;
        end
        e4 = e2(lastinf+1:length(e2));
        e4 = reshape(e4,length(e4),1);
        thresh4 = reshape(linspace(1.02-0.02*length(e4),1,length(e4)),length(e4),1);
%         ft4 = fittype('a*x^2+b*x+c');
%         f4 = fit(thresh4,e4,ft4,'StartPoint',[0 0 0]);
%                 p_seb= polyfit(thresh3,e3,2);
%         p2_seb= polyval(p_seb, thresh3);
        poly4coeff = polyfit(thresh4,e4,2);
        poly4y = polyval(poly4coeff, thresh2);

        %plot(f4)

        %grid
        %xlabel('threshold')
        %ylabel('-log10BER')
        %hold off

%         fun = @(x) f3(x)-f4(x);
        theoyy = e(1001);
%         estxx = fzero(fun,0);
        %estyy = f3(estxx)
%         estyy(k) = -log10(10^(-f3(estxx))+10^(-f4(estxx)));

        %quadratic
        estxx = (-(poly3coeff(2)-poly4coeff(2))+sqrt((poly3coeff(2)-poly4coeff(2))^2-4*(poly3coeff(1)-poly4coeff(1))*(poly3coeff(3)-poly4coeff(3))))/(2*(poly3coeff(1)-poly4coeff(1)));
        if(estxx > A)
            estxx = (-(poly3coeff(2)-poly4coeff(2))-sqrt((poly3coeff(2)-poly4coeff(2))^2-4*(poly3coeff(1)-poly4coeff(1))*(poly3coeff(3)-poly4coeff(3))))/(2*(poly3coeff(1)-poly4coeff(1)));
        end
        
        %cubic
%         diffOfCoeff = [poly3coeff(1)]

%         estyy(k) = -log10(10^(-polyval(poly3coeff,estxx))+10^(-polyval(poly4coeff,estxx)));
        estyy(k) = (10^(-polyval(poly3coeff,estxx))+10^(-polyval(poly4coeff,estxx)));
    end
    theoBER(j+1) = theoyy;
%     estBER(j+1) = mean(estyy);
%     lowerBER(j+1) = estBER(j+1) - 1.96*std2(estyy)/sqrt(numTimes);
%     upperBER(j+1) = estBER(j+1) + 1.96*std2(estyy)/sqrt(numTimes);
%     lowerBER(j+1) = min(estyy);
%     upperBER(j+1) = max(estyy);
    estBER(j+1) = -log10(mean(estyy));
    lowerBER(j+1) = -log10(min(estyy));
    upperBER(j+1) = -log10(max(estyy));

end
% lowerBER2 = smooth(lowerBER,0.5,'loess');
% upperBER2 = smooth(upperBER,0.5,'loess');
plot(SNRspace,theoBER,SNRspace,estBER,'x',SNRspace,lowerBER,'x',SNRspace,upperBER,'x')
grid
xlabel('SNR (dB)')
ylabel('-log10BER')