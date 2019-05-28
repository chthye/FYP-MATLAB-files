%RRC section code from https://uk.mathworks.com/help/comm/gs/compute-ber-for-a-qam-system-with-awgn-using-matlab.html
%changed demod part for different threshold


%root raised cosine filter
rolloff = 0.8;
numSamplesPerSymbol = 8;        
span = 4;
rrcFilter = rcosdesign(rolloff,span,numSamplesPerSymbol);
%fvtool(rrcFilter,'Analysis','Impulse') %to visualise filter


p1 = 0.5;
A = 1;
B = -1;
mu = 0;
EbNostep = 0.5;
EbNopoints = 10;
numTimes = 25;
estBER = zeros(EbNopoints+1,1);
upperBER = zeros(EbNopoints+1,1);
lowerBER = zeros(EbNopoints+1,1);
theoBER = zeros(EbNopoints+1,1);
EbNospace = linspace(10,10+EbNostep*EbNopoints,EbNopoints+1);

for j = 0:EbNopoints
    
    EbNo = 10+EbNostep*j
    sigma = sqrt(1/(2*(10^((EbNo)/10))));
    %using definition SNR = 10*log10(Eb/No), and variance = No/2
    digits(40);
    N = makedist('Normal','mu',mu,'sigma',sigma);
    thresh = linspace(B,A,2001);


    e = (0.5)*(cdf(N,thresh-B,'upper'))+0.5*(cdf(N,thresh-A));
    %note: upper returns complement of cdf, much more accurate than 1-cdf
    e = -vpa(log10(e));
%     plot(thresh,e)
%     hold on

    %bits for transmission
    numBits = 400000;
    %threshold and error matrices
    e2 = zeros(101,1);
    thresh2 = linspace(-1,1,101);
    %estimated BER matrix for each EbNo
    estyy = zeros(1,numTimes);
    
    for j2 = 1:numTimes
        dataIn = randi([0 1], numBits, 1);
        %plot(linspace(1,numBits,numBits),dataIn) %plot bits out for test

        %group bits to integers for QAM. M=2 means no regroup, BPSK
        M = 2;
        k = log2(M);
        dataInMatrix = reshape(dataIn, length(dataIn)/k, k);
        dataSymbolsIn = bi2de(dataInMatrix);
        %plot(linspace(1,length(dataSymbolsIn),length(dataSymbolsIn)),dataSymbolsIn)

        %QAM modulation
        dataMod = qammod(dataSymbolsIn, M);
        %plot(linspace(1,length(dataMod),length(dataMod)),dataMod)



        %upsample and apply root raised cosine filter
        txSignal = upfirdn(dataMod, rrcFilter, numSamplesPerSymbol, 1);
        % figure(1);
        % plot(linspace(1,length(txSignal),length(txSignal)),txSignal)
        % hold on


        %set snr(AWGN)
        %EbNo = 10;
        snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol) + 10*log10(2);

        %note: EbNo,snr in dB. 
        %snr = Psignal/Pnoise = (Esignal*symbolrate)/(noisePSD*bandwidth)
        %Es = Eb*log2(M),noisePSD = No??,  ???bandwidth=fs/2???check simulation

        %pass through channel, to receiver
        rxSignal = awgn(txSignal, snr, 'measured');
        %plot(linspace(1,length(rxSignal),length(rxSignal)),rxSignal)
        %hold off

        % Downsample and filter
        rxFiltSignal = upfirdn(rxSignal,rrcFilter,1,numSamplesPerSymbol);
        % Account for delay
        rxFiltSignal = rxFiltSignal(span+1:end-span);

        %own decoder for bpsk from here

        for n = 0:100
            rxDecodedSignal = rxFiltSignal;
            threshold = -1+0.02*n;
            decodeCondition = (rxDecodedSignal >= threshold);
            rxDecodedSignal(decodeCondition) = 1;
            rxDecodedSignal(~decodeCondition) = -1;
            %demodulate
            demodCondition = (rxDecodedSignal == -1);
            rxDecodedSignal(demodCondition) = 0;
            %e2(n+1) = biterr(dataIn, rxDecodedSignal)/numBits;
            e2(n+1) = sum(abs(dataIn-rxDecodedSignal))/(numBits);
        end
        % figure(2);
        e2 = -log10(e2);
    %     plot(thresh2,e2,'x')
        % hold off

        firstinf = find(e2==inf);
        if(isempty(firstinf))
            firstinf = 50;
        end
        e3 = e2(1:firstinf-1);
        e3 = reshape(e3,length(e3),1);

        thresh3 = reshape(linspace(-1,-1.02+0.02*length(e3),length(e3)),length(e3),1);

%         ft3 = fittype('a*x^2+b*x+c');
%         f3 = fit(thresh3,e3,ft3,'StartPoint',[0 0 0]);
    %     plot(f3)
        poly3coeff = polyfit(thresh3,e3,2);
        poly3y = polyval(poly3coeff, thresh2);

        lastinf = find(e2==inf,1,'last');
        if(isempty(lastinf))
            lastinf = 50;
        end
        e4 = e2(lastinf+1:length(e2));
        e4 = reshape(e4,length(e4),1);
        thresh4 = reshape(linspace(1.02-0.02*length(e4),1,length(e4)),length(e4),1);
%         ft4 = fittype('a*x^2+b*x+c');
%         f4 = fit(thresh4,e4,ft4,'StartPoint',[0 0 0]);
    %     plot(f4)
        poly4coeff = polyfit(thresh4,e4,2);
        poly4y = polyval(poly4coeff, thresh2);

    %     grid
    %     xlabel('threshold')
    %     ylabel('-log10BER')
    %     hold off

%         fun = @(x) f3(x)-f4(x);
        theoyy = e(1001);
%         estxx = fzero(fun,0);
%         estyy(j2) = -log10(10^(-f3(estxx))+10^(-f4(estxx)));
        estxx = (-(poly3coeff(2)-poly4coeff(2))+sqrt((poly3coeff(2)-poly4coeff(2))^2-4*(poly3coeff(1)-poly4coeff(1))*(poly3coeff(3)-poly4coeff(3))))/(2*(poly3coeff(1)-poly4coeff(1)));
        if(estxx > A)
            estxx = (-(poly3coeff(2)-poly4coeff(2))-sqrt((poly3coeff(2)-poly4coeff(2))^2-4*(poly3coeff(1)-poly4coeff(1))*(poly3coeff(3)-poly4coeff(3))))/(2*(poly3coeff(1)-poly4coeff(1)));
        end
        estyy(j2) = -log10(10^(-polyval(poly3coeff,estxx))+10^(-polyval(poly4coeff,estxx)));
    %     disp(theoyy)
    %     disp(estxx)
    %     disp(estyy)
        


    end
    theoBER(j+1) = theoyy;
    estBER(j+1) = mean(estyy);
    lowerBER(j+1) = estBER(j+1) - 1.96*std2(estyy)/numTimes;
    upperBER(j+1) = estBER(j+1) + 1.96*std2(estyy)/numTimes;
end
lowerBER2 = smooth(lowerBER,0.5,'loess');
upperBER2 = smooth(upperBER,0.5,'loess');
plot(EbNospace,theoBER,EbNospace,estBER,'x',EbNospace,lowerBER,'.',EbNospace,upperBER,'.',EbNospace,lowerBER2,EbNospace,upperBER2)
%plot(EbNospace,theoBER,EbNospace,estBER,'x')
grid
xlabel('EbNo (dB)')
ylabel('-log10BER')

