%From QPSK_W_EVT.m and BPSK_SNRvsBER_intervalv2.m

numBits = 1000;
numExtremes = 1000;
SNRstart = 13;
SNRstep = 1;
SNRpoints = 3;
numTimes = 25;
estBER = zeros(SNRpoints+1,1);
upperBER = zeros(SNRpoints+1,1);
lowerBER = zeros(SNRpoints+1,1);
theoBER = zeros(SNRpoints+1,1);
idealBER = zeros(SNRpoints+1,1);%nophaseerror
SNRspace = linspace(SNRstart,SNRstart+SNRstep*SNRpoints,SNRpoints+1);
for j = 0:SNRpoints
    %additive noise
    mu = 0;
    SNR = SNRstart+SNRstep*j



    sigma = sqrt(1/(2*(10^((SNR)/10))));
    %phase noise
    deltanu = 1;
    T = 0.0001;
    %sampling size for loop
    sampleSize = 45;
    if mod(sampleSize,2)==0
        sampleSize = sampleSize+1;
    end
    sigmaPhase = sqrt(2*pi*deltanu*T);
    %alpha for wiener filter
    wienerParm = (sigmaPhase^2+2*sigma^2-sigmaPhase*sqrt(sigmaPhase^2+4*sigma^2))/(2*sigma^2);
    wienerFiltCoeff = zeros(sampleSize,1);
    %calculate coeff as (1-a)*a^k for both sides
    for i = 1:sampleSize
       wienerFiltCoeff(i) = (1-wienerParm)*wienerParm^(abs(i-(sampleSize+1)/2)); 
    end
    wienerFiltCoeff = wienerFiltCoeff/sum(wienerFiltCoeff);
    
    %Ideal Theoretical
    eNoPhaseError = 0.25*erfc(sin(pi/4)/(sqrt(2)*sigma))+0.25*erfc(sin(pi/4)/(sqrt(2)*sigma));
    eNoPhaseError = -log10(eNoPhaseError);
    
    theoyy = zeros(numTimes,1);
    estyy = zeros(numTimes,1);
    %loop for confidence interval
    for j2 = 1:numTimes
        disp(j2)
        %loop to get extremes
        yLower = zeros(numExtremes,1);
        yUpper = zeros(numExtremes,1);
        for i = 1:numExtremes
%             disp(i)
            dataIn = randi([0 1], numBits, 1);
            %4QAM = QPSK
            M = 4;
            k = log2(M);
            dataInMatrix = reshape(dataIn, k, length(dataIn)/k)';
            dataSymbolsIn = bi2de(dataInMatrix);

            dataMod = qammod(dataSymbolsIn, M);

            txSignal = dataMod;
            txSignal = txSignal./sqrt(2);
            noise = normrnd(mu,sigma,length(txSignal),1);
            noise2 = normrnd(mu,sigma,length(txSignal),1);

            %Signal before phase noise
            rxSignal = txSignal + noise+1i*noise2;

            % scatterplot(rxSignal)

            %phase noise
            dphi = randn(numBits/k,1)*sqrt(2*pi*deltanu*T);
            phi = cumsum(dphi);

            rxSignalMag = abs(rxSignal);
            rxSignalAngle = angle(rxSignal);
            rxSignalAngle = rxSignalAngle + phi;

            %Received Signal
            rxSignal2 = rxSignalMag.*cos(rxSignalAngle) + rxSignalMag.*1i.*sin(rxSignalAngle);
            rxSignal3 = zeros(numBits/k,1);
            % pnoise = normrnd(mu,  ,length(txSignal),1);
            % rxSignal3 = rxSignal2.*exp(1i*pnoise);
            rxSignalAngle2 = angle(rxSignal2);


            estdphi = zeros(numBits/k,1);
            estphi = zeros(numBits/k,1);


            rxSignal3(1:sampleSize) = rxSignal2(1:sampleSize);
            %starting iteration
            for iteration = 1:(sampleSize-1)/2
                rxSignal3(1:iteration+(sampleSize-1)/2) = rxSignal2(1:iteration+(sampleSize-1)/2).*exp(-1i*estphi(max(1,iteration-1)));

               estdphi(iteration) =  mean(angle(-1*rxSignal3(1:iteration+(sampleSize-1)/2).^4))/4;
               if iteration>1
                   estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
               else
                   estphi(iteration) = estdphi(iteration);
               end
            end
            %main iteration
            for iteration = ((sampleSize+1)/2):(numBits/k-(sampleSize-1)/2)
            %     rxSignal3(iteration:iteration+(sampleSize-1)/2) = rxSignal2(iteration:iteration+(sampleSize-1)/2).*exp(-1i*estphi(iteration-1));
                rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2) = rxSignal2(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).*exp(-1i*estphi(iteration-1));
            %     estdphi(iteration) = sum(angle(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4).*wienerFiltCoeff)/4;
                estdphi(iteration) = angle((sum(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4.*wienerFiltCoeff))^(1/4));
                estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
            end
            %ending iteration
            for iteration = (numBits/k-(sampleSize-3)/2):numBits/k
            %     rxSignal3(iteration:numBits/k) = rxSignal2(iteration:numBits/k).*exp(-1i*estphi(iteration-1));
                rxSignal3(iteration-(sampleSize-1)/2:numBits/k) = rxSignal2(iteration-(sampleSize-1)/2:numBits/k).*exp(-1i*estphi(iteration-1));
                estdphi(iteration) =  mean(angle(-1*rxSignal3(iteration-(sampleSize-1)/2:numBits/k).^4))/4;
                estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
            end
            rxSignal3 = rxSignal2.*exp(-1i*estphi);
            rxNoiseAngle = angle(-rxSignal3.^4);

            %min and max deviations scaled from (-pi to pi), to (-1 to 1)
            yLower(i) = min(rxNoiseAngle)/pi;
            yUpper(i) = max(rxNoiseAngle)/pi;

        end

        
        %quasi-theoretical
        errorphi = estphi-phi;
        e = zeros(101,1);
        epsilon = linspace(-pi,pi-pi/200,400);
        alpha = 1/var(errorphi);

            threshold = 0;
            %phase error pdf (Tikhonov dist) (shifted by threshold before multiplying with BER given phase error)
        %     ptheta = exp(alpha*cos(linspace(-pi-threshold*pi/4,pi-pi/200-threshold*pi/4,400)))/(2*pi*besseli(0,alpha));
            %Gaussian approx
            ptheta = exp(-linspace(-pi-threshold*pi/4,pi-pi/200-threshold*pi/4,400).^2/(2*var(errorphi)))/(std(errorphi)*sqrt(2));
            ptheta = ptheta/sum(ptheta);%normalise probability
            e = 0;
            for eps = 51:150
                e = e + (0.5-0.25*erfc(sin(-epsilon(eps)-pi/4)/(sqrt(2)*sigma))+0.25*erfc(sin(3*pi/4+epsilon(eps))/(sqrt(2)*sigma))).*ptheta(eps);
            end
            for eps = 151:200
                e = e + (0.25*erfc(sin(pi/4+epsilon(eps))/(sqrt(2)*sigma))+0.25*erfc(sin(pi/4-epsilon(eps))/(sqrt(2)*sigma))).*ptheta(eps);
            end
            for eps = 201:250
                e = e + (0.25*erfc(sin(pi/4-epsilon(eps))/(sqrt(2)*sigma))+0.25*erfc(sin(pi/4+epsilon(eps))/(sqrt(2)*sigma))).*ptheta(eps);
            end
            for eps = 251:350
                e = e + (0.5-0.25*erfc(sin(epsilon(eps)-pi/4)/(sqrt(2)*sigma))+0.25*erfc(sin(3*pi/4-epsilon(eps))/(sqrt(2)*sigma))).*ptheta(eps);
            end

            


        %using gevfit on extremes
        yLower = -1 - yLower;
        yUpper = 1 - yUpper;

        lowerFx = linspace(1/(numExtremes+1),numExtremes/(numExtremes+1),numExtremes);
        lowerSpace = reshape(sort(yLower),1,numExtremes);

        upperFx = flip(linspace(1/(numExtremes+1),numExtremes/(numExtremes+1),numExtremes));
        upperSpace = reshape(sort(yUpper),1,numExtremes);

        upperSpace = 1 - upperSpace;
        [parmhatLower,parmciLower] = gevfit(lowerSpace);
        [parmhatUpper,parmciUpper] = gevfit(upperSpace);

        thresh2 = linspace(-1,1,101);

        lowerVal = vpa(1-gevcdf(thresh2,parmhatLower(1),parmhatLower(2),parmhatLower(3)).^(1/numBits));
        upperVal = vpa(1-gevcdf(1-thresh2,parmhatUpper(1),parmhatUpper(2),parmhatUpper(3)).^(1/numBits));
        BER = (lowerVal + upperVal);

    
        theoyy(j2) = e;
        estyy(j2) = min(BER);
        
    end
    
    idealBER(j+1) = eNoPhaseError;
    theoBER(j+1) = -log10(mean(theoyy));
    estBER(j+1) = mean(estyy);
    estBER(j+1) = -log10(estBER(j+1));
    lowerBER(j+1) = -log10(min(estyy));
    upperBER(j+1) = -log10(max(estyy));

end

%     e2 = -log10(BER);
%     figure(1)
%     plot(thresh2,e2)
%     hold on
%     plot(thresh2,e,'--')
%     plot(thresh2,eNoPhaseError,'.')
%     hold off
%     xlabel('threshold in \pi')
%     ylabel('-log_1_0BER')
%     grid
%     display(e2(51))
plot(SNRspace,idealBER,SNRspace,theoBER,'--k',SNRspace,estBER,'x',SNRspace,lowerBER,'.',SNRspace,upperBER,'.')
grid
xlabel('SNR (dB)')
ylabel('-log10BER')