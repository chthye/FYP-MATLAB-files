%from QPSK_W_CScount_ver2.m
%instead of loop across phase window, loop across SNR

numBits = 1000000;
bitErrorWindow = 10;
CStotal = zeros(1,26);
CSmeantime = zeros(1,26);
CSprobs2 = zeros(1,26);
CSlogprobs2 = zeros(1,26);
for iter = 26
    disp(iter)
    %additive noise
    mu = 0;
    SNR = 10+0.2*(iter-1);
    sigma = sqrt(1/(2*(10^((SNR)/10))));
    %phase noise
    deltanu = 1;
    T = 0.0030;
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


    phaseWindow = 2*pi;

    dataIn = randi([0 1], numBits, 1);
    %4QAM = QPSK
    M = 4;
    k = log2(M);
    dataInMatrix = reshape(dataIn, k, length(dataIn)/k)';
    dataSymbolsIn = bi2de(dataInMatrix);

    dataMod = qammod(dataSymbolsIn, M);

    txSignal = dataMod;
    %normalise amplitude
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

    rxSignalAngle2 = angle(rxSignal2);


    estdphi = zeros(numBits/k,1);
    estphi = zeros(numBits/k,1);


    rxSignal3(1:sampleSize) = rxSignal2(1:sampleSize);
    %starting iteration
    for iteration = 1:(sampleSize-1)/2
        rxSignal3(1:iteration+(sampleSize-1)/2) = rxSignal2(1:iteration+(sampleSize-1)/2).*exp(-1i*estphi(max(1,iteration-1)));
        phases = angle(-1*rxSignal3(1:iteration+(sampleSize-1)/2).^4);
        phases(phases > phaseWindow/2) = phases(phases>phaseWindow/2) - phaseWindow;
        phases(phases < -phaseWindow/2) = phases(phases < -phaseWindow/2) + phaseWindow;

       estdphi(iteration) =  mean(phases)/4;
       if iteration>1
           estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
       else
           estphi(iteration) = estdphi(iteration);
       end
    end
    %main iteration
    for iteration = ((sampleSize+1)/2):(numBits/k-(sampleSize-1)/2)
        rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2) = rxSignal2(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).*exp(-1i*estphi(iteration-1));
        phases = angle(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4);
        phases(phases > phaseWindow/2) = phases(phases>phaseWindow/2) - phaseWindow;
        phases(phases < -phaseWindow/2) = phases(phases < -phaseWindow/2) + phaseWindow;

%         estdphi(iteration) =  mean(phases)/4;
%         estdphi(iteration) = sum(phases.*wienerFiltCoeff)/4;
        estdphi(iteration) = angle((sum(abs(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4).*exp(1i.*phases).*wienerFiltCoeff))^(1/4));

        estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
    end
    %ending iteration
    for iteration = (numBits/k-(sampleSize-3)/2):numBits/k
        rxSignal3(iteration:numBits/k) = rxSignal2(iteration:numBits/k).*exp(-1i*estphi(iteration-1));
        phases = angle(-1*rxSignal3(iteration-(sampleSize-1)/2:numBits/k).^4);
        phases(phases > phaseWindow/2) = phases(phases>phaseWindow/2) - phaseWindow;
        phases(phases < -phaseWindow/2) = phases(phases < -phaseWindow/2) + phaseWindow;

        estdphi(iteration) =  mean(phases)/4;
        estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
    end
    rxSignal3 = rxSignal2.*exp(-1i*estphi);


    
    dataSymbolsOut = qamdemod(rxSignal3,M);
    dataOutMatrix = de2bi(dataSymbolsOut,k);
    % dataOut = dataOutMatrix(:);
    dataOut = reshape(dataOutMatrix',numBits,1);
    disp(biterr(dataIn,dataOut)/numBits);

    %Cycle slip counter ver 2 starts here
    dataSymbolsOut = qamdemod(rxSignal3,M);
    dataOutMatrix = de2bi(dataSymbolsOut,k);
    % dataOut = dataOutMatrix(:);
    dataOut = reshape(dataOutMatrix',numBits,1);

    bitError = abs(dataIn-dataOut);
    bitErrorSum = cumsum(bitError);
    bitError2 = reshape(bitError,bitErrorWindow,numBits/bitErrorWindow);
    bitErrorSum2 = sum(bitError2)';
    CSindicator2 = 0;
    CStotal2 = 0;
    for iteration = 1:length(bitErrorSum2)
        switch CSindicator2
            case 0
                if bitErrorSum2(iteration)>bitErrorWindow*0.4
                    CStotal2 = CStotal2+1;
                    CSindicator2 = 1;
    %                 disp(iteration)
                end
            case 1
                if bitErrorSum2(iteration)>bitErrorWindow*0.9
                    CStotal2 = CStotal2+1;
                    CSindicator2 = 2;
    %                 disp(iteration)
                elseif bitErrorSum2(iteration)<bitErrorWindow*0.1
                    CStotal2 = CStotal2+1;
                    CSindicator2 = 0;
    %                 disp(iteration)
                end
            case 2
                if bitErrorSum2(iteration) < bitErrorWindow*0.6
                    CStotal2 = CStotal2+1;
                    CSindicator2 = 1;
    %                 disp(iteration)
                end
        end
    end
    CSmeantime2 = numBits/(k*CStotal2);
%     disp(CStotal2)
%     disp(CSmeantime2)
%     disp(log10(CSmeantime2))
%     disp(-log10(T));
    CStotal(iter) = CStotal2;
    CSmeantime(iter) = CSmeantime2;
    CSprobs2(iter) = 1/CSmeantime2;
    CSlogprobs2(iter) = log10(CSmeantime2);
end

fittedx = linspace(0,5,101);
ft = fittype('a*x^n+b');
fitcoeff = fit(-log10(CSy(11:20))',CSlogtable(11:20,26),ft,'Start',[1 2 1]);
fittedy = fitcoeff(fittedx);