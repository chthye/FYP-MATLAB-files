%from QPSK_W_CScount_ver3.m
%for different SNR, same linewidth, loop to get CS prob
CSprobs = zeros(1,26);
CSlogprobs = zeros(1,26);
for iter = 26
    disp(iter)
    numBits = 100000;
    %additive noise
    mu = 0;
    SNR = 10+0.2*(iter-1);
    sigma = sqrt(1/(2*(10^((SNR)/10))));
    %phase noise
    deltanu = 1;
    T = 0.0025;
    %sampling size for loop
    sampleSize = 45;
    if mod(sampleSize,2)==0
        sampleSize = sampleSize+1;
    end
    phaseWindow = 2*pi;
    CSrequired = 200;
%     CSdist = zeros(1,CSrequired);
    
    sigmaPhase = sqrt(2*pi*deltanu*T);
    %alpha for wiener filter
    wienerParm = (sigmaPhase^2+2*sigma^2-sigmaPhase*sqrt(sigmaPhase^2+4*sigma^2))/(2*sigma^2);
    wienerFiltCoeff = zeros(sampleSize,1);
    %calculate coeff as (1-a)*a^k for both sides
    for i = 1:sampleSize
       wienerFiltCoeff(i) = (1-wienerParm)*wienerParm^(abs(i-(sampleSize+1)/2)); 
    end
    wienerFiltCoeff = wienerFiltCoeff/sum(wienerFiltCoeff);

    loopNumber = 0;
    CStotal = 0;
    CStimes = zeros(CSrequired,1);
    while CStotal < CSrequired

        loopNumber = loopNumber + 1;
    %     disp(loopNumber)
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

        %     estdphi(iteration) =  mean(phases)/4;
        %     estdphi(iteration) = sum(phases.*wienerFiltCoeff)/4;
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

        % scatterplot(rxSignal3)
        % figure(2);
        % plot(phi)
        % hold on
        % plot(estphi)
        % grid
        % hold off



        %Cycle slip counter ver 2 starts here
        dataSymbolsOut = qamdemod(rxSignal3,M);
        dataOutMatrix = de2bi(dataSymbolsOut,k);
        % dataOut = dataOutMatrix(:);
        dataOut = reshape(dataOutMatrix',numBits,1);

        bitError = abs(dataIn-dataOut);
    %     disp(sum(bitError)/numBits)
    %     bitErrorSum = cumsum(bitError);
        if sum(bitError) > 100
            bitError2 = reshape(bitError,20,numBits/20);
            bitErrorSum2 = sum(bitError2)';
            CSindicator2 = 0;

            for iteration = 1:length(bitErrorSum2)
                if CSindicator2 == 0
                    if bitErrorSum2(iteration)>8
                        CStotal = CStotal+1;
                        disp(CStotal)
    %                     CStime = CStime + 20*iteration;
                        CStimes(CStotal) = ((loopNumber-1)*numBits+20*iteration)/k-sampleSize;
                        disp(CStimes(CStotal))
                        if 20*iteration/k < sampleSize
                            CStotal = CStotal-1;
                            loopNumber = loopNumber-1;
                        else
                            loopNumber = 0;
                        end
                        
                        CSindicator2 = 1;
    %                     disp(CStotal)
                    end
                end
            end
        end
    end
    % CStime = (CStime + (numBits)*(loopNumber-CSrequired))/k;
    CSmeantime = mean(CStimes);

    disp(CSmeantime)

    CSprobs(iter) = 1/CSmeantime;
    CSlogprobs(iter) = log10(CSmeantime);
end
