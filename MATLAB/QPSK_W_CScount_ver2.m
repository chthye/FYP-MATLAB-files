%From QPSK_W_CScount_ver1.m
%use BER to identify CS

numBits = 1000000;
loopTimes = 6;
%additive noise
mu = 0;
SNR = 11;
sigma = sqrt(1/(2*(10^((SNR)/10))));
%phase noise
deltanu = 1;
T = 0.0015;
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


CStotal = zeros(loopTimes,1);
CSmeantime = zeros(loopTimes,1);
for iter = 1:loopTimes
    disp(iter)
    phaseWindow = (1+(iter-1)/(loopTimes-1))*pi;

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

%     scatterplot(rxSignal3)
    % figure(2);
    % plot(phi)
    % hold on
    % plot(estphi)
    % grid
    % hold off

    
    dataSymbolsOut = qamdemod(rxSignal3,M);
    dataOutMatrix = de2bi(dataSymbolsOut,k);
    % dataOut = dataOutMatrix(:);
    dataOut = reshape(dataOutMatrix',numBits,1);
    disp(biterr(dataIn,dataOut)/numBits);

    %Cycle slip counter (ver1)
    % errorphi = 2*(estphi - phi)/pi;
    % errorphi2 = errorphi;
    % for iteration = 1:(sampleSize-1)/2
    %     errorphi2(iteration) = mean(errorphi(1:iteration+(sampleSize-1)/2));
    % end
    % for iteration = ((sampleSize+1)/2):(numBits/k-(sampleSize-1)/2)
    %     errorphi2(iteration) = mean(errorphi(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2));
    % end
    % for iteration = (numBits/k-(sampleSize-3)/2):numBits/k
    %     errorphi2(iteration) = mean(errorphi(iteration-(sampleSize-1)/2:numBits/k));
    % end
    % CSindicator = abs(diff(round(errorphi2)));
    % CStotal = sum(CSindicator);
    % CSmeantime = numBits/(k*CStotal);
    % disp(CStotal)
    % disp(CSmeantime)

    %Cycle slip counter ver 2 starts here
    dataSymbolsOut = qamdemod(rxSignal3,M);
    dataOutMatrix = de2bi(dataSymbolsOut,k);
    % dataOut = dataOutMatrix(:);
    dataOut = reshape(dataOutMatrix',numBits,1);

    bitError = abs(dataIn-dataOut);
    bitErrorSum = cumsum(bitError);
    bitError2 = reshape(bitError,20,numBits/20);
    bitErrorSum2 = sum(bitError2)';
    CSindicator2 = 0;
    CStotal2 = 0;
    for iteration = 1:length(bitErrorSum2)
        switch CSindicator2
            case 0
                if bitErrorSum2(iteration)>8
                    CStotal2 = CStotal2+1;
                    CSindicator2 = 1;
    %                 disp(iteration)
                end
            case 1
                if bitErrorSum2(iteration)>18
                    CStotal2 = CStotal2+1;
                    CSindicator2 = 2;
    %                 disp(iteration)
                elseif bitErrorSum2(iteration)<2
                    CStotal2 = CStotal2+1;
                    CSindicator2 = 0;
    %                 disp(iteration)
                end
            case 2
                if bitErrorSum2(iteration) < 12
                    CStotal2 = CStotal2+1;
                    CSindicator2 = 1;
    %                 disp(iteration)
                end
        end
    end
    CSmeantime2 = numBits/(k*CStotal2);
%     disp(CStotal2)
%     disp(CSmeantime2)
    CStotal(iter) = CStotal2;
    CSmeantime(iter) = CSmeantime2;

end
% figure(1)
% plot(phi)
% hold on
% plot(estphi)
% hold off
% grid
% legend('phi','estimated phi')
% xlabel('i^t^h symbol')
% ylabel('angle in rad')
phaseSpace = linspace(1,2,loopTimes)';

figure(1)
plot(linspace(1,2,loopTimes),CStotal)

figure(2)
plot(linspace(1,2,loopTimes),-log10(CStotal/(numBits/k)))

figure(3)
plot(linspace(1,2,loopTimes),CSmeantime)


% figure(2)
% plot(bitErrorSum)