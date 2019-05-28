%ber threshold in phase
numBits = 5000000;
%additive noise
mu = 0;
SNR = 13;
sigma = sqrt(1/(2*(10^((SNR)/10))));
%phase noise
deltanu = 1;
T = 0.0003;
%sampling size for loop
sampleSizes = [5;7;9;11;13;15;17;19;21;23;25;27;29;31;33];

sigmaPhase = sqrt(2*pi*deltanu*T);
%alpha for wiener filter
wienerParm = (sigmaPhase^2+2*sigma^2-sigmaPhase*sqrt(sigmaPhase^2+4*sigma^2))/(2*sigma^2);

variances = zeros(length(sampleSizes),1);
BER = zeros(length(sampleSizes),1);
estBER = zeros(length(sampleSizes),1);

for iter=1:length(sampleSizes)
    sampleSize = sampleSizes(iter);
    disp(sampleSize)
    
    wienerFiltCoeff = zeros(sampleSize,1);
    %calculate coeff as (1-a)*a^k for both sides
    for i = 1:sampleSize
       wienerFiltCoeff(i) = (1-wienerParm)*wienerParm^(abs(i-(sampleSize+1)/2)); 
    end
    wienerFiltCoeff = wienerFiltCoeff/sum(wienerFiltCoeff);

    
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
%         rxSignal3(iteration:iteration+(sampleSize-1)/2) = rxSignal2(iteration:iteration+(sampleSize-1)/2).*exp(-1i*estphi(iteration-1));
        rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2) = rxSignal2(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).*exp(-1i*estphi(iteration-1));%ver 9/2
%         estdphi(iteration) =  mean(angle(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4))/4;
%         estdphi(iteration) = sum(angle(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4).*wienerFiltCoeff)/4;
        estdphi(iteration) = angle((sum(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4.*wienerFiltCoeff))^(1/4));

        estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
    end
    %ending iteration
    for iteration = (numBits/k-(sampleSize-3)/2):numBits/k
        rxSignal3(iteration-(sampleSize-1)/2:numBits/k) = rxSignal2(iteration-(sampleSize-1)/2:numBits/k).*exp(-1i*estphi(iteration-1));
        estdphi(iteration) =  mean(angle(-1*rxSignal3(iteration-(sampleSize-1)/2:numBits/k).^4))/4;
        estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
    end
    rxSignal3 = rxSignal2.*exp(-1i*estphi);

    errorphi = sort(estphi-phi);
    variances(iter) = var(estphi-phi);
    % scatterplot(rxSignal3)
    % figure(2);
    % plot(phi)
    % hold on
    % plot(estphi)
    % grid
    % hold off

    %changing thresholds
    e2 = zeros(100,1);
    for i = 0:100

        threshold = -pi/4 + 0.005*pi*i;
        pseudoRxSignal3 = rxSignal3.*exp(1i*threshold);

        dataSymbolsOut = qamdemod(pseudoRxSignal3,M);
        dataOutMatrix = de2bi(dataSymbolsOut,k);
        % dataOut = dataOutMatrix(:);
        dataOut = reshape(dataOutMatrix',numBits,1);

        e2(i+1) = biterr(dataIn,dataOut)/numBits;

    end


    % plot(phi)
    % hold on
    % plot(estphi)
    % hold off
    % grid
    % legend('phi','estimated phi')
    % xlabel('i^t^h symbol')
    % ylabel('angle in rad')
    thresh2 = linspace(-1,1,101);
    % plot (thresh2,e2)

    e22 = erfcinv(4*e2);
    e222 = -log10(e2);

%     figure(3);
%     plot(thresh2,e222,'x')
%     hold on

    firstzero = find(e2==0);
    if(isempty(firstzero))
        firstzero = 50;
    end
    e33 = e22(1:firstzero-1);
    e33 = reshape(e33,length(e33),1);

    thresh3 = reshape(linspace(-1,-1.02+0.02*length(e33),length(e33)),length(e33),1);
    poly33coeff = polyfit(thresh3,e33,1);
    sigma33 = abs(1/(sqrt(2)*poly33coeff(1)));
    mu33 = -poly33coeff(2)*sqrt(2)*sigma33;

    lastzero = find(e2==0,1,'last');
    if(isempty(lastzero))
        lastzero = 50;
    end
    e44 = e22(lastzero+1:length(e22));
    e44 = reshape(e44,length(e44),1);

    thresh4 = reshape(linspace(1.02-0.02*length(e44),1,length(e44)),length(e44),1);
    poly44coeff = polyfit(thresh4,e44,1);
    sigma44 = abs(1/(sqrt(2)*poly44coeff(1)));
    mu44 = poly44coeff(2)*sqrt(2)*sigma44;

    e5 = -log10(0.25*(erfc(abs(mu33-thresh2)/(sigma33*sqrt(2)))+erfc(abs(mu44-thresh2)/(sigma44*sqrt(2)))));
%     plot(thresh2,e5)
% 
%     grid
%     xlabel('threshold')
%     ylabel('-log10BER')
%     hold off
    
    BER(iter) = e222(51);
    estBER(iter)=(max(e5));
end
figure(1)
plot(sampleSizes,variances)
figure(2)
plot(sampleSizes,BER,sampleSizes,estBER,'--')