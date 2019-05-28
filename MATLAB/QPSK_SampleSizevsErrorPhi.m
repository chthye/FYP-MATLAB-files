numBits = 1000000;
%additive noise
mu = 0;
SNR = 14;
sigma = sqrt(1/(2*(10^((SNR)/10))));
%phase noise
deltanu = 1;
T = 0.0003;
%sampling size for loop
sampleSizes = 11;%[7;9;11;13;15;17;19;21;23;25];

sigmaPhase = sqrt(2*pi*deltanu*T);
%alpha for wiener filter
wienerParm = (sigmaPhase^2+2*sigma^2-sigmaPhase*sqrt(sigmaPhase^2+4*sigma^2))/(2*sigma^2);


variances = zeros(length(sampleSizes),1);
spread = zeros(length(sampleSizes),1);
spread2 = zeros(length(sampleSizes),1);
for iter = 1:length(sampleSizes)
    sampleSize = sampleSizes(iter);
    disp(sampleSize)
% if mod(sampleSize,2)==0
%     sampleSize = sampleSize+1;
% end

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
        estdphi(iteration) =  mean(angle(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4))/4;
%         estdphi(iteration) = sum(angle(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4).*wienerFiltCoeff)/4;
%         estdphi(iteration) = angle((sum(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4.*wienerFiltCoeff))^(1/4));

        estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
    end
    %ending iteration
    for iteration = (numBits/k-(sampleSize-3)/2):numBits/k
        rxSignal3(iteration-(sampleSize-1)/2:numBits/k) = rxSignal2(iteration-(sampleSize-1)/2:numBits/k).*exp(-1i*estphi(iteration-1));
        estdphi(iteration) =  mean(angle(-1*rxSignal3(iteration-(sampleSize-1)/2:numBits/k).^4))/4;
        estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
    end
    rxSignal3 = rxSignal2.*exp(-1i*estphi);
    
    rxSignal4 = rxSignal3(real(rxSignal3)>0 & imag(rxSignal3)>0);
    rxSignal4 = rxSignal4 - exp(1i*pi/4);
    rxSignal4 = rxSignal4*exp(-1i*pi/4);
    % sum(abs(real(rxSignal4)).^2)
    % sum(abs(imag(rxSignal4)).^2)
    spread(iter) = sum(abs(real(rxSignal4)))-sum(abs(imag(rxSignal4)));
    spread2(iter) = sum(abs(real(rxSignal4)).^2)-sum(abs(imag(rxSignal4)).^2);
    errorphi = sort(estphi-phi);
%     plot(errorphi(1:10:500000),linspace(0,1,50000))
%     hold on
    variances(iter) = var(estphi-phi);
    % scatterplot(rxSignal3)
    % figure(2);
    % plot(phi)
    % hold on
    % plot(estphi)
    % grid
    % hold off
end
figure(4)
plot(sampleSizes,variances)
figure(5)
plot(sampleSizes,spread,sampleSizes,spread2)
grid