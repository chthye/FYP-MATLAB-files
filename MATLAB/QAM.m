%code from https://uk.mathworks.com/help/comm/gs/compute-ber-for-a-qam-system-with-awgn-using-matlab.html

%bits for transmission
numBits = 10000;
dataIn = randi([0 1], numBits, 1);
%plot(linspace(1,numBits,numBits),dataIn) %plot bits out for test

%group bits to integers for QAM. M=16 means regroup 4 bits, 16QAM
M = 16;
k = log2(M);
dataInMatrix = reshape(dataIn, length(dataIn)/k, k);
dataSymbolsIn = bi2de(dataInMatrix);
%plot(linspace(1,length(dataSymbolsIn),length(dataSymbolsIn)),dataSymbolsIn)

%QAM modulation
dataMod = qammod(dataSymbolsIn, M);
%plot(linspace(1,length(dataMod),length(dataMod)),dataMod)
% scatterplot(dataMod)

%root raised cosine filter
rolloff = 0.5;
numSamplesPerSymbol = 8;        
span = 4;
rrcFilter = rcosdesign(rolloff,span,numSamplesPerSymbol);
%fvtool(rrcFilter,'Analysis','Impulse') %to visualise filter

%upsample and apply root raised cosine filter
txSignal = upfirdn(dataMod, rrcFilter, numSamplesPerSymbol, 1);
% figure(1);
% plot(linspace(1,length(txSignal),length(txSignal)),txSignal)
% hold on

%set snr(AWGN)
EbNo = 40;
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);
%note: EbNo,snr in dB. 
%snr = Psignal/Pnoise = (Esignal*symbolrate)/(noisePSD*bandwidth)
%Es = Eb*log2(M),noisePSD = No??,  ???bandwidth=fs/2???check simulation

%pass through channel, to receiver
rxSignal = awgn(txSignal, snr, 'measured');
% plot(linspace(1,length(rxSignal),length(rxSignal)),rxSignal)
% hold off


% Downsample and filter
rxFiltSignal = upfirdn(rxSignal,rrcFilter,1,numSamplesPerSymbol);
% Account for delay
rxFiltSignal = rxFiltSignal(span+1:end-span);
rxFiltSignalMag = abs(rxFiltSignal);
rxFiltSignalAngle = angle(rxFiltSignal);

%phase change code here
mu = 0;
sigma = 0.05;
phaseNoise = normrnd(mu,sigma,length(rxFiltSignalAngle),1);
rxFiltSignalAngle = rxFiltSignalAngle + phaseNoise;
for i = 1:length(rxFiltSignalAngle)
    rxFiltSignalAngle(i) = rxFiltSignalAngle(i) + 0.0003*i;
end



clear('j')
rxFiltSignal2 = rxFiltSignalMag.*cos(rxFiltSignalAngle) + 1i*rxFiltSignalMag.*sin(rxFiltSignalAngle);

scatterplot(rxFiltSignal2)

%PLL and demodulate here
% dataSymbolsOut = qamdemod(rxFiltSignal2, M);

qamAngle = angle(qammod([0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15], M));
dataSymbolsOut = zeros(length(rxFiltSignal2),1);
rxFiltSignal3 = rxFiltSignal2;

angleDiff = zeros(length(rxFiltSignal2),1);
angleChange = zeros(length(rxFiltSignal2),1);
dataSymbolsOut(1) = qamdemod(rxFiltSignal2(1), M);
angleDiff(1) = angle(rxFiltSignal3(1)) - qamAngle(dataSymbolsOut(1)+1);
angleChange(1) = 0 + angleDiff(1);
actualDiff = 0.0000;
angleChange(1) = 0 + actualDiff;
for i = 2:length(rxFiltSignal3)
    rxFiltSignal3(i) = rxFiltSignal3(i)*(cos(angleChange(i-1))-1i*sin(angleChange(i-1)));
    dataSymbolsOut(i) = qamdemod(rxFiltSignal3(i), M);
    angleDiff(i) = angle(rxFiltSignal3(i)) - qamAngle(dataSymbolsOut(i)+1);
    angleChange(i) = angleChange(i-1) + angleDiff(i);
%     angleChange(i) = angleChange(i-1) + actualDiff;
end
scatterplot(rxFiltSignal3)
dataOutMatrix = de2bi(dataSymbolsOut,k);
dataOut = dataOutMatrix(:);
ber = biterr(dataIn,dataOut)/numBits
