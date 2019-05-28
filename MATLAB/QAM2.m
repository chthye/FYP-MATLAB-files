%code from https://uk.mathworks.com/help/comm/gs/compute-ber-for-a-qam-system-with-awgn-using-matlab.html

%bits for transmission
numBits = 5000;
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
EbNo = 20;
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
sigma = 0.1;
phaseNoise = normrnd(mu,sigma,length(rxFiltSignalAngle),1);
rxFiltSignalAngle = rxFiltSignalAngle + phaseNoise;

clear('j')
rxFiltSignal2 = rxFiltSignalMag.*cos(rxFiltSignalAngle) + 1i*rxFiltSignalMag.*sin(rxFiltSignalAngle);

scatterplot(rxFiltSignal2)


dataSymbolsOut = qamdemod(rxFiltSignal2, M);
dataOutMatrix = de2bi(dataSymbolsOut,k);
dataOut = dataOutMatrix(:);
ber = biterr(dataIn,dataOut)/numBits
