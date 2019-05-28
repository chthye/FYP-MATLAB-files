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
numSamplesPerSymbol = 10;        
span = 4;
rrcFilter = rcosdesign(rolloff,span,numSamplesPerSymbol);
%fvtool(rrcFilter,'Analysis','Impulse') %to visualise filter

%upsample and apply root raised cosine filter
txSignal = upfirdn(dataMod, rrcFilter, numSamplesPerSymbol, 1);
figure(1);
plot(linspace(1,length(txSignal),length(txSignal)),txSignal)
hold on

%set snr(AWGN)
EbNo = 10;
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);
%note: EbNo,snr in dB. 
%snr = Psignal/Pnoise = (Esignal*symbolrate)/(noisePSD*bandwidth)
%Es = Eb*log2(M),noisePSD = No??,  ???bandwidth=fs/2???check simulation

%pass through channel, to receiver
rxSignal = awgn(txSignal, snr, 'measured');
plot(linspace(1,length(rxSignal),length(rxSignal)),rxSignal)
hold off



% Downsample and filter
rxFiltSignal = upfirdn(rxSignal,rrcFilter,1,numSamplesPerSymbol);
% Account for delay
rxFiltSignal = rxFiltSignal(span+1:end-span);

scatterplot(rxFiltSignal)

%demod. threshold is changed by shifting filtered signal
ber = zeros(101,101);
for x = 0:100
    for y = 0:100
        rxFiltSignalTemp = rxFiltSignal-1-1i+0.02*(x+y*1i);
        dataSymbolsOut = qamdemod(rxFiltSignalTemp, M);
        dataOutMatrix = de2bi(dataSymbolsOut,k);
        dataOut = dataOutMatrix(:);
        ber(x+1,y+1) = biterr(dataIn,dataOut)/numBits;
       
    end
end
ber2 = -log10(ber);
figure(3);
hold on
for ii = 26:76
plot(linspace(-1,1,101),ber2(:,ii))
end

hold off
