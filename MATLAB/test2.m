%code from https://uk.mathworks.com/help/comm/gs/compute-ber-for-a-qam-system-with-awgn-using-matlab.html
%changed demod part for different threshold



%bits for transmission
numBits = 20000;
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

%root raised cosine filter
rolloff = 0.5;
numSamplesPerSymbol = 8;        
span = 4;
rrcFilter = rcosdesign(rolloff,span,numSamplesPerSymbol);
%fvtool(rrcFilter,'Analysis','Impulse') %to visualise filter

%upsample and apply root raised cosine filter
txSignal = upfirdn(dataMod, rrcFilter, numSamplesPerSymbol, 1);
figure(1);
plot(linspace(1,length(txSignal),length(txSignal)),txSignal)
hold on

%set snr(AWGN)
%EbNo = 10;
%snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);
snr = 10;
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

%own decoder for bpsk from here
ber = zeros(101,1);
for n = 0:100
    rxDecodedSignal = rxFiltSignal;
    threshold = -1+0.02*n;
    decodeCondition = (rxDecodedSignal >= threshold);
    rxDecodedSignal(decodeCondition) = 1;
    rxDecodedSignal(~decodeCondition) = -1;
    %demodulate
    demodCondition = (rxDecodedSignal == -1);
    rxDecodedSignal(demodCondition) = 0;
    ber(n+1) = biterr(dataIn, rxDecodedSignal)/numBits;
end
figure(2);
ber = -log10(ber);
plot(linspace(-1,1,101),ber)
hold off

