%LDPC encode, add bit errors before correction according to SNR, then
%decode and calculate corrected BER
SNRstart = 3.0;
SNRstep = 0.03;
SNRtimes = 11;
%DVB-S2 rates: 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 5/6, 8/9, and 9/10.
p = dvbs2ldpc(9/10);
numLoops = 300;
BER = zeros(SNRtimes,1);
for iter = 1:SNRtimes
    SNR = SNRstart + (iter-1)*SNRstep;
    disp(SNR)
%     sigma = sqrt(1/(2*(10^((SNR)/10))));
%     e = 0.25*(erfc(1/(sigma*sqrt(2)))+erfc(abs(1)/(sigma*sqrt(2))));


    numBits = (size(p,2)-size(p,1))*numLoops;
    % totalErrors = 0;
    hEnc = comm.LDPCEncoder(p);
    hMod = comm.PSKModulator(2, 'BitInput',true);
    hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',SNR);
    hDemod = comm.PSKDemodulator(2, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio','Variance', 1/10^(hChan.SNR/10));
    hDec = comm.LDPCDecoder(p);
    hError = comm.ErrorRate;
    for counter = 1:numLoops
        disp(counter)
        data = logical(randi([0 1], (size(p,2)-size(p,1)), 1));
        encodedData    = step(hEnc, data);
        modSignal      = step(hMod, encodedData);
        receivedSignal = step(hChan, modSignal);
        demodSignal    = step(hDemod, receivedSignal);
        receivedBits   = step(hDec, demodSignal);
        errorStats     = step(hError, data, receivedBits);%(1)=error rate,(2)=num errors,(3)=numbits


    end
    % BER = totalErrors/numBits;
    % disp(BER)
    disp(errorStats(1))
    disp(errorStats(2))
    BER(iter) = log10(errorStats(1));
    disp(BER(iter))
end