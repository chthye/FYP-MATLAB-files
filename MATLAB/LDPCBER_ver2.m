%from LDPCBER_ver1.m
%use thresholds
modValue = 0.92387953251128684950543856757577+0.38268343236508967075693021797633*1i;
thresh2 = linspace(-1,1,11);

SNR = -2.2;
%DVB-S2 rates: 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 5/6, 8/9, and 9/10.
p = dvbs2ldpc(1/2);
numLoops = 30;
numBits = (size(p,2)-size(p,1))*numLoops;
% totalErrors = 0;
hEnc = comm.LDPCEncoder(p);
hMod = comm.PSKModulator(2, 'BitInput',true);
hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',SNR);
hDemod = comm.PSKDemodulator(2, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio','Variance', 1/10^(hChan.SNR/10));
hDec = comm.LDPCDecoder(p);
hError = comm.ErrorRate;


e2 = zeros(length(thresh2),1);
for iter = 1:length(thresh2)
    disp(iter)
    bias = modValue*thresh2(iter);
%     sigma = sqrt(1/(2*(10^((SNR)/10))));
%     e = 0.25*(erfc(1/(sigma*sqrt(2)))+erfc(abs(1)/(sigma*sqrt(2))));


    
    for counter = 1:numLoops
        disp(counter)
        data = logical(randi([0 1], (size(p,2)-size(p,1)), 1));
        encodedData    = step(hEnc, data);
        modSignal      = step(hMod, encodedData);
        receivedSignal = step(hChan, modSignal);
        %add bias threshold here
        receivedSignal = receivedSignal+bias;
        demodSignal    = step(hDemod, receivedSignal);
        receivedBits   = step(hDec, demodSignal);
        errorStats     = step(hError, data, receivedBits);%(1)=error rate,(2)=num errors,(3)=numbits


    end
    % BER = totalErrors/numBits;
    % disp(BER)
    disp(errorStats(1))
    disp(errorStats(2))
    e2(iter) = log10(errorStats(1));
    disp(e2(iter))
end