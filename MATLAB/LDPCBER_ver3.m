%from LDPCBER_ver2.m
%use same transmitted and received signal when changing thresholds

modValue = 0.92387953251128684950543856757577+0.38268343236508967075693021797633*1i;
thresh2 = linspace(-0.25,0.25,21);

SNR = 1.575;
%DVB-S2 rates: 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 5/6, 8/9, and 9/10.
p = dvbs2ldpc(4/5);
numLoops = 80;
numBits = (size(p,2)-size(p,1))*numLoops;
% totalErrors = 0;
hEnc = comm.LDPCEncoder(p);
hMod = comm.PSKModulator(2, 'BitInput',true);
hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',SNR);
hDemod = comm.PSKDemodulator(2, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio','Variance', 1/10^(hChan.SNR/10));
hDec = comm.LDPCDecoder(p);
hError = comm.ErrorRate;
dataTable = logical(randi([0 1], (size(p,2)-size(p,1)), numLoops));
receivedSignalTable = zeros(size(p,2), numLoops);
e2 = zeros(length(thresh2),1);

for counter = 1:numLoops
    disp(counter)
    data = dataTable(:,counter);
    encodedData    = step(hEnc, data);
    modSignal      = step(hMod, encodedData);
    receivedSignal = step(hChan, modSignal);
    receivedSignalTable(:,counter) = receivedSignal;
end
for iter = 1:length(thresh2)
    disp(iter)
    bias = modValue*thresh2(iter);
%     sigma = sqrt(1/(2*(10^((SNR)/10))));
%     e = 0.25*(erfc(1/(sigma*sqrt(2)))+erfc(abs(1)/(sigma*sqrt(2))));


    errorSum = 0;
    for counter = 1:numLoops
        disp(counter)
        data = dataTable(:,counter);
        %add bias threshold here
        receivedSignal = receivedSignalTable(:,counter)+bias;
        demodSignal    = step(hDemod, receivedSignal);
        receivedBits   = step(hDec, demodSignal);
%         errorStats     = step(hError, data, receivedBits);%(1)=error rate,(2)=num errors,(3)=numbits
        errorSum = errorSum+biterr(data,receivedBits);

    end
    % BER = totalErrors/numBits;
    % disp(BER)
    disp(errorSum)
%     disp(errorStats(1))
%     disp(errorStats(2))
    e2(iter) = log10(errorSum/(length(data)*numLoops));
    disp(e2(iter))
end