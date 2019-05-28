%ver4: changing phase noise detection window

numBits = 1000000;
%additive noise
mu = 0;
SNR = 12;
sigma = sqrt(1/(2*(10^((SNR)/10))));
%phase noise
deltanu = 1;
T = 0.0015;
%sampling size for loop
sampleSize = 9;
if mod(sampleSize,2)==0
    sampleSize = sampleSize+1;
end
phaseWindow = 2.0*pi;

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
    rxSignal3(iteration:iteration+(sampleSize-1)/2) = rxSignal2(iteration:iteration+(sampleSize-1)/2).*exp(-1i*estphi(iteration-1));
    phases = angle(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4);
    phases(phases > phaseWindow/2) = phases(phases>phaseWindow/2) - phaseWindow;
    phases(phases < -phaseWindow/2) = phases(phases < -phaseWindow/2) + phaseWindow;
    
    estdphi(iteration) =  mean(phases)/4;
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

% scatterplot(rxSignal3)
% figure(2);
% plot(phi)
% hold on
% plot(estphi)
% grid
% hold off

% %changing thresholds
% e2 = zeros(100,1);
% for i = 0:100
% 
%     threshold = -pi/4 + 0.005*pi*i;
%     pseudoRxSignal3 = rxSignal3.*exp(1i*threshold);
%     
%     dataSymbolsOut = qamdemod(pseudoRxSignal3,M);
%     dataOutMatrix = de2bi(dataSymbolsOut,k);
%     % dataOut = dataOutMatrix(:);
%     dataOut = reshape(dataOutMatrix',numBits,1);
% 
%     e2(i+1) = biterr(dataIn,dataOut)/numBits;
% 
% end
% 
dataSymbolsOut = qamdemod(rxSignal3,M);
dataOutMatrix = de2bi(dataSymbolsOut,k);
% dataOut = dataOutMatrix(:);
dataOut = reshape(dataOutMatrix',numBits,1);
disp(biterr(dataIn,dataOut)/numBits);

%Cycle slip counter
% disp(sum(abs(diff(round(2*(phi-estphi)/pi)))));
errorphi = 2*(estphi - phi)/pi;
errorphi2 = errorphi;
for iteration = 1:(sampleSize-1)/2
    errorphi2(iteration) = mean(errorphi(1:iteration+(sampleSize-1)/2));
end
for iteration = ((sampleSize+1)/2):(numBits/k-(sampleSize-1)/2)
    errorphi2(iteration) = mean(errorphi(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2));
end
for iteration = (numBits/k-(sampleSize-3)/2):numBits/k
    errorphi2(iteration) = mean(errorphi(iteration-(sampleSize-1)/2:numBits/k));
end
CSindicator = abs(diff(round(errorphi2)));
CStotal = sum(CSindicator);
CSmeantime = numBits/(k*CStotal);
disp(CSmeantime)
% if CStotal > 2
%     CSfirst = find(CSindicator==1,1);
%     CSlast = find(CSindicator==1,1,'last');
%     CSmeantime = (CSlast-CSfirst)/(CStotal-1);
%     disp(CSmeantime)
% end
plot(phi)
hold on
plot(estphi)
hold off
grid
legend('phi','estimated phi')
xlabel('i^t^h symbol')
ylabel('angle in rad')
