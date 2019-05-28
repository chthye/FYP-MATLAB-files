
numBits = 10000;
%additive noise
mu = 0;
SNR = 10;
sigma = sqrt(1/(2*(10^((SNR)/10))));
%phase noise
deltanu = 1;
T = 0.0011;
%sampling size for loop
sampleSize = 10;

dataIn = randi([0 1], numBits, 1);
%4QAM = QPSK
M = 4;
k = log2(M);
dataInMatrix = reshape(dataIn, k, length(dataIn)/k)';
dataSymbolsIn = bi2de(dataInMatrix);

dataMod = qammod(dataSymbolsIn, M);

txSignal = dataMod;

noise = normrnd(mu,sigma,length(txSignal),1);
noise2 = normrnd(mu,sigma,length(txSignal),1);

%Signal before phase noise
rxSignal = txSignal + noise+1i*noise2;

%normalise
% rxSignal = rxSignal./abs(rxSignal);

% rxSignalAngle = angle(rxSignal);
% phaseNoise = normrnd(mu,sigma,length(txSignal),1);
% rxSignalAngle = rxSignalAngle + phaseNoise;
% rxSignal = cos(rxSignalAngle) + 1i*sin(rxSignalAngle);

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

%1st iteration
rxSignal3(1:sampleSize) = rxSignal2(1:sampleSize);

estdphi(1:sampleSize) = mean(angle(-1*rxSignal3(1:sampleSize).^4))/4;
estphi(1:sampleSize) = estdphi(1);
%loop
for iteration = (sampleSize+1):(sampleSize):(numBits/k-sampleSize+1)
    rxSignal3(iteration:iteration+sampleSize-1) = rxSignal2(iteration:iteration+sampleSize-1).*exp(-1i*estphi(iteration-1));
    
    estdphi(iteration:iteration+sampleSize-1) = mean(angle(-1*rxSignal3(iteration:iteration+sampleSize-1).^4))/4;
    estphi(iteration:iteration+sampleSize-1) = estphi(iteration-1)+estdphi(iteration);
    
end




dataSymbolsOut = qamdemod(rxSignal3,M);
dataOutMatrix = de2bi(dataSymbolsOut,k);
% dataOut = dataOutMatrix(:);
dataOut = reshape(dataOutMatrix',numBits,1);
ber500 = biterr(dataIn(1:500),dataOut(1:500))/500;
ber5000 = biterr(dataIn(1:5000),dataOut(1:5000))/5000;
ber = biterr(dataIn,dataOut)/numBits;
disp(ber500)
disp(ber5000)
disp(ber)

% for j = 1:100
%     berj = biterr(dataIn(j*100-99:j*100),dataOut(j*100-99:j*100))/100;
%     disp(berj)
%     
% end

plot(phi)
hold on
plot(estphi)
hold off
grid
legend('phi','estimated phi')
xlabel('i^t^h symbol')
ylabel('angle in rad')