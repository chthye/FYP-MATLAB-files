numBits = 1000000;
%additive noise
mu = 0;
SNR = 15;
sigma = sqrt(1/(2*(10^((SNR)/10))));
%phase noise
deltanu = 1;
T = 0.0005;
%sampling size for loop
sampleSize = 3;

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

e2 = zeros(100,1);
for i = 0:100

    threshold = -1 + 0.02*i;
    pseudoRxSignal3 = rxSignal3 - threshold;
    
    dataSymbolsOut = qamdemod(pseudoRxSignal3,M);
    dataOutMatrix = de2bi(dataSymbolsOut,k);
    % dataOut = dataOutMatrix(:);
    dataOut = reshape(dataOutMatrix',numBits,1);

    e2(i+1) = biterr(dataIn,dataOut)/numBits;

end
% disp(ber)

% for j = 1:100
%     berj = biterr(dataIn(j*100-99:j*100),dataOut(j*100-99:j*100))/100;
%     disp(berj)
%     
% end

% plot(phi)
% hold on
% plot(estphi)
% hold off
% grid
% legend('phi','estimated phi')
% xlabel('i^t^h symbol')
% ylabel('angle in rad')
thresh2 = linspace(-1,1,101);
% plot (thresh2,e2)

e22 = erfcinv(4*e2);
e222 = -log10(e2);
plot(thresh2,e222,'x')
hold on

firstzero = find(e2==0);
if(isempty(firstzero))
    firstzero = 50;
end
e33 = e22(1:firstzero-1);
e33 = reshape(e33,length(e33),1);

thresh3 = reshape(linspace(-1,-1.02+0.02*length(e33),length(e33)),length(e33),1);
poly33coeff = polyfit(thresh3,e33,1);
sigma33 = abs(1/(sqrt(2)*poly33coeff(1)));
mu33 = -poly33coeff(2)*sqrt(2)*sigma33;

lastzero = find(e2==0,1,'last');
if(isempty(lastzero))
    lastzero = 50;
end
e44 = e22(lastzero+1:length(e22));
e44 = reshape(e44,length(e44),1);

thresh4 = reshape(linspace(1.02-0.02*length(e44),1,length(e44)),length(e44),1);
poly44coeff = polyfit(thresh4,e44,1);
sigma44 = abs(1/(sqrt(2)*poly44coeff(1)));
mu44 = poly44coeff(2)*sqrt(2)*sigma44;

e5 = -log10(0.25*(erfc(abs(mu33-thresh2)/(sigma33*sqrt(2)))+erfc(abs(mu44-thresh2)/(sigma44*sqrt(2)))));
plot(thresh2,e5)

grid
xlabel('threshold')
ylabel('-log10BER')
hold off

