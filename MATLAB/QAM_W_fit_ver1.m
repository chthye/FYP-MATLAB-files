%from QPSK_W_fit_ver2
%16QAM
%sliding sample window, ber threshold in phase
%Decision Directed Algorithm

numBits = 100000;
%additive noise
mu = 0;
SNR = 15;
sigma = sqrt(1/(2*(10^((SNR)/10))));
%phase noise
deltanu = 1;
T = 0.00001;
%sampling size for loop
sampleSize = 45;
if mod(sampleSize,2)==0
    sampleSize = sampleSize+1;
end
sigmaPhase = sqrt(2*pi*deltanu*T);
%alpha for wiener filter
wienerParm = (sigmaPhase^2+2*sigma^2-sigmaPhase*sqrt(sigmaPhase^2+4*sigma^2))/(2*sigma^2);
wienerFiltCoeff = zeros(sampleSize,1);
%calculate coeff as (1-a)*a^k for both sides
for i = 1:sampleSize
   wienerFiltCoeff(i) = (1-wienerParm)*wienerParm^(abs(i-(sampleSize+1)/2)); 
end
wienerFiltCoeff = wienerFiltCoeff/sum(wienerFiltCoeff);

dataIn = randi([0 1], numBits, 1);
%16QAM
M = 16;
k = log2(M);
qamAngle = angle(qammod([0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15], M));%for easy access to angle values
tempSymbols = zeros(numBits/k,1);%temporary store symbols for calculating estdphi
dataInMatrix = reshape(dataIn, k, length(dataIn)/k)';
dataSymbolsIn = bi2de(dataInMatrix);

dataMod = qammod(dataSymbolsIn, M);

txSignal = dataMod;
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
% pnoise = normrnd(mu,  ,length(txSignal),1);
% rxSignal3 = rxSignal2.*exp(1i*pnoise);
rxSignalAngle2 = angle(rxSignal2);


estdphi = zeros(numBits/k,1);
estphi = zeros(numBits/k,1);


rxSignal3(1:sampleSize) = rxSignal2(1:sampleSize);
%starting iteration
for iteration = 1:(sampleSize-1)/2
    rxSignal3(1:iteration+(sampleSize-1)/2) = rxSignal2(1:iteration+(sampleSize-1)/2).*exp(-1i*estphi(max(1,iteration-1)));
    
    tempSymbols(iteration:iteration+(sampleSize-1)/2) = qamdemod(rxSignal3(iteration:iteration+(sampleSize-1)/2),M);
    estdphi(iteration) = mean(angle(rxSignal3(1:iteration+(sampleSize-1)/2))-qamAngle(tempSymbols(1:iteration+(sampleSize-1)/2)+1)');
%    estdphi(iteration) =  mean(angle(-1*rxSignal3(1:iteration+(sampleSize-1)/2).^4))/4;
   if iteration>1
       estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
   else
       estphi(iteration) = estdphi(iteration);
   end
end
%main iteration
for iteration = ((sampleSize+1)/2):(numBits/k-(sampleSize-1)/2)
%     rxSignal3(iteration:iteration+(sampleSize-1)/2) = rxSignal2(iteration:iteration+(sampleSize-1)/2).*exp(-1i*estphi(iteration-1));
    rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2) = rxSignal2(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).*exp(-1i*estphi(iteration-1));
    
    tempSymbols(iteration:iteration+(sampleSize-1)/2) = qamdemod(rxSignal3(iteration:iteration+(sampleSize-1)/2),M);
    estdphi(iteration) = sum((angle(rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2))-qamAngle(tempSymbols(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2)+1)').*wienerFiltCoeff);

%     estdphi(iteration) = angle((sum(-1*rxSignal3(iteration-(sampleSize-1)/2:iteration+(sampleSize-1)/2).^4.*wienerFiltCoeff))^(1/4));
    estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
end
%ending iteration
for iteration = (numBits/k-(sampleSize-3)/2):numBits/k
%     rxSignal3(iteration:numBits/k) = rxSignal2(iteration:numBits/k).*exp(-1i*estphi(iteration-1));
    rxSignal3(iteration-(sampleSize-1)/2:numBits/k) = rxSignal2(iteration-(sampleSize-1)/2:numBits/k).*exp(-1i*estphi(iteration-1));
    
    tempSymbols(iteration:numBits/k) = qamdemod(rxSignal3(iteration:numBits/k),M);
    estdphi(iteration) = mean(angle(rxSignal3(iteration-(sampleSize-1)/2:numBits/k))-qamAngle(tempSymbols(iteration-(sampleSize-1)/2:numBits/k)+1)');
%     estdphi(iteration) =  mean(angle(-1*rxSignal3(iteration-(sampleSize-1)/2:numBits/k).^4))/4;
    estphi(iteration) = estphi(iteration-1) + estdphi(iteration);
end
rxSignal3 = rxSignal2.*exp(-1i*estphi);

errorphi = estphi-phi;
% scatterplot(rxSignal3)
% figure(2);
plot(phi)
hold on
plot(estphi)
grid
hold off

%changing thresholds
e2 = zeros(100,1);
for i = 50%0:100

    threshold = -pi/4 + 0.005*pi*i;
    pseudoRxSignal3 = rxSignal3.*exp(1i*threshold);
    
    dataSymbolsOut = qamdemod(pseudoRxSignal3,M);
    dataOutMatrix = de2bi(dataSymbolsOut,k);
    % dataOut = dataOutMatrix(:);
    dataOut = reshape(dataOutMatrix',numBits,1);

    e2(i+1) = biterr(dataIn,dataOut)/numBits;

end
disp(e2(51))
% 
% % plot(phi)
% % hold on
% % plot(estphi)
% % hold off
% % grid
% % legend('phi','estimated phi')
% % xlabel('i^t^h symbol')
% % ylabel('angle in rad')
% thresh2 = linspace(-1,1,101);
% % plot (thresh2,e2)
% 
% e22 = erfcinv(4*e2);
% e222 = -log10(e2);
% 
% figure(3);
% plot(thresh2,e222,'x')
% hold on
% 
% firstzero = find(e2==0);
% if(isempty(firstzero))
%     firstzero = 50;
% end
% e33 = e22(1:firstzero-1);
% e33 = reshape(e33,length(e33),1);
% 
% thresh3 = reshape(linspace(-1,-1.02+0.02*length(e33),length(e33)),length(e33),1);
% poly33coeff = polyfit(thresh3,e33,1);
% sigma33 = abs(1/(sqrt(2)*poly33coeff(1)));
% mu33 = -poly33coeff(2)*sqrt(2)*sigma33;
% 
% lastzero = find(e2==0,1,'last');
% if(isempty(lastzero))
%     lastzero = 50;
% end
% e44 = e22(lastzero+1:length(e22));
% e44 = reshape(e44,length(e44),1);
% 
% thresh4 = reshape(linspace(1.02-0.02*length(e44),1,length(e44)),length(e44),1);
% poly44coeff = polyfit(thresh4,e44,1);
% sigma44 = abs(1/(sqrt(2)*poly44coeff(1)));
% mu44 = poly44coeff(2)*sqrt(2)*sigma44;
% 
% e5 = -log10(0.25*(erfc(abs(mu33-thresh2)/(sigma33*sqrt(2)))+erfc(abs(mu44-thresh2)/(sigma44*sqrt(2)))));
% plot(thresh2,e5)
% 
% grid
% xlabel('threshold')
% ylabel('-log10BER')
% hold off
% 
% % disp(max(e5))
% 
% %numerical integration using measured variance
% e = zeros(101,1);
% 
% epsilon = linspace(-pi,pi-pi/200,400);
% alpha = 1/var(errorphi);
% for i = 1:101
%     threshold = -1.02+0.02*i;
%     %phase error pdf (Tikhonov dist) (shifted by threshold before multiplying with BER given phase error)
% %     ptheta = exp(alpha*cos(linspace(-pi-threshold*pi/4,pi-pi/200-threshold*pi/4,400)))/(2*pi*besseli(0,alpha));
%     %Gaussian approx
%     ptheta = exp(-linspace(-pi-threshold*pi/4,pi-pi/200-threshold*pi/4,400).^2/(2*var(errorphi)))/(std(errorphi)*sqrt(2));
%     ptheta = ptheta/sum(ptheta);%normalise probability
%     for eps = 51:150
%         e(i) = e(i) + (0.5-0.25*erfc(sin(-epsilon(eps)-pi/4)/(sqrt(2)*sigma))+0.25*erfc(sin(3*pi/4+epsilon(eps))/(sqrt(2)*sigma))).*ptheta(eps);
%     end
%     for eps = 151:200
%         e(i) = e(i) + (0.25*erfc(sin(pi/4+epsilon(eps))/(sqrt(2)*sigma))+0.25*erfc(sin(pi/4-epsilon(eps))/(sqrt(2)*sigma))).*ptheta(eps);
%     end
%     for eps = 201:250
%         e(i) = e(i) + (0.25*erfc(sin(pi/4-epsilon(eps))/(sqrt(2)*sigma))+0.25*erfc(sin(pi/4+epsilon(eps))/(sqrt(2)*sigma))).*ptheta(eps);
%     end
%     for eps = 251:350
%         e(i) = e(i) + (0.5-0.25*erfc(sin(epsilon(eps)-pi/4)/(sqrt(2)*sigma))+0.25*erfc(sin(3*pi/4-epsilon(eps))/(sqrt(2)*sigma))).*ptheta(eps);
%     end
%     
%     e(i) = -log10(e(i));
% end
% 
% eNoPhaseError = 0.25*erfc(sin(pi/4-pi*thresh2/4)/(sqrt(2)*sigma))+0.25*erfc(sin(pi/4+pi*thresh2/4)/(sqrt(2)*sigma));
% eNoPhaseError = -log10(eNoPhaseError);
% hold on
% plot(thresh2,e,'--')
% plot(thresh2,eNoPhaseError,'.')
% hold off
% %for checking scatterplot distribution variance
% rxSignal4 = rxSignal3(real(rxSignal3)>0 & imag(rxSignal3)>0);
% rxSignal4 = rxSignal4 - exp(1i*pi/4);
% rxSignal4 = rxSignal4*exp(-1i*pi/4);
% % sum(abs(real(rxSignal4)).^2)
% % sum(abs(imag(rxSignal4)).^2)