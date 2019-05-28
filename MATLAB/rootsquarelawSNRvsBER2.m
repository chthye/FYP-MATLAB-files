p1 = 0.5;
A = 1;
A2 = 1;
B = 0;
DC = 0.5;
numThresh = 100;

mu = 0;
digits(40);

estBER = zeros(101,1);
theoBER = zeros(101,1);
SNRspace = linspace(10,20,101);
for j = 0:100
    threshMin = B+DC;
    threshMax = A+DC; %will be overwritten
    threshStep = (threshMax-threshMin)/numThresh;%will be overwritten
    SNR = 10+0.1*j
    sigma = sqrt(0.5*A^2/(2*(10^((SNR)/10))));
    N = makedist('Normal','mu',mu,'sigma',sigma);
    thresh = linspace(B,A,numThresh+1);

    %for Y = (X+N)^2
    %e = (0.5)*(cdf(N,sqrt(thresh)-A))-0.5*(cdf(N,-sqrt(thresh)-A))+0.5*cdf(N,-sqrt(thresh)-B)+0.5*(cdf(N,sqrt(thresh)-B,'upper'));
    %for Y = (X+N+DC)^2/4
    e = (0.5)*(cdf(N,A2*thresh-A-DC))-0.5*(cdf(N,-A2*thresh-A-DC))+0.5*cdf(N,-A2*thresh-B-DC)+0.5*(cdf(N,A2*thresh-B-DC,'upper'));
    %note: upper returns complement of cdf, much more accurate than 1-cdf
    e = -vpa(log10(e));

    %find max -log10BER and set to centre
%     [M,I] = max(e);
%     threshMax = threshMin+2*threshStep*(I-1);
%     threshStep = (threshMax-threshMin)/numThresh;
%     thresh = linspace(threshMin,threshMax,numThresh+1);


    %for Y = (X+N)^2
    %e = (0.5)*(cdf(N,sqrt(thresh)-A))-0.5*(cdf(N,-sqrt(thresh)-A))+0.5*cdf(N,-sqrt(thresh)-B)+0.5*(cdf(N,sqrt(thresh)-B,'upper'));
    %for Y = (X+N+1)^2/4
%     e = (0.5)*(cdf(N,A2*thresh-A-DC))-0.5*(cdf(N,-A2*thresh-A-DC))+0.5*cdf(N,-A2*thresh-B-DC)+0.5*(cdf(N,A2*thresh-B-DC,'upper'));
%     e = -vpa(log10(e));


    %plot(thresh,e)
    %hold on





    thresh2 = linspace(threshMin,threshMax,numThresh+1);
    numBits = 100000;
    e2 = zeros(1,numThresh+1);

    for i = 0:numThresh
        dataIn = randi([0 1],1,numBits);
        conditionX = (dataIn==0);
        dataIn(conditionX) = A;
        dataIn(~conditionX) = B;
        n = normrnd(mu,sigma,1,numBits);%mu,sigma,numrows,numcolumns
        y = dataIn+n;
        % added SJS
        %yh=(y+DC).^2/A2;
        yh=abs(y+DC)/A2;
        %yh = y.^2;
        y=yh;

        threshold = threshMin+threshStep*i;
        y2 = y;
        conditionY2 = (y2 >= threshold);
        y2(conditionY2) = A;
        y2(~conditionY2) = B;
        z = dataIn-y2;
        z = abs(z);
        e2(i+1) = 0.5*sum(z)/numBits;
    end

    e2 = -log10(e2);
    %plot(thresh2,e2,'o')

    firstinf = find(e2==inf);
    if(isempty(firstinf))
        firstinf = 50;
    end
    e3 = e2(1:firstinf-1);
    e3 = reshape(e3,length(e3),1);

    thresh3 = reshape(linspace(threshMin,threshMin-threshStep+threshStep*length(e3),length(e3)),length(e3),1);
    % ft3 = fittype('a*exp(b*x)+c');
    % option3 = fitoptions(ft3);
    % option3.Lower = [1 3 -inf];
    % option3.Upper = [inf inf  0];
    % f3 = fit(thresh3, e3, ft3,'StartPoint',[2 1 5])
    % plot(f3)
    ft3 = fittype('a*x^2+b*x+c');
    f3 = fit(thresh3,e3,ft3,'StartPoint',[0 0 0]);
    %plot(f3)

    lastinf = find(e2==inf,1,'last');
    if(isempty(lastinf))
        lastinf = 50;
    end
    e4 = e2(lastinf+1:length(e2));
    e4 = reshape(e4,length(e4),1);
    thresh4 = reshape(linspace(threshMax+threshStep-threshStep*length(e4),threshMax,length(e4)),length(e4),1);
    ft4 = fittype('a*x^2+b*x+c');
    f4 = fit(thresh4,e4,ft4,'StartPoint',[0 0 0]);
    %plot(f4)

    %grid
    %xlabel('threshold')
    %ylabel('-log10BER')
    %hold off    

    fun = @(x) f3(x)-f4(x);
    [M,I] = max(e);
    theoxx = thresh(I);
    theoyy = e(I);
    estxx = fzero(fun,1.2);
    estyy = f3(estxx);
    theoBER(j+1) = theoyy;
    estBER(j+1) = estyy;
    %disp(theoxx)
    %disp(theoyy)
    disp(estxx)
    disp(estyy)
end

plot(SNRspace,theoBER,SNRspace,estBER,'x')
grid
xlabel('SNR (dB)')
ylabel('-log10BER')