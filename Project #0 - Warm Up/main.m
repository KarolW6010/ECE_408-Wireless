%Karol Wadolowski - Project #0 Warm Up
clear; clc; close all;

%% Part 1: No Channel QAM
%First confirm that the theoretical and simulated BER and SER curves match
%for several QAM constellations.

%Simulation parameters
numIter = 1e2;                          %Number of iterations
numSyms = 1e3;                          %Number of symbols per iteration
SNRdB = 0:2:16;                         %SNRs to use in dB
varN = 10.^(-SNRdB/10);                 %Noise variance
lenSNR = length(SNRdB);                 %Number of SNRs to test
M = [2,4,8,16];                         %QAM constellation size

chan = 1;                               %No channel

%Run the simulation
figure('units','normalized','outerposition',[0,0,1,1])
for ii = 1:length(M)
    bar = waitbar(0,[num2str(M(ii)), ' QAM']);
    
    EbNodB = SNRdB - 10*log10(log2(M(ii)));         %Energy per bit in dB
    %The SNR to EbNo mapping is just dividing the linear SNR by number of 
    %bits per symbol.
    
    %Bit/Symbol error rate for each SNR and each iteration
    berIter = nan(lenSNR,numIter);
    serIter = nan(lenSNR,numIter);
    
    for jj = 1:numIter
        %Use the same message for each SNR
        txBits = randi([0,1],numSyms,log2(M(ii)));  %Generate message bits
        txMSG = bi2de(txBits,'left-msb');           %Transmitted message
        tx = qammod(txMSG,M(ii));                   %Transmitted QAM signal

        for kk = 1:lenSNR    
            txChan = filter(chan,1,tx);                     %Apply the channel
            
            %Note that the built awgn does not add complex noise for BPSK
            %and thus causes issues
            if(M(ii) == 2)
                noise = sqrt(varN(kk)/2)*(randn(numSyms,1)+1j*randn(numSyms,1));
                txNoisy = txChan + noise;                       %Add noise
            else
                txNoisy = awgn(txChan,SNRdB(kk),'measured');    %Add noise
            end
            
            rxMSG = qamdemod(txNoisy,M(ii));                %Received message
            rxBits = de2bi(rxMSG,log2(M(ii)),'left-msb');   %Received bits

            [~,berIter(kk,jj)] = biterr(txBits,rxBits);     %Calculate BER
            serIter(kk,jj) = sum(txMSG ~= rxMSG)/numSyms;   %Calculate SER
        end
        waitbar(jj/numIter,bar)
    end
    
    berSim = mean(berIter,2);       %Average BER for each simulated SNR
    serSim = mean(serIter,2);       %Average SER for each simulated SNR
    
    %Theoretical BER and SER curves
    if(M(ii) == 2)
        %BPSK
        EbNo = 10.^(EbNodB/10);
        berTh = qfunc(sqrt(2*EbNo));
        serTh = berTh;
    else
        %General QAM
        [berTh,serTh] = berawgn(EbNodB,'qam',M(ii));
        %Note: For non rectangular qam configurations the theoretical and
        %simulation BER curves do not match up. I believe this is because
        %either there is a term missing in the SNR to EbNo mapping at the 
        %top of the code or (what I think is more likely) the berawgn 
        %function does not properly account for decision boundaries in the 
        %non-rectangular QAM configurations. So this would be for 32*4*k 
        %QAM sizes where k is a natural number. Note this refers to the QAM
        %constellations used by the qammod function.
    end
    
    subplot(2,length(M),ii)
    semilogy(SNRdB,berTh,'Linewidth',2)
    hold on;
    semilogy(SNRdB,berSim,'--x','Linewidth',2,'MarkerSize',10)
    xlabel('SNR (dB)')
    ylabel('BER')
    title([num2str(M(ii)), ' QAM: Bit Error Rate'])
    legend('Theoretical','Simulation')
    grid on
    
    subplot(2,length(M),length(M)+ii)
    semilogy(SNRdB,serTh,'Linewidth',2)
    hold on;
    semilogy(SNRdB,serSim,'--x','Linewidth',2,'MarkerSize',10)
    xlabel('SNR (dB)')
    ylabel('SER')
    title([num2str(M(ii)), ' QAM: Symbol Error Rate'])
    legend('Theoretical','Simulation')
    grid on
    
    close(bar)
end
sgtitle('BER and SER for Various Size QAM')

%% Part 1: BPSK with moderate ISI channel
%For a moderate ISI channel using BPSK introduce BPSK to bring down the BER
%Note: For BPSK BER and SER are the same, hence one plot

numIter = 5e2;                          %Number of iterations
numSyms = 1e3;                          %Number of symbols per iteration
SNRdB = 0:2:12;                         %SNRs to use in dB
varN = 10.^(-SNRdB/10);                 %Noise variance for BPSK
lenSNR = length(SNRdB);                 %Number of SNRs to test
M = 2;                                  %BPSK
trainLen = 1e2;                         %Number of training symbols
chan = [1, .2, .4];                     %Moderate ISI channel

%Bit error rate for each SNR and iteration
berNone  = nan(lenSNR,numIter);         %No ISI
berModUn = nan(lenSNR,numIter);         %Moderate ISI, unequalized channel
ber2 = nan(lenSNR,numIter);             %Moderate ISI, equalized channel

%Run the simulation
bar = waitbar(0,'Moderate ISI');
for ii = 1:numIter
    %Use the same message for each SNR
    txBits = randi([0,1],numSyms,log2(M));  %Generate message bits
    txMSG = bi2de(txBits,'left-msb');       %Transmitted message
    tx = qammod(txMSG,M);                   %Transmitted QAM signal

    for jj = 1:lenSNR
        %No ISI
        txChan = tx;                                    %Apply the channel
        noise = sqrt(varN(jj)/2)*(randn(numSyms,1)+1j*randn(numSyms,1));
        txNoisy = txChan + noise;                       %Add noise
        rxMSG = qamdemod(txNoisy,M);                    %Received message
        rxBits = de2bi(rxMSG,log2(M),'left-msb');       %Received bits
        
        [~,berNone(jj,ii)] = biterr(txBits,rxBits);     %Calculate BER
        
        %Moderate ISI, unequalized channel
        txChan = filter(chan,1,tx);                     %Apply the channel
        noise = sqrt(varN(jj)/2)*(randn(numSyms,1)+1j*randn(numSyms,1));
        txNoisy = txChan + noise;                       %Add noise
        rxMSG = qamdemod(txNoisy,M);                    %Received message
        rxBits = de2bi(rxMSG,log2(M),'left-msb');       %Received bits
        
        [~,berModUn(jj,ii)] = biterr(txBits,rxBits);    %Calculate BER
        
        %Moderate ISI, equalized channel
        %Use the same txNosiy as in the unequalized channel and use an
        %equalizer. Using the same txNoisy gives a direct performance
        %comparison.
        eq1 = dfe(1,2,lms(1e-2));                       %Equalizer
        eq1.SigConst = qammod(0:(M-1),M);               %Set constellation
        rxEq = equalize(eq1,txNoisy,tx(1:trainLen));    %Received eq signal
        rxMSG = qamdemod(rxEq,M);                       %Received message
        rxBits = de2bi(rxMSG,log2(M),'left-msb');       %Received bits
        
        [~,ber2(jj,ii)] = biterr(txBits(trainLen+1:end,:),...
            rxBits(trainLen+1:end,:));                    %Calculate BER
    end
    waitbar(ii/numIter,bar)
end
close(bar)

%Average the BER curves
berNone  = mean(berNone,2);
berModUn = mean(berModUn,2);
ber2 = mean(ber2,2);

%Theoretical BER for BPSK no ISI
berTh = qfunc(sqrt(2*10.^(SNRdB/10)));              

figure('units','normalized','outerposition',[0,0,1,1])
semilogy(SNRdB,berTh,'Linewidth',2)
hold on
semilogy(SNRdB,berNone,'--x','Linewidth',2,'MarkerSize',10)
semilogy(SNRdB,berModUn,'--x','Linewidth',2,'MarkerSize',10)
semilogy(SNRdB,ber2,'--x','Linewidth',2,'MarkerSize',10)
legend('No ISI Theory','No ISI Simulation','Moderate ISI Un','Moderate ISI Eq')
title('BER Curves for a Moderate ISI Channel')
xlabel('SNR (dB)')
ylabel('BER')
grid on

%% Part 2: BER of 1e-6 at SNR of 12dB in a moderate ISI channel

numIter = 4e3;                          %Number of iterations
numSyms = 1e3;                          %Number of symbols per iteration
SNRdB = 0:2:12;                         %SNRs to use in dB
varN = 10.^(-SNRdB/10);                 %Noise variance for BPSK
lenSNR = length(SNRdB);                 %Number of SNRs to test
k = 1;                                  %Number of info bits
n = 2;                                  %Number of output bits from encoder
M = 2^n;                                %M QAM
trainLen = 3e2;                         %Number of training symbols
fftaps = 1;                             %Number of feed forward taps in eq
fbtaps = 2;                             %Number of feed back raps in eq
trellis = poly2trellis([9],[561,753]);  %Rate 1/2 convolutional encoder
chan = [1, .2, .4];                     %Moderate ISI channel

%Bit error rate for each SNR and iteration
ber2 = nan(lenSNR,numIter);             %BER for part 2

%Run the simulation
bar = waitbar(0,'Moderate ISI');
for ii = 1:numIter
    %Use the same message for each SNR
    txBits = randi([0,1],k,numSyms);        %Generate message bits
    txBits = reshape(txBits,1,[]);
    
    msgBits = convenc(txBits,trellis);      %Encode the data
    msgBits = reshape(msgBits,n,numSyms)';  
    txMSG = bi2de(msgBits,'left-msb');      %Transmitted message
    tx = qammod(txMSG,M);                   %Transmitted QAM signal

    for jj = 1:lenSNR
        %Moderate ISI, unequalized channel
        txChan = filter(chan,1,tx);                     %Apply the channel
        
        if(M == 2)
            noise = sqrt(varN(jj)/2)*(randn(numSyms,1)+1j*randn(numSyms,1));
            txNoisy = txChan + noise;                       %Add noise
        else
            txNoisy = awgn(txChan,SNRdB(jj),'measured');    %Add noise
        end

        eq1 = dfe(fftaps,fbtaps,lms(1e-2));             %Equalizer
        eq1.SigConst = qammod(0:(M-1),M);               %Set constellation
        rxEq = equalize(eq1,txNoisy,tx(1:trainLen));    %Received eq signal
        
        rxMSG = qamdemod(rxEq,M);                       %Received message
        rxBits = de2bi(rxMSG,log2(M),'left-msb');       %Received bits
        rxBits = reshape(rxBits',1,[]);
        decBits = vitdec(rxBits,trellis,numSyms,'trunc','hard');
        
        ind = k*trainLen+1;     %k due to the ECC
        [~,ber2(jj,ii)] = biterr(txBits(ind:end),...
            decBits(ind:end));                          %Calculate BER
    end
    waitbar(ii/numIter,bar)
end
close(bar)

%Average the BER curve
ber2 = mean(ber2,2);

%Theoretical BER for BPSK no ISI
berTh = qfunc(sqrt(2*10.^(SNRdB/10)));              

figure('units','normalized','outerposition',[0,0,1,1])
semilogy(SNRdB,berTh,'Linewidth',2)
hold on
semilogy(SNRdB,ber2,'--x','Linewidth',2,'MarkerSize',10)
legend('No ISI Theory','Moderate ISI: Equalization + ECC')
title('BER Curves for a Moderate ISI Channel')
xlabel('SNR (dB)')
ylabel('BER')
grid on







