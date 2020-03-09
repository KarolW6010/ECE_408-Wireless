%Karol Wadolowski - Project #1: Alamouti Codes
clear; clc; close all;

%Simulation parameters
iters = 5e2;        %Number of iterations to perform 
syms  = 5e4;        %Number of symbols per iteration 
SNRdB1_1 = 0:2:50;  %SNRs for (1 Tx, 1 Rx) MRRC (BPSK)
SNRdB1_2 = 0:2:30;  %SNRs for (1 Tx, 2 Rx) MRRC
SNRdB1_4 = 0:2:16;  %SNRs for (1 Tx, 4 Rx) MRRC
SNRdB1_6 = 0:2:10;  %SNRs for (1 Tx, 6 Rx) MRRC
SNRdB1_8 = 0:2:6;  %SNRs for (1 Tx, 8 Rx) MRRC

%Get MMRC BER curves
mrrcBER1_1 = berRayleighMRRC(1,SNRdB1_1,iters,syms);
mrrcBER1_2 = berRayleighMRRC(2,SNRdB1_2,iters,syms);
mrrcBER1_4 = berRayleighMRRC(4,SNRdB1_4,iters,syms);
mrrcBER1_6 = berRayleighMRRC(6,SNRdB1_6,iters,syms);
mrrcBER1_8 = berRayleighMRRC(8,SNRdB1_8,iters,syms);

figure
semilogy(SNRdB1_1,mrrcBER1_1,'-o','MarkerFaceColor',ones(3,1));
hold on;
semilogy(SNRdB1_2,mrrcBER1_2,'-v','MarkerFaceColor',ones(3,1));
semilogy(SNRdB1_4,mrrcBER1_4,'-s','MarkerFaceColor',ones(3,1));
semilogy(SNRdB1_6,mrrcBER1_6,'-d','MarkerFaceColor',ones(3,1));
semilogy(SNRdB1_8,mrrcBER1_8,'-^','MarkerFaceColor',ones(3,1));
legend('(1 Tx, 1 Rx) No Diversity','(1 Tx, 2 Rx) MRRC',...
    '(1 Tx, 4 Rx) MRRC','(1 Tx, 6 Rx) MRRC','(1 Tx, 8 Rx) MRRC')
xlabel('SNR (dB)')
ylabel('BER')
grid on
ylim([1e-6,1])
xlim([0,50])

%%
%Simulation parameters
iters = 5e2;        %Number of iterations to perform 
syms  = 5e4;        %Number of symbols per iteration 
SNRdB1_1 = 0:2:50;  %SNRs for (1 Tx, 1 Rx) MRRC (BPSK)
SNRdB2_1 = 0:2:32;  %SNRs for (2 Tx, 1 Rx) Alamouti
SNRdB2_2 = 0:2:16;  %SNRs for (2 Tx, 2 Rx) Alamouti
SNRdB2_3 = 0:2:16;  %SNRs for (2 Tx, 3 Rx) Alamouti
SNRdB2_4 = 0:2:10;  %SNRs for (2 Tx, 4 Rx) Alamouti


%Get MMRC BER curves
mrrcBER1_1 = berRayleighMRRC(1,SNRdB1_1,iters,syms);

%Get Alamouti BER curves
alamBER2_1 = berRayleighAlamouti(1,SNRdB2_1,iters,syms,'sum1');
alamBER2_2 = berRayleighAlamouti(2,SNRdB2_2,iters,syms,'sum1');
alamBER2_3 = berRayleighAlamouti(3,SNRdB2_3,iters,syms,'sum1');
alamBER2_4 = berRayleighAlamouti(4,SNRdB2_4,iters,syms,'sum1');

figure
semilogy(SNRdB1_1,mrrcBER1_1,'-o','MarkerFaceColor',ones(3,1));
hold on;
semilogy(SNRdB2_1,alamBER2_1,'-v','MarkerFaceColor',ones(3,1));
semilogy(SNRdB2_2,alamBER2_2,'-s','MarkerFaceColor',ones(3,1));
semilogy(SNRdB2_3,alamBER2_3,'-d','MarkerFaceColor',ones(3,1));
semilogy(SNRdB2_4,alamBER2_4,'-^','MarkerFaceColor',ones(3,1));
legend('(1 Tx, 1 Rx) No Diversity','(2 Tx, 1 Rx) Alamouti',...
    '(2 Tx, 2 Rx) Alamouti','(2 Tx, 3 Rx) Alamouti','(2 Tx, 4 Rx) Alamouti')
xlabel('SNR (dB)')
ylabel('BER')
grid on
ylim([1e-6,1])