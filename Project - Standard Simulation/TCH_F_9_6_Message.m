%Karol Wadolowski, Project: GSM PHY Layer Simulation
%This code uses TCH /F9.6 to send a message over various SNR values
clear; clc; close all;
warning('off')  %For a clearer output. Turned on at the end

SNRdB = -6:1:-2;        %Signal to noise ratios to test

%Prepare message
message = ['\tType your message here!\n\tNote your message will have ',...
    'zeros appended\n\tto it so that it has a length divisible by 240.'];
fprintf(['Original Message:\n',message,'\n\n'])
data = de2bi(double(message),8)';
data = data(:);

%Increase bit length to be a multiple of 240
extra = mod(length(data),240);
if(extra ~= 0)
    data = [data;zeros(240-extra,1)];
end

data = reshape(data,240,[]);
N = size(data,2);

%Create the 1/2 rate convolutional encoder and decoder with puncturing
trellis = poly2trellis([5],[23,33]);
convEncoder = comm.ConvolutionalEncoder(trellis);
convEncoder.PuncturePatternSource = 'Property';
convEncoder.PuncturePattern = 1*(mod((0:487)'-11,15) > 0);

vitDecoder = comm.ViterbiDecoder(trellis,'InputFormat','Hard');
vitDecoder.PuncturePatternSource = 'Property';
vitDecoder.PuncturePattern = convEncoder.PuncturePattern;
vitDecoder.TracebackDepth = 224;
vitDecoder.TerminationMethod = 'Terminated';

%Differential Decoder
diffDecoder = comm.DifferentialEncoder();

%Create modulator and demodulator
sps = 4;        %Samples per symbol
pLen = 4;       %GMSK filter length
gmskMod = comm.GMSKModulator('BitInput',true,'PulseLength',pLen,...
    'SamplesPerSymbol',sps,'InitialPhaseOffset',pi/4);
gmskDem = comm.GMSKDemodulator('BitOutput',true,'PulseLength',pLen,...
    'SamplesPerSymbol',sps,'InitialPhaseOffset',pi/4);

%Interleaving indices
k = (0:455)';
bitLoc = mod(k,19) + 19*mod(k,6) + 1;   %Plus one for indexing (j)
const = mod(k,19) + floor(k/114) + 1;   %Constant portion of (B) with Bo = 1

%Time slot formatting constants
ts = 4*N+18;                %Number of time slows (bursts) needed
tailBits = zeros(3,ts);     %Tail bits
midBits = zeros(28,ts);     %Midamble bits and stealing flags (Fig 11.10 Rappaport)
guardbits = ones(8,ts);     %Bits in between time slots


%The following code sends the message in various SNR channels and displays
%the received message. This code is almost identical to the BER version.

info = [data;zeros(4,N)];           %Append 4 zero bits to clear convolutional coder

%Encode the data (240 bits at a time aka 1 block)
enc = nan(456,N);                   
for jj = 1:N
    enc(:,jj) = convEncoder(info(:,jj));  
end

%Perform interleaving
inter = zeros(114,ts);              %Not nan because still needs to be modulated
for jj = 1:N
    burstLoc = 4*(jj-1) + const;
    for kk = 1:456
        inter(bitLoc(kk),burstLoc(kk)) = enc(kk,jj);
    end
end

%Now each column of 114 bits can be formatted in a time slot
formatted = vertcat(tailBits,inter(1:57,:),midBits,inter(58:end,:),...
    tailBits,guardbits);
form1 = [formatted(:);zeros(gmskDem.TracebackDepth,1)];
%^Zeros are added to the end to deal with the demodulating delay.

%Differentially encode
diffEnc = xor(form1,[0;form1(1:(end-1))]);

%Modulate the data
sent = gmskMod(diffEnc);

%Send the modulated data in different SNR channels
for jj = 1:length(SNRdB)
    %Create noise channel
    awgnChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
        'SNR',SNRdB(jj));

    %Apply noise to sent signal and then perform GMSK demodulation
    noisy = awgnChan(sent);                         %Add noise
    rec1 = gmskDem(noisy);                          %Demodulate
    rec = rec1((gmskDem.TracebackDepth+1):end);     
    %^Remove the zeros added to account for demodulating delay

    %Differentially decode
    diffDec = diffDecoder(rec);
    
    %Undo formatting. Keep only encoded data
    unform1 = reshape(diffDec,156,[]);
    unformatted = vertcat(unform1(4:60,:),unform1(89:145,:));

    %Inverse Interleaving
    Ndetec = (size(unformatted,2)-18)/4;    %Amount of blocks detected (same as N)
    deInt = nan(456,Ndetec);
    for kk = 1:Ndetec
        burstLoc = 4*(kk-1) + const;
        for mm = 1:456
            deInt(mm,kk) = unformatted(bitLoc(mm),burstLoc(mm));
        end
    end

    %Decode
    dec = nan(244,Ndetec);
    for kk = 1:Ndetec
        dec(:,kk) = vitDecoder(deInt(:,kk));
    end
    dec1 = dec(1:240,:);    %Remove the 4 added zeros

    retMess = reshape(dec1,8,[]);
    retMess = char(bi2de(retMess')');

    fprintf('SNR = %d dB\n',SNRdB(jj))
    fprintf(retMess)
    fprintf('\n\n')

    %Reset the following for next iteration
    gmskMod.reset();            %Modulator
    gmskDem.reset();            %Demodulator
    convEncoder.reset();        %Encoder
    vitDecoder.reset();         %Decoder
    diffDecoder.reset();        %Differential Decoder
end
    
warning('on')