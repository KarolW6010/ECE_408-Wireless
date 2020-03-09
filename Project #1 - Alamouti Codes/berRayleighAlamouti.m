function ber = berRayleighAlamouti(Rx,SNRdB,iters,syms,power)
%{
Returns the BER curve for the Alamouti Coding Scheme in a rayleigh fading 
channel for an uncoded coherent BPSK signal. It is assumed that the channel
is perfectly known by the receiver.

Inputs:
    Rx          : (1, ) Number of receiver antennas (integer >= 1)
    SNRdB       : (N, ) N signal to noise ratios in dB 
    iters       : (1, ) Run this many iterations
    syms        : (1, ) Send this many symbols per iteration
    power       : (1, ) Power level at each transmitter ('sum1','sum2')
                  'sum1' corresponds to sending equal power as MRRC
                  'sum2' corresponds to sending double the power of MRRC

Outputs:
    ber         : (N, ) BER perfromance curve obtained from the various
                  SNRs
%}

noiseVar = 10.^(-SNRdB/10);         %Convert SNRdB to noise variance
N = length(SNRdB);                  %Number of SNRs to test

berInd = nan(iters,N);              %Will hold individual BER values

bar = waitbar(0,'Iterations Progress');
for ii = 1:iters
    %Generate the bits to be sent
    bits = randi([0,1],syms,1);
    tx = 1 - 2*bits;            %Sent signal
    
    %Separate into 2 Tx streams (Alamouti Scheme)
    temp = reshape(tx,2,[]);    %Top row is s0, Bottom row is s1
    temp0 = temp;           
    temp0(2,:) = -temp(2,:);    %Conjugate and negate second row
    temp1 = flipud(temp);       %Conjugate second row and flip
    %Note conjugation is not actually performed as we are doing BPSK(-1,+1)
    tx0 = temp0(:);             %Transmitter 0 signal
    tx1 = temp1(:);             %Transmitter 1 signal
    
    %Generate Rx rayleigh fading channels
    rayChan1 = cell(2,1);       %Will hold the channel values
    rayChan2 = cell(2,1);       %Same as above just values appear twice
    rayChanConj = cell(2,1);    %Conjuagte of rayleigh fading channel 
    for jj = 1:2
        rayChan1{jj} = nan(syms/2,Rx);
        for kk = 1:Rx
            rayChan1{jj}(:,kk) = rayleighFading(syms/2);
        end
        
        rayChan2{jj} = repelem(rayChan1{jj},2,1);
        rayChanConj{jj} = conj(rayChan1{jj});
    end
    
    %Apply the channel to the sent signal
    txChan0 = rayChan2{1}.*repmat(tx0,1,Rx);
    txChan1 = rayChan2{2}.*repmat(tx1,1,Rx);
    
    %Get the components of each received signal (after applying channel)
    rx0_0 = downsample(txChan0,2,0);        %h_(2k)    s0
    rx0_1 = downsample(txChan0,2,1);        %h_(2k)  (-s1*)
    rx1_0 = downsample(txChan1,2,0);        %h_(2k+1)  s1
    rx1_1 = downsample(txChan1,2,1);        %h_(2k+1)  s0*
    
    rxEven = rx0_0 + rx1_0;
    rxOdd  = rx0_1 + rx1_1;
    
    %Measure signal power for each channel so that noise is properly scaled
    switch power
        case{'sum1'}
            sigPowerEven = mean(abs(rxEven).^2,1);
            sigPowerOdd  = mean(abs(rxOdd).^2,1);
        case{'sum2'}
            sigPowerEven = mean(abs(rxEven).^2,1)/2;
            sigPowerOdd  = mean(abs(rxOdd).^2,1)/2;
        otherwise
            error('Choose valid power option!');
    end
    
    %For each SNR
    for jj = 1:N 
        %Add AWGN noise to signal
        noiseEven = sqrt(noiseVar(jj)/2)*(randn(syms/2,Rx) + 1j*randn(syms/2,Rx));
        noiseEven = repmat(sqrt(sigPowerEven),syms/2,1).*noiseEven;
        rxEvenNoisy = rxEven + noiseEven;       %Received signal at time t
        
        noiseOdd = sqrt(noiseVar(jj)/2)*(randn(syms/2,Rx) + 1j*randn(syms/2,Rx));
        noiseOdd = repmat(sqrt(sigPowerOdd),syms/2,1).*noiseOdd;
        rxOddNoisy = rxOdd + noiseOdd;          %Received signal at time t+T
        
        %Apply channel estimator and employ ML detector
        sHat0 = sum(rayChanConj{1}.*rxEvenNoisy + rayChan1{2}.*conj(rxOddNoisy),2);
        sHat1 = sum(rayChanConj{2}.*rxEvenNoisy - rayChan1{1}.*conj(rxOddNoisy),2);
        
        %Format
        sHatTemp = [sHat0.';sHat1.'];
        sHat = sHatTemp(:);
        
        %Use estimate to guess sent symbol. For BPSK this comes down to
        %taking the sign of the real part.
        dec = real(sHat) < 0;
        
        %Sum up bit errors
        berInd(ii,jj) = sum(dec ~= bits);     
    end
    
    waitbar(ii/iters,bar)
end
berInd = berInd/syms;       %Divide by number of symbols per iteration
ber = mean(berInd,1);       %Take the mean BER curve

close(bar)
