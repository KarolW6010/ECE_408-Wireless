function ber = berRayleighMRRC(Rx,SNRdB,iters,syms)
%{
Returns the BER curve for the Maximal-Ratio Receive Combing Scheme in a
rayleigh fading channel for an uncoded coherent BPSK signal. It is assumed 
that the channel is perfectly known by the receiver.

Inputs:
    Rx          : (1, ) Number of receiver antennas (integer >= 1)
    SNRdB       : (N, ) N signal to noise ratios in dB 
    iters       : (1, ) Run this many iterations
    syms        : (1, ) Send this many symbols per iteration

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
    
    %Map to BPSK constellation
    tx = 1 - 2*bits;                %Sent signal
    
    %Generate Rx rayleigh fading channels
    rayChan = nan(syms,Rx);         %Will hold the channel values
    for jj = 1:Rx
        rayChan(:,jj) = rayleighFading(syms);
    end
    
    %Conjugate of rayleigh fading channel for receiving
    rayChanConj = conj(rayChan); 
    
    %Apply the channel to the sent signal
    txChan = rayChan.*repmat(tx,1,Rx);
    
    %Measure signal power for each channel so that noise is properly scaled
    sigPower = mean(abs(txChan).^2,1);
    
    %For each SNR
    for jj = 1:N
        %Add AWGN noise to signal
        noise = sqrt(noiseVar(jj)/2)*(randn(syms,Rx) + 1j*randn(syms,Rx));
        noise = repmat(sqrt(sigPower),syms,1).*noise;
        rx = txChan + noise;        %Received signal
        
        %Multiply by channel conjugate and sum to get signal estimate
        sHat = sum(rayChanConj.*rx,2);
        
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









