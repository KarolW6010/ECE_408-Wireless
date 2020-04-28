function [M,Rst,DRs,mmseBERs,zfBERs] = bestOFDMconfig(channel,N,SNRdB,...
    targBER,prnt)
%{
Finds the highest data rate OFDM system configuration using either a MMSE
or Zero-Forcing equalizer.

Inputs:
    channel     : (D, ) Channel that the signal will go through
    N           : (1, ) Generate N OFDM symbols
    SNRdB       : (1, ) SNR of the channel
    targBER     : (1, ) Target BER value
    prnt        : (1, ) Print all configurations (true or false)

Outputs:
    M           : (4, ) The 4 possible QAMs
    Rst         : (3, ) Strings of the 3 possible rates
    DRs         : (4,3) Possible data rates for each combo of M and R
    mmseBERs    : (4,3) MMSE BER for each combo of M and R
    zfBERs      : (4,3) Analogous to above
    
    if(prnt)
    Lists all configurations that can reach the target BER at the given
    SNR. Also displays the best configuration.
%}

SPS = 2.5e5;                %Symbols per second in 802.11a

D = length(channel);        %Channel length
M = [2,4,16,64];            %M-ary QAM types
R = [1/2,2/3,3/4];          %ECC rates
Rst = {'1/2','2/3','3/4'};

DRs = nan(4,3);
mmseBERs = nan(4,3);
zfBERs = nan(4,3);

berMMSEbest = 1;
berZFbest = 1;
DRbestMMSE = 0;
DRbestZF = 0;
mmseBest = false;
zfBest = false;

if(prnt)
fprintf('Target BER = %.0e, SNR = %.1f dB, %.1e Symbols used\n',...
    targBER,SNRdB,N)
end
for ii = 1:length(M)
    for jj = 1:length(R)
        %Generate OFDM symbol
        [DataBits,sym] = OFDMsymGen(N,R(jj),M(ii));
        
        %Apply the channel
        txChan = filter(channel,[1],sym);
        %txChan = txChan(D:end);
        
        %Add noise
        sigPower = mean(abs(txChan).^2);
        noiseVar = sigPower*10^(-SNRdB/10);
        noise = sqrt(noiseVar/2)*...
            (randn(length(txChan),1) + 1j*randn(length(txChan),1));
        rx = txChan + noise;
        
        %Perform equalization
        symMMSE = filter([1],channel+noiseVar.*eye(size(channel)),rx);
        symZF = filter([1],channel,rx);
        
        %Decode OFDM Symbol
        DecBitsMMSE = OFDMdemod(symMMSE,R(jj),M(ii));
        DecBitsZF = OFDMdemod(symZF,R(jj),M(ii));
        
        %Calculate BER
        berMMSE = mean(DecBitsMMSE ~= DataBits);
        berZF = mean(DecBitsZF ~= DataBits);
        
        %Calculate Data Rate in Mbps
        DR = 48*R(jj)*log2(M(ii))*SPS/1e6;      
        
        DRs(ii,jj) = DR;
        mmseBERs(ii,jj) = berMMSE;
        zfBERs(ii,jj) = berZF;
        
        %Check if the target is met
        if(berMMSE < targBER)
            if(DR > DRbestMMSE)
                statMMSE = 'Target Met';
                berMMSEbest = berMMSE;
                MbestMMSE = M(ii);
                RbestMMSE = Rst{jj};
                DRbestMMSE = DR;
                mmseBEST = true;
            end
        else
            if((berMMSE < berMMSEbest) && ~mmseBest)
                statMMSE = 'Target Not Met';
                berMMSEbest = berMMSE;
                MbestMMSE = M(ii);
                RbestMMSE = Rst{jj};
                DRbestMMSE = DR;
            end
        end
        if(berZF < targBER)
            if(DR > DRbestZF)
                statZF = 'Target Met';
                berZFbest = berZF;
                MbestZF = M(ii);
                RbestZF = Rst{jj};
                DRbestZF = DR;
                zfBest = true;
            end
        else            
            if((berZF < berZFbest) && ~zfBest)
                statZF = 'Target Not Met';
                berZFbest = berZF;
                MbestZF = M(ii);
                RbestZF = Rst{jj};
                DRbestZF = DR;
            end
        end
        
        %Print Information
        if(prnt)
            fprintf(['%2d QAM, Rate %s, MMSE BER = %.2e, ZF BER = %.2e',...
            ', Data Rate = %2.0f Mbps\n'],M(ii),Rst{jj},berMMSE,berZF,DR)
        end
        
    end
end

if(prnt)
    fprintf('\n')
    fprintf('Best Configurations\n')
    fprintf(['MMSE:\t %2d QAM, Rate %s, BER = %.2e, Data Rate = %2.0f Mbps',...
        ', %s\n'],MbestMMSE,RbestMMSE,berMMSEbest,DRbestMMSE,statMMSE)
    fprintf(['ZF:\t\t %2d QAM, Rate %s, BER = %.2e, Data Rate = %2.0f Mbps',...
        ', %s\n\n'],MbestZF,RbestZF,berZFbest,DRbestZF,statZF)
end



