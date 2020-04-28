%Karol Wadolowski - Project #3: MIMO (Part 1)
clear; clc; close all;

N = 3e5;                %Number of symbols to generate
SNRdB = 35;             %SNR in dB
M = [2,4,16,64];        %M-ary QAM
R = [1/2,2/3,3/4];      %ECC code rates
Rst = {'1/2','2/3','3/4'};
fD = 10.^(2:4);         %Max Doppler Frequency for Rayleigh Fading
targBER = 10.^(-(2:4)); %Target System BER

%Will hold the BERs for each configuration
BER_SVD  = cell(length(fD),1);
BER_MMSE = cell(length(fD),1);
BER_ZF   = cell(length(fD),1);

%Will hold the possible data rates
DR = nan(length(M),length(R));     

%Calculate the BERs for all the scenarios
for ii = 1:length(fD)
    fprintf('Calculating BERs for fD = %.e Hz....', fD(ii))
    
    BER_SVD{ii}  = nan(length(M),length(R));
    BER_MMSE{ii} = nan(length(M),length(R));
    BER_ZF{ii}   = nan(length(M),length(R));
    
    for jj = 1:length(M)
        for kk = 1:length(R)
            %Generate Data
            [DataBits,data,H,dr] = genMIMOdata(2,M(jj),R(kk),N,fD(ii));
            DR(jj,kk) = dr;
            
            %Send and Decode data various ways
            berSVD = precodingMIMO(2,M(jj),R(kk),SNRdB,data,H,DataBits);
            [berMMSE,berZF] = mmse_zfMIMO(2,M(jj),R(kk),SNRdB,data,H,...
                DataBits);
            
            %Store BER
            BER_SVD{ii}(jj,kk)  = berSVD;
            BER_MMSE{ii}(jj,kk) = berMMSE;
            BER_ZF{ii}(jj,kk)   = berZF;
        end
    end
    fprintf('Complete!\n')
end
fprintf('\n')

DR = DR/1e6;    %Data rate in Mbps

%For each target BER and channel type print the best configurations. If
%some configuration met the BER target then choose the highest Data Rate.
%If no configuration met the BER choose the best BER configuration.
for ii = 1:length(fD)
    fprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',...
        '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n'])
    
    fprintf(['Flat fading Rayleigh channel with max Doppler frequency ',...
        '%.0e Hz\n'],fD(ii))
    for jj = 1:length(targBER)
        %For SVD
        if(min(BER_SVD{ii},[],'all') < targBER(jj))
            mask = BER_SVD{ii} < targBER(jj);
            temp = mask.*DR;
            [bestSVD,svdR] = max(temp,[],1);
            [~,svdC] = max(bestSVD);
            statSVD = 'Target Met';
        else
            [bestSVD,svdR] = min(BER_SVD{ii},[],1);
            [~,svdC] = min(bestSVD);
            statSVD = 'Target Not Met';
        end
           
        %For MMSE
        if(min(BER_MMSE{ii},[],'all') < targBER(jj))
            mask = BER_MMSE{ii} < targBER(jj);
            temp = mask.*DR;
            [bestMMSE,mmseR] = max(temp,[],1);
            [~,mmseC] = max(bestMMSE);
            statMMSE = 'Target Met';
        else
            [bestMMSE,mmseR] = min(BER_MMSE{ii},[],1);
            [~,mmseC] = min(bestMMSE);
            statMMSE = 'Target Not Met';
        end
        
        %For ZF
        if(min(BER_ZF{ii},[],'all') < targBER(jj))
            mask = BER_ZF{ii} < targBER(jj);
            temp = mask.*DR;
            [bestZF,zfR] = max(temp,[],1);
            [~,zfC] = max(bestZF);
            statZF = 'Target Met';
        else
            [bestZF,zfR] = min(BER_ZF{ii},[],1);
            [~,zfC] = min(bestZF);
            statZF = 'Target Not Met';
        end
        
        fprintf('Target BER = %.0e, SNR = %.1f dB, %.1e Symbols used\n',...
            targBER(jj),SNRdB,N)
        fprintf('Best Configurations\n')
        
        fprintf(['%4s:\t %2d QAM, Rate %s, BER = %.2e, Data Rate = ',...
            '%3.0f Mbps, %s\n'],'SVD',M(svdR(svdC)),Rst{svdC},...
            BER_SVD{ii}(svdR(svdC),svdC),DR(svdR(svdC),svdC),statSVD)
        
        fprintf(['%4s:\t %2d QAM, Rate %s, BER = %.2e, Data Rate = ',...
            '%3.0f Mbps, %s\n'],'MMSE',M(mmseR(mmseC)),Rst{mmseC},...
            BER_MMSE{ii}(mmseR(mmseC),mmseC),DR(mmseR(mmseC),mmseC),...
            statMMSE)
        
        fprintf(['%4s:\t %2d QAM, Rate %s, BER = %.2e, Data Rate = ',...
            '%3.0f Mbps, %s\n\n'],'ZF',M(zfR(zfC)),Rst{zfC},...
            BER_ZF{ii}(zfR(zfC),zfC),DR(zfR(zfC),zfC),statZF)
    end
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    







