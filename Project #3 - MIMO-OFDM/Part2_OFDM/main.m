%Karol Wadolowski - Project #3: OFDM (Part 2)
clear; clc; close all;

N = 4e4;                            %Number of OFDM symbols to generate
SNRdB = 10;                 
targBER = 10.^(-(2:4));             %Target System BER
channel = cell(1,1);                %Channels
channel{1} = 1;                     %No ISI
channel{2} = [1,-.3];               %More ISI
channel{3} = [1,.2,.4];             %Even More ISI
channel{4} = [1,.2,.4,-.3,.1,-.7];  %Tons of ISI

%For each BER and channel print the best configuration. If some
%configuration met the BER target then choose the highest Data Rate. If no 
%configuration met the target BER choose the best BER configuration.
for ii = 1:length(channel)
    current_Channel = channel{ii}
    [M,R,DRs,mmseBERs,zfBERs] = bestOFDMconfig(channel{ii},N,SNRdB,...
        1,false);
    for jj = 1:length(targBER)
        if(min(mmseBERs,[],'all') < targBER(jj))
            mask = mmseBERs < targBER(jj);
            temp = mask.*DRs;
            [mmseMax,mmseR] = max(temp,[],1);
            [~,mmseC] = max(mmseMax);
            statMMSE = 'Target Met';
        else
            [mmseMax,mmseR] = min(mmseBERs,[],1);
            [~,mmseC] = min(mmseMax);
            statMMSE = 'Target Not Met';
        end
        
        if(min(zfBERs,[],'all') < targBER(jj))
            mask = zfBERs < targBER(jj);
            temp = mask.*DRs;
            [zfMax,zfR] = max(temp,[],1);
            [~,zfC] = max(zfMax);
            statZF = 'Target Met';
        else
            [zfMax,zfR] = min(zfBERs,[],1);
            [~,zfC] = min(zfMax);
            statZF = 'Target Not Met';
        end
        
        fprintf('Target BER = %.0e, SNR = %.1f dB, %.1e Symbols used\n',...
        targBER(jj),SNRdB,N)
        fprintf('Best Configurations\n')
        fprintf(['%4s:\t %2d QAM, Rate %s, BER = %.2e, Data Rate = ',...
            '%2.0f Mbps, %s\n'],'MMSE',M(mmseR(mmseC)),R{mmseC},...
            mmseBERs(mmseR(mmseC),mmseC),DRs(mmseR(mmseC),mmseC),statMMSE)
        fprintf(['%4s:\t %2d QAM, Rate %s, BER = %.2e, Data Rate = ',...
            '%2.0f Mbps, %s\n\n'],'ZF',M(zfR(zfC)),R{zfC},...
            zfBERs(zfR(zfC),zfC),DRs(zfR(zfC),zfC),statZF)
       
    end
    fprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',...
        '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n'])
end