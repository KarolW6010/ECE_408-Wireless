function [DataBits,data,H,DR] = genMIMO_OFDMdata(ants,M,R,N,fD)
%{
Generates MIMO-OFDM data for the system given the parameters. Use a
sampling rate of 20 MHz (in line with 802.11a).

Inputs:
    ants        : (1, ) Number of Tx and Rx antennas (ants by ants system)
    M           : (1, ) M-ary QAM (2,4,16,64)
    R           : (1, ) ECC code rate (1/2, 2/3, 3/4)
    N           : (1, ) Number of OFDM symbols to generate 
    fD          : (1, ) Max doppler frequency for flat fading channel

Outputs:
    DataBits    : (log2(M)*N*48*R,ants) Binary data for each of the antennas
    data        : (N*80,ants) The symbols that are being sent
    H           : {N*80}(ants,ants) Cell array containing the MIMO channels
    DR          : (1, ) The data rate of the system (bits per second)
%}

fs = 2e7;           %Sampling rate of 20 MHz
T = (N*80)/fs;      %Duration of the samples

%Generate Data
DataBits = nan(log2(M)*N*48*R,ants);
data = nan(N*80,ants);
for ii = 1:ants
    [DataBits1,data1] = OFDMsymGen(N,R,M);
    DataBits(:,ii) = DataBits1;
    data(:,ii) = data1;
end

%Construct the MIMO flat fading channels
channels = cell(ants,ants);
N1 = N*80;
H = cell(1,N1);
for ii = 1:ants
    for jj = 1:ants
        if(2*T*fD < 2)
            k = fs/(N1*fD);
            r = getRayleighFading(k*T,fD,k*N1);
            r = r(1:N1);
        else
            if(mod(T*fD,1) == 0)
                r = getRayleighFading(T,fD,N1);
            else
                k = mod(T*fD,1);
                r = getRayleighFading(round(T/k,5),fD,round(N1/k,5));
                r = r(1:N1);
            end
        end
        channels{ii,jj} = r;
    end
end
for ii = 1:N1
    for jj = 1:ants
        for kk = 1:ants
            H{ii}(jj,kk) = channels{jj,kk}(ii);
        end
    end
end

%Calculate Data Rate
DR = 48*R*log2(M)*fs*ants/80;

    