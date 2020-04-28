function [DataBits,data,H,DR] = genMIMOdata(ants,M,R,N,fD)
%{
Generates MIMO data for the system with the given parameters. Assume the
sampling rate is 20 MHz. Meaning 2e7 symbols per second

Inputs:
    ants        : (1, ) Number of Tx and Rx antennas (ants by ants system)
    M           : (1, ) M-ary QAM (2,4,16,64)
    R           : (1, ) ECC code rate (1/2, 2/3, 3/4)
    N           : (1, ) Number of symbols to generate (multiple of 1e2)
    fD          : (1, ) Max doppler frequency for flat fading channel

Outputs:
    DataBits    : (log2(M)*R*N,ants) Binary data for each of the antennas
    data        : (N,ants) The symbols that are being sent
    H           : {N}(ants,ants) Cell array containing the N MIMO channels
    DR          : (1, ) The data rate of the system (bits per second)
%}

fs = 2e7;           %Sampling rate of 20 MHz
T = N/fs;           %Duration of the samples

if(mod(N,3e2)~=0)
    error('Choose N to be divisible by 300')
end

if(sum(M == [2,4,16,64]) ~= 1)
    error('Choose valid modulation 2, 4, 16, 64 - QAM')
end

%Select Proper Convolutional Encoder
trellis = poly2trellis([7],[133,171]);
convEncoder = comm.ConvolutionalEncoder(trellis);
switch R
    case(1/2)
        
    case(2/3)
        convEncoder.PuncturePatternSource = 'Property';
        convEncoder.PuncturePattern = [1;1;1;0];
        
    case(3/4)
        convEncoder.PuncturePatternSource = 'Property';
        convEncoder.PuncturePattern = [1;1;1;0;0;1];

    otherwise
        error('Choose one of available encoding rates!')
end

nCodedBits = log2(M)*N;         %Number of coded bits for each antenna
nDataBits = nCodedBits*R;       %Number of data bits for each antenna

%Generate and encode the data for each antenna
DataBits = randi([0,1],nDataBits,ants);
CodedBits = nan(nCodedBits,ants);
for ii = 1:ants
    CodedBits(:,ii) = convEncoder(DataBits(:,ii));
    convEncoder.reset();
end

%Interleave Bits (Block Interleaver)
int = reshape(1:100,10,10).';
int = int(:);
Inter = repmat(int,nCodedBits/1e2,1) + repelem((0:(nCodedBits/1e2-1)).'*100,100,1);

InterBits = CodedBits(Inter,:);

%Convert the bits into symbols
Syms = nan(N,ants);
for ii = 1:ants
    temp = InterBits(:,ii);
    temp = reshape(temp,[],log2(M)).';
    Syms(:,ii) = qammod(temp,M,'InputType','bit','UnitAveragePower',true).';
end

data = Syms;

%Construct the MIMO flat fading channels
channels = cell(ants,ants);
H = cell(1,N);
for ii = 1:ants
    for jj = 1:ants
        %
        if(2*T*fD < 2)
            k = fs/(N*fD);
            r = getRayleighFading(k*T,fD,k*N);
            r = r(1:N);
        else
            r = getRayleighFading(T,fD,N);
        end
        channels{ii,jj} = r;
    end
end
for ii = 1:N
    for jj = 1:ants
        for kk = 1:ants
            H{ii}(jj,kk) = channels{jj,kk}(ii);
        end
    end
end

%Calculate the data rate
DR = fs*log2(M)*R*ants;










