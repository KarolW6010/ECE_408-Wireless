function [berMMSE,berZF] = mmse_zfMIMO(ants,M,R,SNRdB,data,H,DataBits)
%{
Sends the data through a flat fading channel and uses an MMSE equalizer and
a Zero Forcing equalizer separately.  Returns the BER for both cases

Inputs:
    ants        : (1, ) Number of Tx and Rx antennas (ants by ants system)
    M           : (1, ) M-ary QAM (2,4,16,64)
    R           : (1, ) ECC code rate (1/2, 2/3, 3/4)
    SNRdB       : (1, ) SNR in dB
    data        : (N,ants) Data to be sent over the channel
    H           : {m}(ants,ants) Cell array containing N MIMO channels
    DataBits    : (log2(M)*R*N,ants) Binary data (original)

Outputs:
    berMMSE     : (1, ) BER for the data using MMSE equalizer
    berZF       : (1, ) BER for the data using ZF equalizer
%}

N = size(data,1);               %Number of sent symbols

%Select Proper Convolutional Decoder
trellis = poly2trellis([7],[133,171]);
vitDecoder = comm.ViterbiDecoder(trellis,'InputFormat','Hard');
vitDecoder.TerminationMethod = 'Truncated';
switch R
    case(1/2)
        
    case(2/3)
        vitDecoder.PuncturePatternSource = 'Property';
        vitDecoder.PuncturePattern = [1;1;1;0];
        
    case(3/4)
        vitDecoder.PuncturePatternSource = 'Property';
        vitDecoder.PuncturePattern = [1;1;1;0;0;1];

    otherwise
        error('Choose one of available encoding rates!')
end

%Apply the channel
appChan = nan(ants,N);
for ii = 1:N 
    appChan(:,ii) = H{ii}*data(ii,:).';
end

%Add AWGN
noise = nan(ants,N);
for ii = 1:ants
    sigPower = mean(abs(appChan(ii,:)).^2);
    noiseVar = sigPower*10^(-SNRdB/10);
    noise(ii,:) = sqrt(noiseVar/2)*(randn(1,N) + 1j*randn(1,N));
end
rx = appChan + noise;       %Received signal

%Apply the equalizers
rxMMSE = nan(N,ants);
rxZF = nan(N,ants);
for ii = 1:N
    temp = H{ii}'*H{ii};
    temp2 = H{ii}'*rx(:,ii);
    rxMMSE(ii,:) = (((temp + 10^(-SNRdB/10)*eye(ants))^-1)*temp2).';
    rxZF(ii,:) = ((temp^-1)*temp2).';
end

%Demodulate the symbols
deModMMSE = nan(log2(M)*N,ants);
deModZF = nan(log2(M)*N,ants);
for ii = 1:ants
    temp = rxMMSE(:,ii).';
    temp1 = qamdemod(temp,M,'OutputType','bit','UnitAveragePower',true).';
    deModMMSE(:,ii) = temp1(:);
    
    temp = rxZF(:,ii).';
    temp1 = qamdemod(temp,M,'OutputType','bit','UnitAveragePower',true).';
    deModZF(:,ii) = temp1(:);
end

%Deinterleave to get the encoded bits
nCodedBits = log2(M)*N;
int = reshape(1:100,10,10).';
int = int(:);
Inter = repmat(int,nCodedBits/1e2,1) + repelem((0:(nCodedBits/1e2-1)).'*100,100,1);
[~,deInter] = sort(Inter);
CodedBitsMMSE = deModMMSE(deInter,:);
CodedBitsZF = deModZF(deInter,:);

%Deocde
decBitsMMSE = nan(size(DataBits));
decBitsZF = nan(size(DataBits));
for ii = 1:ants
    decBitsMMSE(:,ii) = vitDecoder(CodedBitsMMSE(:,ii));
    vitDecoder.reset();
    
    decBitsZF(:,ii) = vitDecoder(CodedBitsZF(:,ii));
    vitDecoder.reset();
end

berMMSE = mean(decBitsMMSE ~= DataBits,'all');
berZF = mean(decBitsZF ~= DataBits,'all');




















