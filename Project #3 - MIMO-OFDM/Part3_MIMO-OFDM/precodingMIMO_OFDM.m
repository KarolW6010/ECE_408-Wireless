function ber = precodingMIMO_OFDM(ants,M,R,SNRdB,data,H,DataBits)
%{
Sends the data through a flat fading channel and uses precoding and
receiver shaping. Returns the BER.

Inputs:
    ants        : (1, ) Number of Tx and Rx antennas (ants by ants system)
    M           : (1, ) M-ary QAM (2,4,16,64)
    R           : (1, ) ECC code rate (1/2, 2/3, 3/4)
    SNRdB       : (1, ) SNR in dB
    data        : (N,ants) Data to be sent over the channel
    H           : {N}(ants,ants) Cell array containing MIMO channels
    DataBits    : (log2(M)*R*N*48,ants) Binary data (original)

Outputs:
    ber         : (1, ) BER for the data
%}

N = size(data,1);               %Number of sent symbols (N/80 OFDM symbols)

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

%First perform the SVD of each channel
U = cell(1,N);
S = cell(1,N);
V = cell(1,N);
for ii = 1:N
    [u,s,v] = svd(H{ii});
    U{ii} = u;
    S{ii} = diag(s);
    V{ii} = v;
end

%Precode the data and send through the channel
preData = nan(ants,N);
appChan = nan(ants,N);
for ii = 1:N
    preData(:,ii) = V{ii}*data(ii,:).';
    appChan(:,ii) = H{ii}*preData(:,ii);
end

%Add AWGN
noise = nan(ants,N);
for ii = 1:ants
    sigPower = mean(abs(appChan(ii,:)).^2);
    noiseVar = sigPower*10^(-SNRdB/10);
    noise(ii,:) = sqrt(noiseVar/2)*(randn(1,N) + 1j*randn(1,N));
end
rx = appChan + noise;       %Received signal

%Perform receiver shaping
recModData = nan(size(data));
for ii = 1:N
    %Note the diag term here is scaling things so that the symbols can be 
    %properly decoded but the noise is being appropriately scaled to 
    %reflect the channel gain. Basically the signal isn't getting
    %stronger but the noise is getting weaker.
    recModData(ii,:) = (diag(1./S{ii})*U{ii}'*rx(:,ii)).';
end

%Decode
decBits = nan(size(DataBits));
for ii = 1:ants
    decBits(:,ii) = OFDMdemod(recModData(:,ii),R,M);
end
    
ber = mean(decBits  ~= DataBits,'all');











