function [DataBits,sym] = OFDMsymGen(N,R,M)
%{
Generates an N OFDM symbols according to the 802.11a standard with error
correction code rate R and M-ary QAM.

Inputs:
    N       : (1, ) Number of OFDM symbols to generate
    R       : (1, ) ECC code rate (1/2, 2/3, 3/4)
    M       : (1, ) M-ary QAM (2,4,16,64)

Outputs:
    binData : (log2(M)*N*48*R, ) Binary Data that was sent
    sym     : (N*80, ) OFDM symbol
%}

if(sum(M == [2,4,16,64]) ~= 1)
    error('Choose valid modulation 2, 4, 16, 64 - QAM')
end

nCodedBits = log2(M)*N*48;
nDataBits = nCodedBits*R;

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

%Generate and encode data 
DataBits = randi([0,1],nDataBits,1);
CodedBits = convEncoder(DataBits);

%Interleave
nBits = log2(M)*48;
s = log2(M)/2;
init = (0:(nBits-1)).';
int1 = (nBits/16)*mod(init,16) + floor(init/16);
int2 = s*floor(int1/s) + mod((int1+nBits-floor(16*int1/nBits)),s) + 1;
Inter = repmat(int2,N,1) + repelem((0:(N-1)).'*nBits,nBits,1);

InterBits = CodedBits(Inter); 

%Seperate bits into groups of log2(M) and modulate
SymBits = reshape(InterBits,log2(M),[]);
Syms = qammod(SymBits,M,'InputType','bit','UnitAveragePower',true);
Syms = reshape(Syms,48,[]);

%Now do the IFFT
preIFFT = nan(64,N);
preIFFT( 1   ,:) = 0;                  %Zero 
preIFFT( 2: 7,:) = Syms( 1: 6,:);      %Data
preIFFT( 8   ,:) = 1;                  %Pilot
preIFFT( 9:21,:) = Syms( 7:19,:);      %Data
preIFFT(22   ,:) = 1;                  %Pilot
preIFFT(23:27,:) = Syms(20:24,:);      %Data
preIFFT(28:38,:) = 0;                  %Zero
preIFFT(39:42,:) = Syms(25:28,:);      %Data
preIFFT(43   ,:) = 1;                  %Pilot
preIFFT(44:56,:) = Syms(29:41,:);      %Data
preIFFT(57   ,:) = 1;                  %Pilot
preIFFT(58:64,:) = Syms(42:48,:);      %Data

postIFFT = ifft(preIFFT,[],1);

%Add the cyclic prefix
OFDMsym = [postIFFT((end-15):end,:);postIFFT];
sym = OFDMsym(:);


















