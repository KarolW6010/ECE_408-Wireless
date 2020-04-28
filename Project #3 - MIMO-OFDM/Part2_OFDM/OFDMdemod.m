function DecBits = OFDMdemod(sym,R,M)
%{
Takes a stream of OFDM symbols and decodes it.

Inputs:
    sym         : (80*N, ) N OFDM symbols of length 80
    R           : (1, ) ECC code rate (1/2, 2/3, 3/4)
    M           : (1, ) M-ary QAM (2,4,16,64)

Outputs:
    DecBits     : (log2(M)*N*48*R, ) Decoded bits
%}

if(sum(M == [2,4,16,64]) ~= 1)
    error('Choose valid modulation 2, 4, 16, 64 - QAM')
end

N = length(sym)/80;

%Remove cyclic prefix and Perform FFT
temp = reshape(sym,80,[]);
preFFT = temp(17:end,:);
postFFT = fft(preFFT,[],1);

%Remove the pilots and zeros to get the symbols
Syms = nan(48,N);
Syms( 1: 6,:) = postFFT( 2: 7,:);
Syms( 7:19,:) = postFFT( 9:21,:);
Syms(20:24,:) = postFFT(23:27,:);
Syms(25:28,:) = postFFT(39:42,:);
Syms(29:41,:) = postFFT(44:56,:);
Syms(42:48,:) = postFFT(58:64,:);


%Demodulate
SymBits = qamdemod(Syms,M,'OutputType','bit','UnitAveragePower',true);
SymBits = SymBits(:);

%Deinterleave to get the encoded bits
nBits = log2(M)*48;
s = log2(M)/2;
init = (0:(nBits-1)).';
int1 = (nBits/16)*mod(init,16) + floor(init/16);
int2 = s*floor(int1/s) + mod((int1+nBits-floor(16*int1/nBits)),s) + 1;
Inter = repmat(int2,N,1) + repelem((0:(N-1)).'*nBits,nBits,1);
[~,deInter] = sort(Inter);

CodedBits = SymBits(deInter);

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

DecBits = vitDecoder(CodedBits);





