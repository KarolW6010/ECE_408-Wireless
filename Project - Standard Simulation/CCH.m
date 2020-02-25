%Karol Wadolowski, Project: GSM PHY Layer Simulation
%This code generates the control channel data for one multiframe. No noise
%is added for this demonstration. This is just to show that the encoding,
%interleaving, formatting, etc. where implemented correctly.
clear; clc; close all;

%Create the 1/2 rate convolutional encoder and decoder
trellis = poly2trellis([5],[23,33]);
convEncoder = comm.ConvolutionalEncoder(trellis);

vitDecoder = comm.ViterbiDecoder(trellis,'InputFormat','Hard');
vitDecoder.TracebackDepth = 39;
vitDecoder.TerminationMethod = 'Terminated';

%Differential Decoder
diffDecoder = comm.DifferentialEncoder();

%Create modulator and demodulator
sps = 4;        %Samples per symbol
pLen = 4;       %GMSK filter length
gmskMod = comm.GMSKModulator('BitInput',true,'PulseLength',pLen,...
    'SamplesPerSymbol',sps,'InitialPhaseOffset',pi/4);
gmskDem = comm.GMSKDemodulator('BitOutput',true,'PulseLength',pLen,...
    'SamplesPerSymbol',sps,'InitialPhaseOffset',pi/4);

%The control control multiframe burst order
const = 'FSCCCCCCCC';
bursts = ['FSBBBBCCCC',const,const,const,const];

%Generate some data
bsic = randi([0,2^6-1],1,1);    %Base transciever station identity code
fprintf('BSIC is %d\n\n',bsic)
BCbursts = randi([0,1],10,184); %Generate data for BC bursts


%Generate a random frame number (at the beginning of a multiframe)
%These are the initial values
FN = floor(randi([0,2715647],1,1)/2048)*2048;    
T1 = floor(FN/1326);            %Which superframe
T2 = mod(FN,26);                %Which multiframe
T3 = floor((mod(FN,51)-1)/10);  %Where in the multiframe
if(T3 < 0)
    T3 = 0;
end
fprintf('The starting frame number is %d\n',FN)
fprintf('\tCorresponding T1 %d\n',T1)
fprintf('\tCorresponding T2 %d\n',T2)
fprintf('\tCorresponding T3 %d\n\n',T3)

%SCH training sequence
trainSeq = [1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,...
            0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1,...
            0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1,...
            1, 0, 1, 1]';

len = length(bursts);

%Generator polynomials for SCH, BCCH, and CCCH bursts
divSCH = zeros(1,11);   divSCH([1,3,5,6,7,9,11]) = 1;
divBC  = zeros(1,41);   divBC([1,4,18,25,27,41]) = 1;

%BC interleaving
k = 0:455;
B = mod(k,4) + 1;
ind = 2*mod(49*k,57) + floor(mod(k,8)/4) + 1;

bitMessage = nan(156,len);

%Convert the letters and data to an actual bit stream
countBC = 0;
ii = 1;
while(ii <= len)
    convEncoder.reset();        %Encoder reset
    
    switch bursts(ii)
        case{'B','C'}
            %BCCH or CCCH burst (treated the same but info is technically
            %supposed to be different.
            countBC = countBC + 1;
            temp = [zeros(1,40),BCbursts(countBC,:)];
            
            %Systematically encode the data
            [~,r] = gfdeconv(temp,divBC,2);
            parity = [1-r,ones(1,40-length(r))];
            temp1 = temp;
            temp1(1:40) = parity;
            
            %Do a check to see if meets the criteria
            [~,rr] = gfdeconv(temp1,divBC,2);
            if(~isequal(rr,ones(1,40)))
                error('Parity incorrect!')
            end
            
            %Encoded data
            encBits = [temp1,zeros(1,4)];
            
            %Conv Encoded data
            convEncBits = convEncoder(encBits');
            
            %Interleaving
            inter = nan(114,4);
            for jj = 1:456
                inter(ind(jj),B(jj)) = convEncBits(jj);
            end
            
            bitMessage(1:3,ii:(ii+3)) = zeros(3,4);             %Tail bits
            bitMessage(4:60,ii:(ii+3)) = inter(1:57,:);         %Encoded data
            bitMessage(61:88,ii:(ii+3)) = zeros(28,4);          %Stealing flags + training sequence
            bitMessage(89:145,ii:(ii+3)) = inter(58:end,:);     %Encoded data
            bitMessage(146:148,ii:(ii+3)) = zeros(3,4);         %Tail bits
            bitMessage(149:156,ii:(ii+3)) = ones(8,4);          %Guard bits
            %Note that the training sequence is not set because I couldn't
            %find which one to use. There where 8 possible sequences but it
            %was specified which one to use. 
            
            ii = ii + 4;
            
        case{'F'}
            %FCCH burst
            bitMessage(1:148,ii) = zeros(148,1);
            bitMessage(149:156,ii) = ones(8,1);
            
            ii = ii + 1;
            
        case{'S'}
            %SCH burst
            curFN = FN + ii -1;         %Current frame number
            curT1 = floor(FN/1326); 
            curT2 = mod(curFN,26);
            curT3 = floor((mod(curFN,51)-1)/10);
            if(curT3 < 0)
                curT3 = 0;
            end
            
            fprintf('Sending frame number is %d\n',curFN)
            fprintf('\tCorresponding T1 %d\n',curT1)
            fprintf('\tCorresponding T2 %d\n',curT2)
            fprintf('\tCorresponding T3 %d\n\n',curT3)
            
            %Systematically encode the data
            temp = [zeros(1,10),...
                    de2bi(bsic,  6,'left-msb'),...
                    de2bi(curT1,11,'left-msb'),...
                    de2bi(curT2, 5,'left-msb'),...
                    de2bi(curT3, 3,'left-msb')];
            
                
            [~,r] = gfdeconv(temp,divSCH,2);
            parity = [1-r,ones(1,10-length(r))];
            temp1 = temp;
            temp1(1:10) = parity;
            
            %Do a check to see if meets the criteria
            [~,rr] = gfdeconv(temp1,divSCH,2);
            if(~isequal(rr,ones(1,10)))
                error('Parity incorrect!')
            end
                    
            %Encoded data
            encBits = [temp1,zeros(1,4)];
            
            %Conv Encoded data
            convEncBits = convEncoder(encBits');
            
            %Now place in the frame
            bitMessage(1:3,ii) = zeros(3,1);                %Tail bits
            bitMessage(4:42,ii) = convEncBits(1:39);        %Encoded data
            bitMessage(43:106,ii) = trainSeq;               %Training sequence
            bitMessage(107:145,ii) = convEncBits(40:end);   %Encoded data
            bitMessage(146:148,ii) = zeros(3,1);            %Tail bits
            bitMessage(149:156,ii) = ones(8,1);             %Guard bits
            
            ii = ii + 1;
            
        otherwise
            error('Something went wrong!')
    end
end

%Modulate
formatted = [bitMessage(:);zeros(gmskDem.TracebackDepth,1)];
diffEnc = xor(formatted,[0;formatted(1:(end-1))]);
sent = gmskMod(diffEnc);

%Demodulate
rec = gmskDem(sent);
diffDec = diffDecoder(rec);
recBits = reshape(diffDec((gmskDem.TracebackDepth+1):end),156,[]);

fprintf('*******************Receiving*******************\n\n')
%Can't print the received frame number. This is because MATLAB doesn't have
%a function that performs the Chinese Remainder theorem and I am unfamiliar
%with how to implement such a function. Should have taken discrete =(

%Recover the original data
retBCbursts = nan(10,184);
countBC = 0;
ii = 1;
while(ii <= len)
    vitDecoder.reset();         %Decoder reset
    
    switch bursts(ii)
        case{'B','C'}
            %BCCH or CCCH burst
            countBC = countBC + 1;
            
            rel = recBits(:,ii:(ii+3));     %Relevant bits
            
            encData = [rel(4:60,:);rel(89:145,:)];
            
            %Inverse interleaving
            unInter = nan(456,1);
            for jj = 1:456
                unInter(jj) = encData(ind(jj),B(jj));
            end
            
            %Undo the convolutional code
            decData = vitDecoder(unInter(:))';
            
            %Take the systematic bits
            retBCbursts(countBC,:) = decData(41:(end-4));
              
            ii = ii + 4;
            
        case{'F'}
            %FCCH burst
            %No data here so do nothing
            ii = ii + 1;

        case{'S'}
            %SCH burst
            rel = recBits(:,ii);            %Relevant bits
            
            encData = [rel(4:42);rel(107:145)];
            decData = vitDecoder(encData)';
            
            %Retrieve the BSCI, T1, T2, and T3 values. Remember data is
            %encoded systematically and the first 10 bits are the parity.
            retBSCI = bi2de(decData((1:6)+10),'left-msb');
            retT1 = bi2de(decData((7:17)+10),'left-msb');
            retT2 = bi2de(decData((18:22)+10),'left-msb');
            retT3 = bi2de(decData((23:25)+10),'left-msb');
            
            fprintf('Receiving frame\n')
            fprintf('\tCorresponding T1 %d\n',retT1)
            fprintf('\tCorresponding T2 %d\n',retT2)
            fprintf('\tCorresponding T3 %d\n\n',retT3)
            
            ii = ii + 1;
            
        otherwise
            error('Should literally be impossible to get here!')
    end
   
end

fprintf('Received BSCI is %d\n\n', retBSCI)

if(isequal(retBCbursts,BCbursts))
    fprintf('All BCCH and CCCH data has been recovered successfully!\n')
end
    






















