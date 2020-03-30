%Karol Wadolowski, Project #2: CDMA
clear; clc; close all;

%%
load('Rcvd_Wadolowski.mat')

%Get the Walsh orthogonal codes
walsh = hadamard(8);
w0 = walsh(1,:);        %Walsh orthogonal code 0
w5 = walsh(5,:);        %Walsh orthogonal code 5

lenR = length(Rcvd);                    %Number of samples
frames = lenR/(4*255);                  %Get the number of frames sent
rrcFilt = rcosdesign(.75,6,4,'sqrt');   %Get the Root Raised Cosine filter
rFiltered = conv(Rcvd,rrcFilt,'same');  %Filter the received signal
rcvd4 = rFiltered(1:4:lenR);            %Take every fourth sample

%Find and generate the correct PN sequence. To do this the initial state of
%the LFSR has to be found. Applying the PN of the correct initial state
%will maximize mean distance from the origin of the initial frame
maxVal = -inf;                          %Current maxVal
maxState = zeros(1,8);                  %State corresponding to this value
pilotDist = nan(1,255);                 %Mean distance from origin for each
                                        %init state
for ii = 1:255
    %Loop through all possible initial states of LFSR
    state = de2bi(ii,8,'left-msb');
    seq = lfsrSeq([8,7,6,1],state,255);
    seq = 2*seq - 1;
    
    %Apply the PN sequence to the pilot samples
    temp = rcvd4(1:255).*seq;
    
    %Check if this initial state is the best one so far
    pilotDist(ii) = abs(mean(temp));
    
    if(pilotDist(ii) > maxVal)
        %If best so far save value and state
        maxVal = pilotDist(ii);
        maxState = state;
    end
end
       
init = maxState;                        %Initial state of LFSR
seq = lfsrSeq([8,7,6,1],init,lenR/4);   %PN sequence
seq = 2*seq - 1;                        %Correspond bits to BPSK

rPN = rcvd4.*seq;                       %Undo the PN sequence

%Get rid of the zero fill samples
temp = reshape(rPN,255,[]);             %Reshape for convenience 
infoSamps = temp(1:192,:).';            %Get rid of zero fill
%InfoVals is all the samples containing the pilot and data information.
%Each row contains the 192 information samples from a frame.

%Obtain the spread pilot and data components from the information samples
pilot = repmat(w0,frames,24).*infoSamps;
data  = repmat(w5,frames,24).*infoSamps;

%Need to aquire the phase of the pilot bits for each frame in order to
%decode the data stream correctly. This will also handle the frequency
%offset.
phases = angle(mean(pilot,2));      

%Adjust the phase of the information values 
corInfoSamps = infoSamps.*repmat(exp(-1j*phases),1,192);
corInfoSamps = reshape(corInfoSamps.',8,[]);

%Perform inner products of groups of 8 bits with the Walsh codes to get the
%pilot and data streams. Divide by 8 to normalize.
iwPilot = w0*corInfoSamps/8;
iwData  = w5*corInfoSamps/8;

%Decode the data stream and convert to ascii
decData = reshape(real(iwData),8,[])' < 0;
decData = bi2de(decData,'right-msb')';

%Don't include first and last frame as those have no actaual data
message = char(decData(4:(end-3)))

%Now to correspond the phase changes to a frequency offset
chipRate = 1e6;
frameRate = chipRate/255;
unPhases = unwrap(phases.');    %Unwrapped phase
diffPh = diff(unPhases);        %Phase difference among adjacent frames

%Approximate frequency using df/dt = 1/(2pi) dtheta/dt
freqOffset = mean(diffPh)*frameRate/(2*pi);
fprintf(['\nFor a chip rate of 1e6 MHz, the frequency offset is ',...
    num2str(freqOffset),' Hz\n']);

%% Now time for some plots
figure('units','normalized','outerposition',[0,0,1,1])
subplot(1,2,1)
scatter(real(rcvd4),imag(rcvd4))
xlim([-2,2]);   ylim([-2,2]);   grid on;    daspect([1 1 1]);
title('Received Samples Before Undoing PN Sequence')

subplot(1,2,2)
scatter(real(rPN),imag(rPN))
xlim([-2,2]);   ylim([-2,2]);   grid on;    daspect([1 1 1]);
title('Received Samples After Undoing PN Sequence')

figure('units','normalized','outerposition',[0,0,1,1])
stem(1:255,pilotDist)
xlim([1,255]);  ylim([0,ceil(max(pilotDist))]); grid on;
xlabel('LFSR Initial State (Left-MSB)')
ylabel('Mean Distance from Origin')
title('Mean Distance of Initial Frame after Applying PN Sequence')

p1 = mean(reshape(pilot.',8,[]));
d1 = mean(reshape(data.',8,[]));
figure('units','normalized','outerposition',[0,0,1,1])
subplot(1,2,1)
scatter(real(p1),imag(p1))
xlim([-2,2]);   ylim([-2,2]);   grid on;    daspect([1 1 1]);
title('Pilot')

subplot(1,2,2)
scatter(real(d1),imag(d1))
xlim([-2,2]);   ylim([-2,2]);   grid on;    daspect([1 1 1]);
title('Data')
sgtitle('Before Phase Adjustment')

figure('units','normalized','outerposition',[0,0,1,1])
subplot(1,2,1)
scatter(real(iwPilot),imag(iwPilot))
xlim([-2,2]);   ylim([-2,2]);   grid on;    daspect([1 1 1]);
title('Pilot')

subplot(1,2,2)
scatter(real(iwData),imag(iwData))
xlim([-2,2]);   ylim([-2,2]);   grid on;    daspect([1 1 1]);
title('Data')
sgtitle('After Phase Adjustment')

figure('units','normalized','outerposition',[0,0,1,1])
plot(1:frames,phases,'-x','Linewidth',2)
hold on;
plot(1:frames,unPhases,'-x','Linewidth',2)
legend('Wrapped','Unwrapped','Location','NorthWest')
xlabel('Frame')
ylabel('Phase Offset (rad)')
title('Phase Offset Present In Each Frame')
grid on


