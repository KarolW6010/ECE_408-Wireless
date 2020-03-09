function pltRayleighFading(T,fD,Npt)
%{
Plots a random rayleigh fading envelope of duration T generated using Fig.
5.24 from Rapaport and the 7 associated steps.

Inputs:
    T   : (1, ) Duration of Rayleigh fading envelope
    fD  : (1, ) Maximum doppler frequency
    Npt : (1, ) Multiply IFFT size by this amount

Outputs:
    A graph of the rayleigh fading envelope.
%}

%First we need to get the appropiate amount of points
delF = 1/T;             %Frequency spacing
N  = 2*fD/delF + 1;     %N point rayleigh fading envelope

if(mod(N,1) ~= 0)
    error('Choose fD and T such that 2*fD*T is an integer')
end

%Clarke spectrum of length N
[~,S] = clarkeSpec(N,fD);

%Generate 2 sets of random gaussian samples
gSamps = sqrt(1/2)*(randn(ceil(N/2),2)+1j*randn(ceil(N/2),2));

%Apply the conjugate symmetry
if(mod(N,2) == 1)
    %N odd
    g = [conj(flipud(gSamps(2:end,:))); gSamps];
else
    %N even
    g = [conj(flipud(gSamps)); gSamps];
end

b4ifft = g.*repmat(sqrt(S),1,2);                %Multiply by square root of Clarke Spectrum
postifft = ifft(ifftshift(b4ifft),N*Npt,1).^2;  %Take IFFT and square it
r = abs(sqrt(sum(postifft,2)));                 %Add the two branches and take absolute value to get envelope
r = r/mean(r);                                  %Divide by mean so average 0 dB power

t = linspace(0,T,N*Npt);                        %Time axis

%Plotting
figure
plot(t,10*log10(r),'Linewidth',2)
title(['Rayleigh Fading Envelope (f_D = ', num2str(fD),' Hz)'])
xlabel('Elapsed Time (s)')
ylabel('Signal Level (dB)')
xlim([0,T])
grid on


