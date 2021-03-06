function [r,t] = rayleighFading2(T,fD,Npt)
%{
Generateas a random rayleigh fading envelope of duration T generated using 
Fig. 5.24 from Rapaport and the 7 associated steps.

Inputs:
    T   : (1, ) Duration of Rayleigh fading envelope
    fD  : (1, ) Maximum doppler frequency
    Npt : (1, ) Multiply IFFT size by this amount

Outputs:
    r   : (N*Npt, ) Rayleigh fading channel
%}

%First we need to get the appropiate amount of points
delF = 1/T;             %Frequency spacing
N  = 2*fD/delF + 1;     %N point rayleigh fading envelope
N = round(N,10);

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
r = sqrt(sum(postifft,2));                      %Add the two branches 
r = r/mean(r);                                  %Divide by mean so average 0 dB power

t = linspace(0,T,N*Npt);                        %Time axis
