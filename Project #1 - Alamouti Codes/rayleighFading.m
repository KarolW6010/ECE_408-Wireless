function r = rayleighFading(N)
%{
Returns a time signal of length N of a rayleigh fading channel generated 
using Fig. 5.24 from Rapaport and the 7 associated steps. Note that the 
maximum doppler frequency does't matter since the IFFT is inherently 
unitless. If a doppler frequency was specified then we would have time 
info. But that is not needed for the task at hand.

Inputs:
    N   : (1, ) Desired length of Rayleigh Fading envelope

Outputs:
    r   : (N,1) Rayleigh fading signal

Note: Taking the absolute value of r will give the rayleigh fading
envelope.
%}

%Clarke spectrum of length N.
[~,S] = clarkeSpec(N,1);

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

b4ifft = g.*repmat(sqrt(S),1,2);            %Multiply by square root of Clarke Spectrum
postifft = ifft(ifftshift(b4ifft),N,1).^2;  %Take IFFT and square it
r = sqrt(sum(postifft,2));                  %Add the two branches
r = r/mean(r);                              %Divide by mean so average 0 dB power


