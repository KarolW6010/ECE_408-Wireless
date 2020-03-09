function [freq,S] = clarkeSpec(N,fD)
%{
Returns an N point Clarke spectrum with maximum doppler frequency fD.

Inputs:
    N       : (1, ) Number of points to generate from Clarke spectrum
    fD      : (1, ) Maximum doppler frequency

Outputs:
    freq    : (N,1) Frequencies over which the spectrum was evaluated over
    S       : (N,1) Clarke spectrum evaluated over freq
%}

freq = linspace(-fD,fD,N)';                 %Generate frequencies
S = 1.5./(pi*fD*sqrt(1-(freq/fD).^2));      %Obtain Clarke spectrum

%Now handle the edge cases (+- fD leads to infinity)
%Use a linearization to obtain spectrum estimate at +- fD
fLin = [freq(2);freq(end-1)];               %Linearization frequencies
sPrime = 1.5*fLin./((fD^3)*(1-(fLin/fD).^2).^1.5);  %Derivative of S

%Now use the linearization to get an approximation of the Clarke Spectrum
%at +- fD
S(1) = sPrime(1)*(-fD-fLin(1)) + S(2);      %-fD
S(end) = sPrime(2)*(fD-fLin(2)) + S(end-1); %+fD
