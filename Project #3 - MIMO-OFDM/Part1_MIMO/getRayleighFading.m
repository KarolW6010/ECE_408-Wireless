function r = getRayleighFading(T,fD,N)
%{
Gets N samples of a rayleigh fading waveform with max doppler frequency
fD and duration T.

Inputs:
    T       : (1, ) Duration
    fD      : (1, ) Max Doppler Frequency
    N       : (1, ) Number of samples

Outputs:
    r       : (N, ) Rayleigh fading samples
%}

%Samples generated from rayleighFading.m with multipier 1
defSamps = 2*fD*T+1;        

%Multiplier needed to get desired number of samples
mult = ceil(N/defSamps);

%Get rayleigh Fading Samples
[rTemp,~] = rayleighFading2(T,fD,mult);

%Only the amount necessary
r = rTemp(1:N);
