function [f, I] = HizhAll(w, D, G, T, P)
w = cast(w, 'double'); % QLM mode frequency, cm-1
D = cast(D, 'double'); % QLM frequency change, cm-1
G = 0.5*cast(G, 'double'); % QLM mode decay width, cm-1, !! IT IS HALF OF THE LOW T-Model Gamma !!
T = cast(T, 'double'); % Model temperature, K
P = cast(P, 'double'); % Emission line position, cm-1
%----------------------------------------------------------------

nbar = 1/(exp(1.44*w/T)-1); % QLM population; 1.44 is for cm-1 -> K conversion
b = sqrt( G - D^2/4 - 1i*G*D*(2*nbar+1) );
if ( real(b) > 0 )
    beta = b;
else
    beta = -b;
end

lambda = G - 1i*D/2 + beta;
alpha = -1i*nbar*D/lambda;

%==========Natural linewidth: gamma_0 FWHM ===============
g0 = 0.003; % FWHM, cm-1 or 0.90 GHz;

% generates Spectrum according to the All-Temperatures model
Fs = 6000; %GHz sampling frequency
tmax = 200; %ns integration time
numPoints = 1 + tmax*Fs; % number of data points in time domain
t = linspace(0, tmax, numPoints);
k = 30; % cm-1 to GHz=1/ns conversion
%--Spectrum using new, all-T model:
y = exp( -pi*k*t.*g0 + 2*pi*1i*t.*P*k+2*pi*t.*1i*k*D*(nbar+1.0/2)-2*pi*t.*1i*k*D*(nbar+1)*alpha/(1+alpha)+(1i*D*(nbar+1)/(lambda*(1+alpha))*log(1+alpha-alpha*exp(-2*pi*t.*k*lambda))));

NFFT = 2^nextpow2(numPoints);
Y = fft(y,NFFT);
Y = Y/max(Y);
%f = 1/k*Fs/2*linspace(0,1,NFFT/2 + 1); % frequency domain spectrum
f = 1/k*Fs/2*linspace(-1,1,NFFT); % frequency domain spectrum
%I = Y(1:NFFT/2 + 1).*conj(Y(1:NFFT/2 + 1));
I = Y(1:NFFT).*conj(Y(1:NFFT));
