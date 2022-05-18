%
%		example.m
%		Douglas L. Jones
%		University of Illinois
%		September 7, 2006
%
%	example.m: Creates a stochastic signal and analyzes spectrum
%		for spectrum analysis example.
%
%

%   SETUP

datalen = 10000;
postdatalen = 8192;
fftlen = 1024;

freqsamps = 2*pi*[0:fftlen/2-1]/fftlen;
dBlimit = -80;


s1freq = pi*0.2;
s1amp = 0.3;

noiseamp = 1.0;

signal = s1amp*sin(s1freq*[0:datalen-1]) + noiseamp*randn(size([0:datalen-1]));
[B,A] = butter(10,0.5);
signal = filtfilt(B,A,signal);
signal = signal((datalen-postdatalen)/2:(datalen-postdatalen)/2+postdatalen-1);
datalen = postdatalen;

figure(1)
plot([0:length(signal)-1],signal)
xlabel('Sample number')
ylabel('Amplitude')
print -dpng examplesig

winlen = 64;
boxcar = ones(1,winlen);
figure(2)
freqmag = abs(fft(boxcar.*signal(1:winlen),fftlen));
%freqmag = freqmag/max(freqmag);
freqmag = 20*log10(freqmag(1:fftlen/2));
freqmag = max(freqmag,dBlimit*ones(size(freqmag)));
plot(freqsamps,freqmag)
xlabel('DTFT frequency')
ylabel('magnitude (dB)')
print -dpng stoch64

winlen = 1024;
boxcar = ones(1,winlen);
figure(3)
freqmag = abs(fft(boxcar.*signal(1:winlen),fftlen));
%freqmag = freqmag/max(freqmag);
freqmag = 20*log10(freqmag(1:fftlen/2));
freqmag = max(freqmag,dBlimit*ones(size(freqmag)));
plot(freqsamps,freqmag)
xlabel('DTFT frequency')
ylabel('magnitude (dB)')
print -dpng stoch1024


winlen = 64;
boxcar = ones(1,winlen);
figure(4)

spec = zeros(1,fftlen);
avcnt = 0;
for ii=1:winlen/4:datalen-winlen,  % average spectra of overlapped blocks
  spec = spec + abs(fft(boxcar.*signal(ii:ii+winlen-1),fftlen)).^2;
  avcnt = avcnt + 1;
end

spec = spec/avcnt;
freqmag = 10*log10(spec(1:fftlen/2));
freqmag = max(freqmag,dBlimit*ones(size(freqmag)));
plot(freqsamps,freqmag)
xlabel('DTFT frequency')
ylabel('magnitude (dB)')
print -dpng stochPSD


%  DONE

