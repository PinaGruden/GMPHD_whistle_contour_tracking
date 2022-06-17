function [Zset] = preprocess_getZset(win_width_s,dt,fs,x,freqrange,peak_thr)
%function preprocess_getZset is a function that applies several
%pre-processing steps to an audio file in order to obtain the spectral peak
% peak measurements (Zset). Preprocessing steps include click removal, FFT,
% median filter, and exponential MA filter (see Gillespie et al., 2013, 
%Gruden & White, 2016 for details)

% Inputs:
% win_width = window length in seconds,
% dt = time increment between windows in seconds
% fs = sampling frequency,
% x = raw data,
% freqrange = a 1 x 2 vector of lower and higher frequency range limits in
% Hz (for example freqrange=[300,6000] will search all frequencies between
% 300 Hz and 6kHz. This range should be selected to encompass all
% frequencies across which the signals of interest can occur).
% peak_thr = a threshold in dB used for detecting spectral peaks from a
% smoothed spectrogram.

% Outputs:
% Zset = cell array containing spectral peak frequency measurements; each
% cell corresponds to a time step

%Pina Gruden, 2016, Institute of Sound and Vibration Research (ISVR),
%University of Southampton


win_width = round(win_width_s*fs); %window length in samples
w=hanning(win_width);
slide_incr= round(dt*fs); %time increment between windows in samples


tresh=5; %threshold for click removal
p=6; %power for click removal
lambda=0.02; %used in exponential moving average
%(higher lambda discounts older observations faster and vice versa)

numstps = ceil((length(x)-win_width)/slide_incr); %Number of overlapping windows

df=fs/(win_width); %frequency bin width
freq = 0:df:fs/2-df;
%for peak detection only process frequencies specified in freqrange
range=freq>freqrange(1) & freq<freqrange(2);
freqn=freq(range);

start=1;Ytnorm=0;
Zset=cell(1,numstps);Ynorm = zeros(size(freq,2),numstps); 

for i=1:numstps %do sliding window
    
    xi= x(start:start+win_width-1);
    xi=xi(:);
    w=w(:);
    
    %----------------Remove clicks from the signal---------------------
    sd=std(xi);
    mu=mean(xi);
    wi=1./(1+((xi-mu)/(tresh*sd)).^p);% Weighting function
    xn = xi.*wi;  % weighted signal (declicked)
    
    %-----------FFT, median filter and exponential MA filter-----------
    Y = fft(xn.*w,win_width); % FFT with hanning window
    L=length(Y);
    Y=Y(1:(L/2),:); %take fft up to Nyquist
    Ymag = 20*log10(abs(Y)); % Magnitude squared on a dB scale
    Ymedfilt = Ymag - medfilt1(Ymag,61); % Normalization across frequency
    Ytnorm = (1-lambda)*Ytnorm+lambda*Ymedfilt;%Normalize across time with exponential moving average
    Ynorm(:,i) = Ymedfilt - Ytnorm;
    start=start+slide_incr;
end

%--------- Gaussian smoothing (convolve spectrogram with kernel G) ------
% G = [1,2,1;2,4,2;1,2,1]/16;
% Ynorm_sm1=conv2(Ynorm,G,'same'); 


%----------------- Find spectral peaks -----------------------
for n=1:size(Ynorm,2)
    Ynormn=Ynorm(range,n);
    warning('off', 'signal:findpeaks:largeMinPeakHeight') % to turn off warning of no peaks greater than MinPeakHeight
    [~,locs]=findpeaks(Ynormn,'MinPeakHeight',peak_thr);
    frq=zeros(size(locs'));
    for ind=1:length(locs)
        x1=freqn(locs(ind)-1:locs(ind)+1);
        y1=Ynormn(locs(ind)-1:locs(ind)+1)';
        frq(ind)=quad_pk(x1,y1); %fit quadratic polynomial through max peak and
        %adjacent frequencies and get max of that quadratic
    end
    if ~isempty(locs) %if there are detections (frequency peaks)
        Zset{n} = frq;
    else
        Zset{n} = [];
    end
end


end

function [peak]=quad_pk(x,y)

x=x(:);
y=y(:);

L=length(y);

[q1,q2]=max(y); %q1=max value; q2=index

if (q2==L) 
	peak=q1;
	return;
end

if (q2==1)
	peak=q1;
	return;
end

qstart=max([q2-1,1]);
qstop=min([q2+1,L]);
args=[qstart:qstop];

y1=y( args );
x1=x( args );

c=quadfit(x1,y1);

peak=-c(2)/(2*c(1));

end

function c=quadfit(x,y)

X=[x.^2 x ones(size(x))];
c=X\y;
end