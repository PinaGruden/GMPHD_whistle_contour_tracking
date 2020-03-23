%% GMPHD DETECTOR
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This script runs a GMPHD detector to track frequency contours of whistles.
% The output of the detector is stored in a struct called DT which has 
% 3 fields for each detected whistle contour: freq (frequency), time, label
% (a unique ID of a given contour).

% Pina Gruden, Institute of Sound and Vibration Research (ISVR), 2014-2018

%-------------------------------------------------------------------------
% For detailed explanation see: Gruden, P. and White, P. (2016). Automated
% tracking of dolphin whistles using Gaussian mixture probability
% hypothesis density filters. Journal of the Acoustical Society of America,
% 140(3),1981-1991.
%-------------------------------------------------------------------------

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear, close all

%% ////////////// Specify PARAMETERS //////////////

% ------------ Parameters for obtaining the measurements ------------------
win_width=2048; %specify hanning window length (samples)
overlap=0.5; %percent overlap 
slide_incr= round((1-overlap)*win_width); %advancment of the window (samples)
freqrange = [2000,50000]; % lower and higher frequency range limits in
% Hz in which your signals occur 
peak_thr =8; %threshold in dB for finding spectral peaks 
%(these become your measurements that will go to the detector)

% ---------  PARAMETERS for the GMPHD detector ----------------------
nClutter=10; %number of clutter points per time step
parameters.wth=0.009; %threshold used in state estimation
parameters.pdet = 0.85; %probability of detection (used in Update calculation)
parameters.psurv = 0.994;%probability of survival
parameters.clutter = nClutter/(50000-2000); % Clutter intensity
parameters.Jmax = 100; %max number of Gaussian componenets used in pruning
parameters.U = 10; %pruning threshold used in Pruning & merging
parameters.Tr = 0.001; %truncation threshold used in Pruning & merging 
load('birthpdf.mat');
models.birthpdf=birthpdf;


tl=10; %specify minimal track length in time steps - this is how long a 
%detection needs to be in order to become a true detection- tl=10 gives you
%tl*dt = 0.053 s (dt=time increment of your spectrogram). Using shorter tl
%will give you more false alarms, but potentially better recall (more
%detected sounds), and using higer tl will give you high precision (less
%false alrams) but consequently worse recall (less sounds that you expected
%to retrieve).


%% /////////////// RUN GMPHD DETECTOR ///////////////////


% ~~~~~~~~~~~~~~~~~~~~~~ Read audio data ~~~~~~~~~~~~~~~~~~~~~~~~~~~

[x,fs]=audioread('TestFile.wav');

%Note: the following only handles one channel data. Adjust as necessary.

%% ~~~~~~~~~~~  Initialize remaining PARAMETERS ~~~~~~~~~~~~~~~~
%--------------Parameters for System and observation Models --------------
% System model: xk =F*xk_1 + v; 
% Measurement model: zk = H*xk + w;

dt=slide_incr/fs;%time increment
models.F = [1, dt; 0, 1]; %state transition matrix
models.Q =[5013,166770;166770, 53973000]; %system noise covariance matrix
models.H = [1, 0]; % measurement matrix
models.R= round((fs/win_width)^2/12) ;%measurement noise covariance matrix 
  
       
%% ~~~~~~~~~~~~~ Make MEASUREMENTS (Zsets) ~~~~~~~~~~~~~~~~~~~
[Zset] = preprocess_getZset(win_width,overlap,fs,x,freqrange,peak_thr);
%Zset is a cell array where each cell contains spectral peaks (frequencies)
% from a particular time step that were above threshold (peak_thr - in dB).


%% ~~~~~~~~~~~~~~~~~ Run GMPHD detection ~~~~~~~~~~~~~~~~~~~

[Xk_m,XTag]=gmphd_freqonly_adaptive(Zset,parameters,models);

%% ~~~~~~~~~~~~~~~~ Extract Tracks (whistle contours) ~~~~~~~~~~~~~~~~~~

[DT,~,~] = tracktarget(XTag,Xk_m,dt,tl);
%DT is a structure of detected signals with 3 fields for each (freq x time x label)

%save([folder,'GMPHD_',file(1:end-4),'.mat'],'DT')

%% ~~~~~~~~~~~~~~~~~~ PLOT Detections ~~~~~~~~~~~~~~~~~~~~~~~~~
figure
t=(0:(size(Zset,2)-1)).*dt;
for m=1:size(Zset,2)
plot(t(m),Zset{m},'k.'),hold on
end
for m=1:size(DT,2)
plot(DT(m).time,DT(m).freq,'LineWidth',1.5),hold on
end

