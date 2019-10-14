% script DeltaFs

% a script for photometry signal processing

% to compute dF/F
% 1. linear least fit 405nm signals to 470nm signals
% 2. F470 - fitted F405 = dF
% 3. compute median value as baseline
% 4. compute dF - baseline
%


clearvars;
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for dat
nChs = 10; % the number of channels
Fs = 1000; % sampling rate
Ofs = 40; % optical stim frequency

maxVolts = 5;
sampsPerVolt = double(intmax('int16')/maxVolts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
Fname = 'Raw.dat'; % raw data file
fid = fopen(Fname);
Data = fread(fid, [nChs, Inf], 'int16')/sampsPerVolt;
fclose(fid);

%% channel configration
Sync405 = Data(1,:); % 405 nm sync pulses
Sync470 = Data(2,:); % 470 nm sync pulses
Fsig = Data(3,:); % photodector signals

clear Data;

T = [1:length(Sync405)]/Fs; % in s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% signal processing for 405 nm signals
[~, EvtPoints, EvtAmps] = MyEvtDetection(Sync405, Fs);
[MyFs405, MyTs] = MyFextraction(Fsig, EvtPoints, EvtAmps, T);

%% Signal processing for 470nm signals 
[~, EvtPoints, EvtAmps] = MyEvtDetection(Sync470, Fs);
[MyFs470, MyTs] = MyFextraction(Fsig, EvtPoints, EvtAmps, T);

%% signal alignment
if length(MyFs405) ~= length(MyFs470)
    MySize = min([length(MyFs405), length(MyFs470)]);
    MyFs405 = MyFs405(1:MySize);
    MyFs470 = MyFs470(1:MySize);
    MyTs = MyTs(1:MySize);
end

clear Fsig T Sync405 Sync470 EvtPoints EvtAmps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dF (power fit) -- CORE
% filter setting
LowCutFreq = 4; 
[b,a] = butter(5, LowCutFreq/(Ofs/2),'low');

%% power2 fit for 405 nm signals
MyF = fit(MyTs', MyFs405, 'power2'); % fit to a*x^b+c
MyCoeff = coeffvalues(MyF); % MyF is a cfit object variable
ExpFit405 = MyCoeff(1)*MyTs.^MyCoeff(2)+MyCoeff(3); % MyCoeff ... a, b, c
MyFs405norm = MyFs405 - ExpFit405'; % subtracted photobleaching component
% + filtering
MyFs405norm = filtfilt(b, a, MyFs405norm);

%% power2 fit for 470 nm signals
MyF = fit(MyTs', MyFs470, 'power2'); % fit to a*x^b+c
MyCoeff = coeffvalues(MyF); % MyF is a cfit object variable
ExpFit470 = MyCoeff(1)*MyTs.^MyCoeff(2)+MyCoeff(3); % MyCoeff ... a, b, c
MyFs470norm = MyFs470 - ExpFit470'; % subtracted photobleaching component
% + filtering
MyFs470norm = filtfilt(b, a, MyFs470norm);

%% linear fit
x = MyFs405norm\MyFs470norm;
MyFs405_Fit = x*MyFs405norm;
% filtering (a bit strong filter)
LowCutFreq = 1; 
[b,a] = butter(12, LowCutFreq/(Ofs/2),'low');
MyFs405_Fit = filtfilt(b, a, MyFs405_Fit);

%% final signals
dFs = MyFs470norm - MyFs405_Fit;

%% save signals
save('Signals.mat','dFs','MyTs','MyFs405norm','MyFs470norm', 'Ofs');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% event detection
function [EvtTimes, EvtPoints, EvtAmps] = MyEvtDetection(MySync, Fs)

%% params
HamWin = 8;
LookAhead = 4*(Fs/1000);
syncThresh = 0.15;         % exclude events on sync chan (trigCh) less than 0.1V
syncRefrac = 4*(Fs/1000);

%% filter
b = diff(hamming(HamWin),1); 
b = b-mean(b);   % diff of smoothed sig
dEeg=Filter0(b, MySync);                          % find event times & amplitudes
minima = LocalMinima(-abs(dEeg), syncRefrac, -syncThresh)-1; % To find LocalMinima, see /ken's code/General/

%% outcomes
EvtPoints = minima((minima>LookAhead & minima<length(dEeg)-LookAhead));
EvtTimes=EvtPoints/(Fs/1000);  % ms-based
EvtAmps = MySync(EvtPoints + LookAhead) - MySync(EvtPoints-LookAhead);

if rem(length(EvtPoints),2) ~= 0 % make sure whether all pulses are complete
    EvtPoints(end) = [];
    EvtTimes(end) = [];
    EvtAmps(end) = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract Fs
function [MyFs, MyTs] = MyFextraction(OrgFs, EvtPoints, EvtAmps, T)

%% prep
nPls = length(EvtPoints)/2;
MyFs = zeros(nPls,1);

%% time points
ON = EvtPoints(EvtAmps > 0);
OFF = EvtPoints(EvtAmps < 0);

for p = 1:nPls
    MyFs(p) = median(OrgFs(ON(p):OFF(p)));
end

MyTs = T(ON);

end
