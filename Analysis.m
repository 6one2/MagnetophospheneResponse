close all

savedata=false;

Fs=10000;

% Protocol = importdata ('/Users/Cadence/OneDrive/PhospheneProject/DataCollection/P48/Cadence_P48.txt');

%% EEG File
EEG.Properties.VariableNames{4} = 'O2';
EEG.Properties.VariableNames{5} = 'O1';
EEG.Properties.VariableNames{6} = 'Oz';
EEG.Properties.VariableNames{7} = 'MF';
EEG.Properties.VariableNames{8} = 'Laser';

EEG_MF = EEG.MF; 
EEG_MF = EEG_MF / 1000000; % microvolts to volts
EEG_MF = EEG_MF * 20; % probe voltage to mT conversion
EEG_MF = EEG_MF * 1.2; % distance conversion
EEG_MF = EEG_MF * 2; % probe connected to one source 
EEG_MF = EEG_MF / sqrt(2);
EEG_Time=linspace((1/Fs),(length(EEG_MF)/Fs),(length(EEG_MF)));

figure(), plot(EEG_Time, EEG_MF)
xlabel ('Time (s)', 'FontSize', 20) 
ylabel('MF Generated (mT)', 'FontSize', 20)
hold on
plot(EEG_Time, EEG.Laser/100000)
set(gca,'FontSize',20);
set(gcf,'color','w');

[M,f]=FFTspectrum(EEG_MF, Fs);
figure(), plot(f, M)
xlim([0 305])
xlabel 'Frequency (Hz)'
set(gca,'FontSize',20);
set(gcf,'color','w');

EEG_Laser = EEG.Laser/100;
figure(), plot(EEG_Laser) % plotting exposure square wave of exposures
idxField = find(EEG_Laser>9800); % Determine when the square wave is greater than 0.9
hold on
idxEnd = find(diff(idxField)>2.5E5); % Diff finds the difference between adjacent numbers greater than 0.9, find is when the diff is greater than 2.5E5
plot(idxField,EEG_Laser(idxField),'k.') % Plot when the square wave is greater than 0.9
plot(idxField(idxEnd),EEG_Laser(idxField(idxEnd)),'ro', 'MarkerSize', 10) % Identify the end of the square wave
plot(idxField(idxEnd+1),EEG_Laser(idxField(idxEnd+1)),'bo', 'MarkerSize', 10) % Identify the start of the square wave
plot(EEG.O2, 'g') 
title 'Raw O2 Electrode'
set(gca,'FontSize',20);
set(gcf,'color','w');


%% Filter
% SELECT DATA
% Exposure_End = idxField(idxEnd(17)) % put in whatever number exposure I want the end of 

Section = ((Exposure_End-(5*Fs)):Exposure_End);
DataSource = EEG.O1((Section),1);

% GLOBAL PARAMETERS
window = 8092*4;
overlap = 4096*4;
Fs = 10000;
alpha_begin_idx = ceil(((window/2)/(Fs/2))*8)+1; %14; % index at which 8 Hz begins, 34 with window =4096 and overlap=2048; 17 for window=2048; 28 for beta
alpha_end_idx = ceil(((window/2)/(Fs/2))*12)+1;

Hifilt = 3;
Lowfilt = 17; 
ramp1 = 0;
ramp2 = 0;
    
current_electrode_raw = DataSource-mean(DataSource);
current_electrode_raw = current_electrode_raw (1:50000);

window2 = hann (length (current_electrode_raw));
current_electrode_win = current_electrode_raw.*window2;

figure (), plot (current_electrode_raw,'r');
hold on
plot (smooth(current_electrode_raw,10),'b');


% FILTER
% HIGHPASS and LOWPASS FILTERING - raw section 
hp = 1-makefilt(Hifilt/Fs,length(current_electrode_raw),ramp1/Fs);  % Filtre les données passe haut
var_four_raw = datafilt (current_electrode_raw, hp); 

lf = makefilt (Lowfilt/Fs,length(var_four_raw), ramp2/Fs);   %  filtre passe  bas
var_four_raw = datafilt (var_four_raw,lf);

% HIGHPASS and LOWPASS FILTERING - windowed section  
hp = 1-makefilt(Hifilt/Fs,length(current_electrode_win),ramp1/Fs);  % Filtre les données passe haut
var_four_win = datafilt (current_electrode_win, hp); 

lf = makefilt (Lowfilt/Fs,length(var_four_win), ramp2/Fs);   %  filtre passe  bas
var_four_win = datafilt (var_four_win,lf);


% PLOTS
figure
subplot (3,1,1) 
plot (current_electrode_raw/max(current_electrode_raw)),
hold on 
plot (window2, 'r', 'LineWidth', 4), 
plot (current_electrode_win/max(current_electrode_win), 'g')
ylim ([-1 1])
xlim ([0 5E4])
set(gca,'FontSize',20);
set(gcf,'color','w');
% title ('Oz Electrode', 'FontSize', 20)
ylabel 'Raw Electrode'
xlabel 'Time (s)'
set(gca,'fontsize',40)
set(gcf,'color','w');
set(gca, 'FontName', 'Times New Roman')

subplot (3,1,2)
plot (var_four_raw, 'LineWidth', 4)
hold on
plot (var_four_win, 'g', 'LineWidth', 4)
set(gca,'FontSize',20);
set(gcf,'color','w');
ylim ([-45 45])
xlim ([0 5E4])
ylabel 'Filtered Electrode'
xlabel 'Time (s)'
% legend ('Filtered Data', 'Windowed Filtered Data')
set(gca,'fontsize',40)
set(gcf,'color','w');
set(gca, 'FontName', 'Times New Roman')

[spec_raw,f] = pwelch(var_four_raw,window,overlap,window,Fs);
[spec_win,f] = pwelch(var_four_win,window,overlap,window,Fs);

subplot (3,1,3)
plot (f, spec_raw, 'LineWidth', 4)
hold on,
plot (f, spec_win, 'g', 'LineWidth', 4)
xlim ([0 20])
ylim ([0 7])
% legend ('Filtered Data', 'Windowed Filtered Data')
xlabel 'Frequency (Hz)'
set(gca,'FontSize',20);
set(gcf,'color','w');
ylabel 'Alpha Power'
set(gca,'fontsize',40)
set(gcf,'color','w');
set(gca, 'FontName', 'Times New Roman')

%% Alpha Power

left_bound = 8; %Chose the left frenquency bound
right_bound =12; %Choose the right frequency bound

%find when the absolute difference between f and left_bound is minimum (i.e
%when f is closest to 8 and store that location in I_left
[M,I_left] = min(abs(f(:,1)-left_bound)); 

%find when the absolute difference between f and right_bound is minimum (i.e
%when f is closest to 8 and store that location in I_right

[M,I_right] = min(abs(f(:,1)-right_bound));

a(:,1) = f(I_left:I_right);
a(:,2) = spec_win(I_left:I_right);

figure
plot(a(:,1),a(:,2))

Int=trapz(a(:,2))


%% labview file
Labview.Properties.VariableNames{1} = 'Time';
Labview.Properties.VariableNames{2} = 'MF';
Labview.Properties.VariableNames{3} = 'Dial';
Labview.Properties.VariableNames{4} = 'Current';
Labview.Properties.VariableNames{5} = 'Command';

Labview_MF = Labview.MF / sqrt(2); % Converts probe MF data to rms 

time = Labview.Time; %
idxGood = ~isnan(time);
Time = time(idxGood); % time for graph but not used for finding start

idxSLab = find(time==0);
idxSLab = idxSLab(2:51); % Should be 51 in a full experiment

figure()
plot(Labview_MF);hold on % Plot Probe in rms values
plot(idxSLab-1,Labview_MF(idxSLab-1),'r.','MarkerSize',20) %Finds the end of the exposure
ylabel 'MF (mT)'
set(gca,'FontSize',20);
set(gcf,'color','w');

figure(), plot(Labview.Dial);hold on % Plot dial voltage peak
plot(idxSLab-1,Labview.Dial(idxSLab-1),'r.','MarkerSize',20) % Find end of exposure
ylabel 'Magnetophosphene Threhsold (mT) - peak'
set(gca,'FontSize',20);
set(gcf,'color','w');
title 'Dial'

%% Make Threshold File
for i = 1:50 % Should say 50 in a full experiment
    row_num = idxSLab(i,1);
    dial(i,1) = Labview(row_num-1,3);
end

Dial = dial{:,1:end}/(sqrt(2));
Threshold_mT = Dial.*18.04; % mT rms values of the MF thresholds based on the voltage sent to the amplifier (from the dial)

figure()
plot(Threshold_mT)
xlabel 'Iteration Number'
ylabel('Magnetophosphene Threshold (mT_{rms})', 'FontSize', 20)

%% Make Threshold Frequency Response Graph

% Load the thresholds60.txt file as column vectors -- values in this file
% are peak threshold values


err=SEM;
% P = Percentage.*100;

figure
bar(Frequency, Percentage,'k','FaceAlpha',1.0)
ylabel 'Magnetophosphene Perceived (%)'
hold on
%yyaxis left, plot(Frequency (1:25), Average (1:25)/sqrt(2), '.k','MarkerSize',20, 'LineWidth', 4)
hold on
%errorbar(Frequency(1:25), Average (1:25)/sqrt(2), err (1:25),'vertical','k', 'LineWidth', 4)
xlim ([0 303])
%ylim([0 70])
%yticks([0 10 20 30 40 50 60 70])
xticks([10 20 30 40 50 60 70 80 90 100 150 200 250 300])
xlabel ('Frequency (Hz)') 
%ylabel ({'Magnetophosphene Threhsold (mT_{rms})'}, 'Color','k')
set(gca,'fontsize',35)
set(gcf,'color','w');

figure
bar(Percentage(:,1))
xlim([0 26])
set(gca,'xtick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25])
set(gcf,'color','w');
xticklabels({'0', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70', '75', '80', '85', '90', '95', '100', '150', '200', '250', '300'});
xlabel ('Frequency (Hz)', 'FontSize', 30)
ylabel ({'Magnetophosphenes', 'Perceived (%)'}, 'FontSize', 30)

figure
plot(Frequency (1:25), Average (1:25)/sqrt(2), '.k','MarkerSize',20)
hold on
errorbar(Frequency(1:25), Average (1:25)/sqrt(2), err (1:25)/sqrt(2),'vertical','k')
xlim ([0 95])
ylim([0 70])
yticks([0 10 20 30 40 50 60 70 80 90 100])
xticks([10 20 30 40 50 60 70 80 90 100 150 200 250 300])
xlabel ('Frequency (Hz)') 
ylabel ('Magnetophosphene Threhsold (mT_{rms})')
%title ('Frequency Response of Magnetophosphene Perception', 'FontSize', 16, 'FontWeight', 'bold')
set(gca,'fontsize',30)
set(gcf,'color','w');

% Load short.txt file
figure, plot(Mean/sqrt(2), '.k','MarkerSize',20)
hold on
errorbar(Mean (1:4)/sqrt(2), SEM1, 'vertical','k')
ylim([0 70])
xlim ([0.9 4.1])
xticks([1 2 3 4])
xticklabels({'Trial 1','Trial 2','Trial 3','Trial 4'})
xlabel ('20 Hz Estimation') 
ylabel ({'Magnetophosphene', 'Threhsold (mT_{rms})'})
set(gca,'fontsize',30)
set(gcf,'color','w');

%dbdt graph

figure
plot(Frequency(1:9,:),dbdt(1:9,:)/sqrt(2), '.b', 'MarkerSize',20)
hold on
errorbar(Frequency(1:9,:), dbdt(1:9,:)/sqrt(2), SEMdbdt(1:9,:)/sqrt(2),'vertical','b', 'LineWidth', 4)
xlim ([0 95.2])
xticks([10 20 30 40 50 60 70 80 90 100 150 200 250 300])
xlabel ('Frequency (Hz)') 
ylabel ({'Magnetophosphene Threhsold (dB/dt)'})
set(gca,'fontsize',22)
set(gcf,'color','w');

figure
plot(Frequency(9:19,:),dbdt(9:19,:)/sqrt(2), '.r', 'MarkerSize',20)
hold on
errorbar(Frequency(9:19,:), dbdt(9:19,:)/sqrt(2), SEMdbdt(9:19,:)/sqrt(2),'vertical','r', 'LineWidth', 4)
xlim ([0 95.2])
xticks([10 20 30 40 50 60 70 80 90 100 150 200 250 300])
xlabel ('Frequency (Hz)') 
ylabel ({'Magnetophosphene Threhsold (dB/dt)'})
set(gca,'fontsize',22)
set(gcf,'color','w');

figure
plot(Frequency(1:20,:),dbdt(1:20,:)/sqrt(2), '.k', 'MarkerSize',20)
hold on
errorbar(Frequency(1:20,:), dbdt(1:20,:)/sqrt(2), SEMdbdt(1:20,:)/sqrt(2),'vertical','k', 'LineWidth', 4)
xlim ([0 95.2])
%xticks([10 20 30 40 50 60 70 80 90 100 150 200 250 300])
xlabel ('Frequency (Hz)') 
ylabel ({'Magnetophosphene Threhsold (dB/dt)'})
set(gca,'fontsize',22)
set(gcf,'color','w');

%% This script splits the Labview data into 25 second chunks by exposure

readfolder = '/Users/Cadence/Documents/OneDrive/PhospheneProject/DataCollection/P14/';
fname = 'Labview';

% resdata = importdata([readfolder fname],'\t',23);

Fs=10000;
RestT = 35;

time = Labview(:,2);
ExpTime = nan(size(time));
SecTime = nan(size(time));
idxG = ~isnan(time);

idxS = find(time==0);
idxE = idxS(2:end)-1;idxE=[idxE;length(time)];

for k=1:size(idxS,1)
    duration = (idxE(k)-idxS(k))/Fs;
    SecTime(idxS(k):idxE(k)) = 0:1/Fs:duration; %% time from 0:25s for each section
    
    if k==1
        ExpTime(idxS(k):idxE(k)) = idxS(k)/Fs:1/Fs:idxE(k)/Fs; %% time in real time from beginnig to endq
    else
        timeShift = idxE(k-1)/Fs+RestT;
        ExpTime(idxS(k):idxE(k)) = idxS(k)/Fs+timeShift:1/Fs:idxE(k)/Fs+timeShift;
    end
end

MF = resdata.data(:,3);
Dial = resdata.data(:,4);
Current = resdata.data(:,5);
Command = resdata.data(:,6);


%% Graph a single exposure from the labview file 
trial = 5;% choose the section (iteration#)
section = idxS(trial):idxE(trial);

figure()
title ('Trial # '+string(trial))
subplot(211)
plot(SecTime(section), MF(section))
xlim([SecTime(section(1)),SecTime(section(end))])
xlabel('Time (s)')
legend 'MF'
subplot(212)
plot(ExpTime(section)/60, MF(section))
xlim([ExpTime(section(1))/60,ExpTime(section(end))/60])
xlabel('Time (min)')
legend 'MF'

%% Plot all variables from the Labview file
figure()
title ('Trial # '+string(trial))
subplot(411)
plot(SecTime(section), Dial(section))
xlim([SecTime(section(1)),SecTime(section(end))])
xlabel('Time (s)')
legend 'Dial'

subplot(412)
plot(SecTime(section), Command(section))
xlim([SecTime(section(1)),SecTime(section(end))])
xlabel('Time (s)')
legend 'Command'

subplot(413)
plot(SecTime(section), Current(section))
xlim([SecTime(section(1)),SecTime(section(end))])
xlabel('Time (s)')
legend 'Current'

subplot(414)
plot(SecTime(section), MF(section))
xlim([SecTime(section(1)),SecTime(section(end))])
xlabel('Time (s)')
legend 'MF'

%% Make a Graph of the Dial
Section = Labview.Dial(4.79E6:5.0274E6); %P44
Section = Section.*18.04;
Section = Section./sqrt(2);
Time = linspace(0,31,length(Section));
plot(Time,Section, 'LineWidth', 4,'Color', 'k')
xlim([0 30.5])
ylim([0 40])
set(gca,'fontsize',40)
set(gcf,'color','w');
xlabel 'Time (s)'
ylabel ({'Magnetophosphene', 'Threhsold (mT_{rms})'})
set(gca, 'FontName', 'Times New Roman')

%% EEG plots 
subplot (3,1,1)
plot (var_four_win, 'g')
set(gca,'FontSize',20);
set(gcf,'color','w');
ylabel 'Oz Electrode'
legend '35 Hz Exposure', 'Sham (0 Hz) Exposure'

[spec_raw,f] = pwelch(var_four_raw,window,overlap,window,Fs);
[spec_win,f] = pwelch(var_four_win,window,overlap,window,Fs);

subplot (2,1,2)
plot (f, spec_win, 'g')
xlim ([0 25])
xlabel 'Frequency (Hz)'
set(gca,'FontSize',20);
set(gcf,'color','w');
ylabel 'Alpha Power'
legend '35 Hz Exposure', 'Sham (0 Hz) Exposure'

