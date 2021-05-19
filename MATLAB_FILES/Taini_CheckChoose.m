%%%%%%TAINI CHANNELS CHECK & CHOOSE CHANNELS FOR SLEEP SCORING
%%% OPTINAL SWDs auto-detection included
%%% Ingrid Buller UoE 2021- ingrid.buller@ed.ac.uk
%%%Taini_Check&Choose

clearvars -except d;

day= 'BL2'; %rat_index(:,6)
rat='1959'; %rat_index(:,3)
line='CDKL5'; %rat_index(:,4)
ratname=[line '_' rat];
rat_ind=['rat_day_index_' line '.csv']
pn='C:\Users\ibull\Documents\POSTDOC\ANALISYS\EEG RESULTS\For Igor'; % Set dir for index data .csv file
cd(pn);
T = readtable(rat_ind);
rat_index=table2cell(T);
i = strcmp(rat_index(:,2), [ratname '-' day]);
i2=double(i);
k=find(i2(:,1) == 1); %row index
p='E:\TAINI DATA'; %Main directory source for .dat files
pn = [p '\' line '\' ratname]
dn= char(rat_index(k,8));

Num_index = table2array(T(:,9:15));
if (k>=2)
l = strcmp(rat_index(k-1,8),dn);
m = l;
end

x=1;
while x==1
	if (k==1) || (m==0)
       	  dat_filename = [pn '\' dn];
         [d, t] = parse_dat(dat_filename);
         %save('d.mat','d');
        else
         x=0
        end
      x=exist('d')-1;  
end

x=1;
while x==1
if (day=='BL3') & (Num_index(k,3)>=1)
      startA = Num_index(k,1);
      endA = Num_index(k,2);
      c1A=d([startA:endA],4);
      c2A=d([startA:endA],13);
      c3A=d([startA:endA],2);
      c4A=d([startA:endA],15);
      startB = Num_index(k,3);
      endB = Num_index(k,4);
      c1B=d([startB:endB],4);
      c2B=d([startB:endB],13);
      c3B=d([startB:endB],2);
      c4B=d([startB:endB],15);
      cn1=[c1A;c1B];
      clear c1A c1B;
      cn2=[c2A;c2B];
      clear c2A c2B;
      cn3=[c3A;c3B];
      clear c3A c3B;
      cn4=[c4A;c4B];
      clear c4A c4B;
   else
      samp_start=Num_index(k,1);
      samp_end=Num_index(k,2);
      cn1=d([samp_start:samp_end],4);
      cn2=d([samp_start:samp_end],13);
      cn3=d([samp_start:samp_end],2);
      cn4=d([samp_start:samp_end],15);
end
      x=exist('cn4')-1;
end

%Remove DC components
x=1;
while x==1
 dcoff_1=cn1-mean(cn1);
 dcoff_2=cn2-mean(cn2);
 dcoff_3=cn3-mean(cn3);
 dcoff_4=cn4-mean(cn4);
 x=exist('dcoff_4')-1;
end

%Remove outliers on signal and filter
x=1;
while x==1
 outlier_off1 = filloutliers(dcoff_1,'linear','ThresholdFactor',4.5);
 filt1 = highpass(outlier_off1,0.02,'Steepness',0.5,'StopbandAttenuation',120);
 outlier_off2 = filloutliers(dcoff_2,'linear','ThresholdFactor',4.5);
 filt2 = highpass(outlier_off2,0.02,'Steepness',0.5,'StopbandAttenuation',120);
 outlier_off3 = filloutliers(dcoff_3,'linear','ThresholdFactor',4.5);
 filt3 = bandpass(outlier_off3,[0.48 0.72],'Steepness',0.5,'StopbandAttenuation',60);
 outlier_off4 = filloutliers(dcoff_4,'linear','ThresholdFactor',4.5);
 filt4 = bandpass(outlier_off4,[0.48 0.72],'Steepness',0.5,'StopbandAttenuation',60);
 x=exist('filt4')-1;
end

%PLOT & CHOOSE---------------------
% Plot EEGs to compare
x=1;
while x==1
figure(1);
set(gcf,'position',[51.4 456.2 1408 252.8])
subplot(2,1,1);
plot(filt1)
ylabel ('eeg1')

subplot(2,1,2);
plot(filt2)
ylabel ('eeg2')

% Choose EEG channel
prompt = {'1 or 2 and press Enter'};
dlgtitle = 'Choose EEG channel';
dims = [1 35];
definput = {'1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
choice = cell2mat(answer);
x= str2double(choice);
disp(x);

%%% Set EEG channel
if x == 1
	ch1=cn1;
	eeg=4;
    disp('eeg1 chosen');
    close
      else
	ch1=cn2;
	eeg=13;
    disp('eeg2 chosen');
    close
end
Num_index(k,6)=eeg;
x=exist('eeg')-1;
end

% Plot EMGs to compare
x=1;
while x==1
figure(1);
set(gcf,'position',[51.4 456.2 1408 252.8])
subplot(2,1,1);
plot(filt3)
ylabel ('emg1')

subplot(2,1,2);
plot(filt4)
ylabel ('emg2')

%%% Choose EMG channel
prompt = {'1 or 2 and press Enter'};
dlgtitle = 'Choose EMG channel';
dims = [1 35];
definput = {'1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
choice = cell2mat(answer);
x= str2double(choice);
disp(x);

%%% Set EMG channel
if x == 1
    ch2=cn3;
    emg=2;
    disp('emg1 chosen');
    close
      else
    ch2=cn4;      
    emg=15;
    disp('emg2 chosen');
    close
end
Num_index(k,7)=emg;
x=exist('emg')-1;
end

pna=['C:\Users\ibull\Documents\POSTDOC\ANALISYS\EEG RESULTS\For Igor\' line '\' ratname];
if ~exist(pna, 'dir')
       mkdir(pna)
end
Chanls = [ch1 ch2];
Channels = array2table(Chanls);
cd(pna);
outputdir=pwd;  
writetable(Channels,[ratname '_' day '_Channels.csv'])    

%%%% CHECK SWDs?------------------------------
done=1;
x=0;
while done==1 
	choice = menu('Check SWD?','Yes','No');
	if choice==1 | choice==0
		x=1;
	        done=0;
        else
	        x=0;
	        close
	        done=0;
                break
        end
done=0;
end

if x==1
   cn1=ch1; %Chosen eeg ch for SWD detection
   TAINI_autodetect(cn1,pn,pna,day,ratname,line,T,k,rat_index,Num_index)
end
