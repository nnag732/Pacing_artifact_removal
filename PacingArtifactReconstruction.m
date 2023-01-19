%% SCRIPT INFORMATION
% Author: Nipuni D. Nagahawatte
% Affiliation: Auckland Bioengineering Institute, University of Auckland, New Zealand
% Date: 18.01.2023
% Email: nnag732@aucklanduni.ac.nz

% With reference to the journal paper:
% Title: A novel framework for the removal of pacing artifacts from bio-electrical recordings
% Authors: Nipuni D. Nagahawatte; Niranchan Paskaranandavadivel; Laura R. Bear; Recep Avci; Leo K Cheng
% Journal: Computers in Biology and Medicine

%Cite As:
%N. D. Nagahawatte, N. Paskaranandavadivel, L. R. Bear, R. Avci, and L. K. Cheng, 
%“A novel framework for the removal of pacing artifacts from bio-electrical recordings,” 
%Comput. Biol. Med., vol. Under revi, 2023.

% The script creates a sample data matrix using the sampleData.mat and add
% pacing artifacts
% The added artifacts are detected using a Hampel filter and reconstructed
% using autoregression, weighted means and linear interpolation
% The signals reconstructed using the 3 methods are visualized on top of
% the original trace and displayed at the end
%% Script to detect pacing artifacts using Hampel filter and reconstruct using Autoregression, Weighted means, and Linear interpolation
clear;
clc;
close all;

%test for communications toolbox
TB=contains(struct2array(ver), 'Communications Toolbox');
if TB==0
    disp(' Please install the Communications Toolbox')
    keyboard
else
end

%creating synthetic slow waves with pacing artefacts

%call for synthetic slow wave function
[sigArray,t_Synth,X,Y,T] = CreateSynthMaps();

%create a new matrix of 100 slow wave channels
syn_array=sigArray.*0.15;

%define the time axis
t=[1/30:1/30:length(sigArray(20,:))/30]; 
syn_array=syn_array;

%matrix with stim off
Soff=syn_array; 
Snoisy=awgn(Soff,14,'measured');

for i=1:size(Snoisy,1)
    Snoisy(i,:)=Snoisy(i,:)+sin(2*pi*t*0.2)*15;    
end
Soff=Snoisy; syn_array=Snoisy; %COMMENT TO REMOVE NOISE

%introduce synthetic pacing artefacts
pulsewidth = 0.2;
fs=30;
dt=1/fs;
pulseperiods(1,:) = [50:10:150];
pulseperiods(2,:) = [200:10:300];
pulseperiods(3,:) = [350:10:450];
pulseperiods(4,:) = [500:10:600];

amp=[2000 2000 2000 2000];
amp=amp*2;

for i=1:size(syn_array,1)
        stim(i,:) = amp(1)*pulstran(t,pulseperiods(1,:),@rectpuls,pulsewidth);
        stim(i,:) = stim(i,:)+amp(2)*pulstran(t,pulseperiods(2,:),@rectpuls,pulsewidth);
        stim(i,:) = stim(i,:)+amp(3)*pulstran(t,pulseperiods(3,:),@rectpuls,pulsewidth);
        stim(i,:) = stim(i,:)+amp(4)*pulstran(t,pulseperiods(4,:),@rectpuls,pulsewidth);
end

pulseperiods = [210:20:510];
amp=4e3;
for i=1:size(syn_array,1)
    stim(i,:) = amp*pulstran(t,pulseperiods,@rectpuls,pulsewidth);
    i=i+1;
end

% stim=repmat(stim,10,1);
stim(end:numel(syn_array))=0;

%combine slow wave and artefact
syn_array=syn_array+stim;
Son=syn_array; %matrix with stim on

%% detecting artefacts using hampel outliers (stim ON)
[A1, postStim_FG, stFG,endFG]=artRem_linreg(Son,75,fs,fs,0);  %inputs=dataset,k,sigma,samplingfreq(downsampled),samplingfreq (original-acquired),type (0=synthetic, 1=experiment)
[A1, postStim_WM, stWM,endWM]=artRem_weightmean3(Son,75,fs,fs,0); %sigma was 10
[A1, postStim_LI, stLI, endLI]=artRem_linint(Son,75,fs,fs,0);

t=[1/30:1/30:length(sigArray(20,:))/30];
figure,plot(t,A1(50,:));hold on;plot(t,(postStim_FG(50,:)));hold on;plot(t,(postStim_WM(50,:)));hold on;plot(t,(postStim_LI(50,:)));
legend('Original Data','Fillgaps method','Weighted mean method','Linear interpolation method',...
       'Location','SW') 
   
%% Autoregression Method
function [A1, artRem, startPos_fs, endPos_fs]=artRem_linreg(Data,k,efs,sr,type)
%Data:Electrode data matrix from the bdf file
%k:k for the hampel filter(number of neighbouring data points)
%sigma:sigma for the hampel filter
%efs:sampling frequency
%sr:native sampling rate of bdf
%type=0 if gastric synthetic,1=GI bdf, 2=cardiac synthetic,3=cardiac exp

%Created on 23.11.2021
%Pacing or stimulation artefacts are identified by the hampel filter and
%are removed using a linear regression approach


%choosing electrode channels
if type==1
    dl=size(Data,1)-1;  %all recorded channels (minus trigger)
else
    dl=size(Data,1);
end

A1 = Data(1:dl,:); 

%time length
sigTest = Data(3,:);
t = 1/sr:1/sr:length(sigTest)/sr;

%remove baseline
if type==1
    for i = 1:size(A1,1)
        baseMed = medfilt1(A1(i,:),(sr*3)+1);
        base =sgolayfilt(baseMed,2,(sr*2)+1);
        A1(i,:)=A1(i,:)-base;
    end
elseif type==3
    for i = 1:size(A1,1)
        baseMed = medfilt1(A1(i,:),round((sr*1.5))+1);
        base =sgolayfilt(baseMed,2,round((sr*0.7))+1);
        A1(i,:)=A1(i,:)-base;
    end
else
end

%replacing the first and last few samples of the recording
if type==2||3 %cardiac
    startInt=round(efs*0.01); 
else %GI
    startInt=round(efs*0.1); 
end
lastSample=size(A1,2);
endInt=lastSample-startInt;
    
for s=1:startInt
    A1(:,s)= A1(:,(startInt+1)); %replacing with the sample immediately after 0.1s
    s=s+1;
end

for e=endInt:lastSample
    A1(:,e)=A1(:,(endInt-1)); %replacing with the sample 0.1s before the end of the recording
    e=e+1;
end

%taking the mean absolute of all signals 
sigAll =abs(A1);
sigAllMean = mean(sigAll,1);

% sigma=1.5*(range(sigAllMean)/std(sigAllMean));
sigma=40;

%Apply Hampel 
[yAll,iAll,xmedianAll,xsigmaAll] = hampel(sigAllMean,round(efs*k),sigma); %k=0.5s, sigma=10


%Detecting outliers
outly=sigAllMean(iAll);
outlt=t(iAll);

%find outliers that are more than 1s apart to find the start and end of stimuli 
endPos=[];
startPos=[];

%position of stimuli in time
if type==2
    endPt = find(diff(outlt)>0.5);
elseif type==3
    endPt = find(diff(outlt)>0.3);
else
    endPt = find(diff(outlt)>1);
end
endPos = [outlt(endPt) outlt(end)];
startPos = [outlt(1) outlt(endPt+1)];

%amplitude of the start and end stimuli points
endY = [outly(endPt) outly(end)];
startY = [outly(1) outly(endPt+1)];

%adding a buffer
% bufferst=0.01;
if type==2%ecg synthetic
    bufferst=0.005;
    bufferend=0.005;
elseif type==3
    bufferst=0.02;
    bufferend=0.1;   %was 0.02
else %GI synthetic
    bufferst=0.1; %was 0.02
    bufferend=0.2; %was 0.2 %longer buffer at the end because of the ringing
end
startPos=startPos-bufferst;
endPos=endPos+bufferend;

%% removing artefact by linear regression (23.11.2021)

%set stim artefact to NaNs
startPos_fs=round(efs*startPos);
endPos_fs=round(efs*endPos);

postStim_FG=A1;

for w=1:size(A1,1) 
    for k=1:length(startPos_fs) 
        postStim_FG(w,startPos_fs(k):endPos_fs(k))=NaN; 
        k=k+1;
    end
    w=w+1;
end

%increase the stim window size
if type==2
    win=0.5*efs;
elseif type==3
    win=efs*0.02;
else
    win=0.5*efs;
end
stimwinS=startPos_fs-win;
stimwinE=endPos_fs+win;
% stimwinE=endPos_fs;

%fillgaps for stimuli window
if stimwinS(1)<0 %if the artefact window at the start is going beyond the start of the signal
    if stimwinE(end)>(length(Data(1,:))) %if the artefact window at the end is 'also' going beyond the end of the signal
        for w=1:size(A1,1) 
            for k=2:length(startPos_fs)-1 
                postStim_FG(w,stimwinS(k):stimwinE(k))=fillgaps(postStim_FG(w,stimwinS(k):stimwinE(k))); %set stim artefact to NaNs
                k=k+1;
            end
            postStim_FG(w,0:stimwinE(1))=fillgaps(postStim_FG(w,0:stimwinE(1)));
            postStim_FG(w,stimwinS(end):(length(Data(w,:))))=fillgaps(postStim_FG(w,stimwinS(end):(length(Data(w,:)))));
            w=w+1;
        end
    else
        for w=1:size(A1,1) 
            for k=2:length(startPos_fs)
                postStim_FG(w,stimwinS(k):stimwinE(k))=fillgaps(postStim_FG(w,stimwinS(k):stimwinE(k))); %set stim artefact to NaNs
                k=k+1;
            end
            w=w+1;
            postStim_FG(w,0:stimwinE(1))=fillgaps(postStim_FG(w,0:stimwinE(1)));
        end
    end
elseif stimwinE(end)>(length(Data(1,:))) %if the artefact window at the end is going beyond the end of the signal
    for w=1:size(A1,1) 
        for k=1:length(startPos_fs)-1 
            postStim_FG(w,stimwinS(k):stimwinE(k))=fillgaps(postStim_FG(w,stimwinS(k):stimwinE(k))); %set stim artefact to NaNs
            k=k+1;
        end
        postStim_FG(w,stimwinS(end):(length(Data(w,:))))=fillgaps(postStim_FG(w,stimwinS(end):(length(Data(w,:)))));
    end
else
    for w=1:size(A1,1) 
        for k=1:length(startPos_fs) 
            postStim_FG(w,stimwinS(k):stimwinE(k))=fillgaps(postStim_FG(w,stimwinS(k):stimwinE(k))); %set stim artefact to NaNs
            k=k+1;
        end
        w=w+1;
    end
end
artRem=postStim_FG;

end





%% Weighted mean Method

function [A1,artRem, startPos_fs, endPos_fs]=artRem_weightmean3(Data,k,efs,sr,type)
%Data:Electrode data matrix from the bdf file
%k:k for the hampel filter(number of neighbouring data points)
%sigma=sigma for the hampel filter
%efs=sampling frequency
%sr:native sampling rate of bdf
%type=0 if gastric synthetic,1=GI bdf, 2=cardiac synthetic,3=cardiac exp

%Created on 23.11.2021
%Pacing or stimulation artefacts are identified by the hampel filter and
%are removed using a weighted mean approach


%choosing electrode channels
if type==1
    dl=size(Data,1)-1;  
else
    dl=size(Data,1);
end

A1 = Data(1:dl,:); %all recorded channels (minus trigger)

%time length
sigTest = Data(3,:);
t = 1/sr:1/sr:length(sigTest)/sr;

% A2=A1;

%remove baseline
if type==1
    for i = 1:size(A1,1)
        baseMed = medfilt1(A1(i,:),(sr*3)+1);
        base =sgolayfilt(baseMed,2,(sr*2)+1);
        A1(i,:)=A1(i,:)-base;
    end
elseif type==3
    for i = 1:size(A1,1)
        baseMed = medfilt1(A1(i,:),round((sr*1.5))+1);
        base =sgolayfilt(baseMed,2,round((sr*0.7))+1);
        A1(i,:)=A1(i,:)-base;
    end
else
end

%replacing the first and last few samples of the recording
if type==2||3 %cardiac
    startInt=round(efs*0.01); 
else %GI
    startInt=round(efs*0.1); 
end %0.1s to eliminate
lastSample=size(A1,2);
endInt=lastSample-startInt;
for s=1:startInt
    A1(:,s)= A1(:,(startInt+1)); %replacing with the sample immediately after 0.1s
    s=s+1;
end

for e=endInt:lastSample
    A1(:,e)=A1(:,(endInt-1)); %replacing with the sample 0.1s before the end of the recording
    e=e+1;
end

%taking the mean absolute of all signals 
sigAll =abs(A1);
sigAllMean = mean(sigAll,1);

% sigma=1.5*(range(sigAllMean)/std(sigAllMean));
sigma=40;

%Apply Hampel 
[yAll,iAll,xmedianAll,xsigmaAll] = hampel(sigAllMean,round(efs*k),sigma); %k=0.5s, sigma=10


%Detecting outliers
outly=sigAllMean(iAll);
outlt=t(iAll);

%find outliers that are more than 1s apart to find the start and end of stimuli 
endPos=[];
startPos=[];

%position of stimuli in time
if type==2
    endPt = find(diff(outlt)>0.5);
elseif type==3
    endPt = find(diff(outlt)>0.3);
else
    endPt = find(diff(outlt)>1);
end
endPos = [outlt(endPt) outlt(end)];
startPos = [outlt(1) outlt(endPt+1)];

%amplitude of the start and end stimuli points
endY = [outly(endPt) outly(end)];
startY = [outly(1) outly(endPt+1)];

%adding a buffer
% bufferst=0.01;
if type==2 %ecg synthetic
    bufferst=0.005;
    bufferend=0.005;
elseif type==3
    bufferst=0.02; %(0.02)
    bufferend=0.02;   %was 0.02 (0.1)
else %GI synthetic
    bufferst=0.1; %was 0.1
    bufferend=0.2; %was 0.2 %longer buffer at the end because of the ringing
end
startPos=startPos-bufferst;
endPos=endPos+bufferend;

%% removing artefact by weighted mean (23.11.2021)

%compute weighted mean
weightedMean=bsxfun(@rdivide, A1, sigAllMean);
weightedMean(isnan(weightedMean)) = 0;

startPos_fs=round(efs*startPos);
endPos_fs=round(efs*endPos);

%subtracting the stimuli from the signal
postStim_WM=A1;

%CASE 2: replace artefact with WM*baseline estimate from Mean sig
postStim_WM=A1;

artEst=sigAllMean;
% i=1;
for i=1:length(startPos_fs)
    artLen=round(length(startPos_fs(i):endPos_fs(i))/2);
    prevS=startPos_fs(i)-artLen;
    prevS=artEst(prevS:startPos_fs(i));
    aftS=endPos_fs(i)+(length(startPos_fs(i):endPos_fs(i))-artLen);
    aftS=artEst(endPos_fs(i):aftS);
    artEst(startPos_fs(i):(startPos_fs(i)+artLen))=prevS;
    artEst((endPos_fs(i)-length(aftS)+1):endPos_fs(i))=aftS;
    artEst(startPos_fs(i):endPos_fs(i))=mean(artEst(startPos_fs(i):endPos_fs(i)));
end

for w=1:size(A1,1) 
    for k=1:length(startPos_fs)
        postStim_WM(w,startPos_fs(k):endPos_fs(k))=weightedMean(w,startPos_fs(k):endPos_fs(k)).*artEst(startPos_fs(k):endPos_fs(k)); %each signal-(weighted mean*avg signal)
        k=k+1;
    end
    w=w+1;
end

artRem=postStim_WM;

end




%% Linear interpolation Method

function [A1, artRem, startPos_fs, endPos_fs]=artRem_linint(Data,k,efs,sr,type)
%Data:Electrode data matrix from the bdf file
%k:k for the hampel filter(number of neighbouring data points)
%sigma:sigma for the hampel filter
%efs:sampling frequency
%sr:native sampling rate of bdf
%type=0 if gastric synthetic,1=GI bdf, 2=cardiac synthetic,3=cardiac exp

%Created on 23.11.2021
%Pacing or stimulation artefacts are identified by the hampel filter and
%are removed using a linear interpolation approach

%choosing electrode channels
if type==1
    dl=size(Data,1)-1;  %all recorded channels (minus trigger)
else
    dl=size(Data,1);
end

A1 = Data(1:dl,:); 

%time length
sigTest = Data(3,:);
t = 1/sr:1/sr:length(sigTest)/sr;

%remove baseline
if type==1
    for i = 1:size(A1,1)
        baseMed = medfilt1(A1(i,:),(sr*3)+1);
        base =sgolayfilt(baseMed,2,(sr*2)+1);
        A1(i,:)=A1(i,:)-base;
    end
elseif type==3
    for i = 1:size(A1,1)
        baseMed = medfilt1(A1(i,:),round((sr*1.5))+1);
        base =sgolayfilt(baseMed,2,round((sr*0.7))+1);
        A1(i,:)=A1(i,:)-base;
    end
else
end

%replacing the first and last few samples of the recording
if type==2||3 %cardiac
    startInt=round(efs*0.01); 
else %GI
    startInt=round(efs*0.1); 
end
lastSample=size(A1,2);
endInt=lastSample-startInt;
    
for s=1:startInt
    A1(:,s)= A1(:,(startInt+1)); %replacing with the sample immediately after 0.1s
    s=s+1;
end

for e=endInt:lastSample
    A1(:,e)=A1(:,(endInt-1)); %replacing with the sample 0.1s before the end of the recording
    e=e+1;
end

%taking the mean absolute of all signals 
sigAll =abs(A1);
sigAllMean = mean(sigAll,1);

% sigma=2*(range(sigAllMean)/std(sigAllMean));
sigma=40;

%Apply Hampel 
[yAll,iAll,xmedianAll,xsigmaAll] = hampel(sigAllMean,round(efs*k),sigma); %k=0.5s, sigma=10


%Detecting outliers
outly=sigAllMean(iAll);
outlt=t(iAll);

%find outliers that are more than 1s apart to find the start and end of stimuli 
endPos=[];
startPos=[];

%position of stimuli in time
if type==2
    endPt = find(diff(outlt)>0.5);
elseif type==3
    endPt = find(diff(outlt)>0.3);
else
    endPt = find(diff(outlt)>1);
end
endPos = [outlt(endPt) outlt(end)];
% endPos=endPos./efs;
startPos = [outlt(1) outlt(endPt+1)];
% startPos=startPos./efs;

%amplitude of the start and end stimuli points
endY = [outly(endPt) outly(end)];
startY = [outly(1) outly(endPt+1)];

%adding a buffer
% bufferst=0.01;
if type==2%ecg synthetic
    bufferst=0.005;
    bufferend=0.005;
elseif type==3
    bufferst=0.02;
    bufferend=0.1;   %was 0.02
else %GI synthetic
    bufferst=0.1; %was 0.02
    bufferend=0.2; %was 0.2 %longer buffer at the end because of the ringing
end
startPos=startPos-bufferst;
endPos=endPos+bufferend;
  
startPos=round(startPos*efs);
endPos=round(endPos*efs);

postStim_LI=A1;

%% Remove (Loop through all channels and sections)

for i = 1:size(A1,1)
% for i = 193:256
    i
    for j = 1:length(startPos)

        %Interpolate and replace 
        x = [startPos(j) endPos(j)];
        y = [postStim_LI(i,startPos(j)) postStim_LI(i,endPos(j))]; 
        
        xInt = [startPos(j):endPos(j)];
        yInt = interp1(x,y,xInt);
            
        postStim_LI(i,startPos(j):endPos(j)) = yInt;
    end    
end

startPos_fs=startPos;
endPos_fs=endPos;

artRem=postStim_LI;

end
