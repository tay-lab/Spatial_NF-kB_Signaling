% Function which compares results between different experiments.
function ChamberComparison()
clear all; close all; clc

%%%%% Options %%%%%
mov = 1;                % Show movie: on:1, off:0;
analysis = 1;           % Analysis.
handletime = 'dt';      % Time.

handS1    = '7_5';      % TNF concentration.
ens       = '10';       % Number of experiment in ensemble.
Threshold = 0.0436;     % Threshold in NFkB where the system is considered to be activated.
peakChange = 1.25;      % Difference between max and min for a peak to be registered.
fig_count = 1;          % Figure counter.
preChamber = 30;
%colormap('jet')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data to be analysed. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Data,TNF,NFkB,Gene1,Gene2,Gene3,TNFtotal]=loadData(handletime,handS1,ens);
Xcells = Data(1); Ycells = Data(2);

y1d = (-preChamber:(Ycells-1-preChamber))*sqrt(3/4)*15/1000;

% Calculate the positions of the hexagonal grid.
x = zeros(1,Xcells*Ycells);
y = zeros(1,Xcells*Ycells);
b = (0:(Xcells-1))*15/1000;

% Create x and y position for the cells.
for i = 0:(Ycells-1)
    if (mod(i,2) == 0)
        x(i*Xcells+1:(i+1)*Xcells) = b + .5*15/1000;
    else
        x(i*Xcells+1:(i+1)*Xcells) = b;
    end
    y(i*Xcells+1:(i+1)*Xcells) = (i-preChamber)*sqrt(3/4)*15/1000;
end

%%%%%%%%%%%%%%%%%
% Movie section %
%%%%%%%%%%%%%%%%%

if mov
    movie(x,y,Data,TNF,NFkB,handletime,handS1);
    close(figure(1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single experiment analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if analysis
    
    [avgNF,stdNF] = AvgMaxNFkB(NFkB,Xcells,0,Threshold);
    % Includes a plot of the average value only for activated cells,
    figure(fig_count); fig_count=fig_count+1;   % understood as a maximum level above a certain background.
    errorbar(y1d,avgNF,stdNF,'o','linewidth',2)
    set(gca,'box','off','linewidth',2); xlabel('Distance [mm]'); ylabel('Average NF-\kappaB [N_{tot}]')
    axis([0,1.3,0,1])
    
    figure(fig_count); fig_count=fig_count+1;
    xcell = 1:10;
    ycell = preChamber+20*ones(1,10);
    plot_cell(TNF,NFkB,xcell,ycell,Xcells)
    
    [distance,peak_time,oscillations,act,act_avg,act_std,osc_avg,osc_std,osc_med,...
        osc_25,osc_75,per_avg,per_std,per2_avg,per2_std,peakNum,peakPos,peakTime,...
        area_avg,area_std,area_med,area_25,area_75]=FirstPeak(NFkB,Xcells,Ycells,Threshold,preChamber,peakChange);
    
%     figure(fig_count); fig_count=fig_count+1;
%     plot(peak_time,distance,'o');
    
    figure(fig_count); fig_count=fig_count+1;
    errorbar(y1d,act_avg,act_std,'o','linewidth',2)
    set(gca,'box','off','linewidth',2); axis([0,1.4,-.05,1.05])
    xlabel('Distance [mm]');    ylabel('Fraction of activation');

%     figure(fig_count); fig_count=fig_count+1;
%     errorbar(y1d,area_avg,area_std,'o','linewidth',2)
%     set(gca,'box','off','linewidth',2); axis([0,1.4,-.5,3])
%     xlabel('Distance [mm]');    ylabel('Oscillations');

    figure(fig_count); fig_count=fig_count+1;
    errorbar(y1d,area_med,(area_med-area_25),(area_75-area_med),'o','linewidth',2)
    set(gca,'box','off','linewidth',2); axis([0,1.4,-.5,4])
    xlabel('Distance [mm]');    ylabel('Oscillations');

%     figure(fig_count); fig_count=fig_count+1;
%     errorbar(y1d,osc_avg,osc_std,'o','linewidth',2)
%     set(gca,'box','off','linewidth',2); axis([0,1.4,-.5,4])
%     xlabel('Distance [mm]');    ylabel('Area [N_{tot}hr]');

    figure(fig_count); fig_count=fig_count+1;
    errorbar(y1d,osc_med,(osc_med-osc_25),(osc_75-osc_med),'o','linewidth',2)
    set(gca,'box','off','linewidth',2); axis([0,1.4,-.5,4])
    xlabel('Distance [mm]');    ylabel('Oscillations');


    figure(fig_count); fig_count=fig_count+1;
    errorbar(y1d,per_avg,per_std,'o','linewidth',2)
    hold on
    errorbar(y1d,per2_avg,per2_std,'ro','linewidth',2)
    set(gca,'box','off','linewidth',2); axis([0,1.4,-.25,2.5])
    xlabel('Distance [mm]');    ylabel('Period [hr]');
    legend(['Time 1. to 2. peak';'Time 2. to 3. peak'])
    
    hfig = figure(fig_count); fig_count=fig_count+1;
    set(hfig,'Position',[200,200,250,550])
    h_act = act == 1;   h_act = reshape(h_act,[1,Xcells*Ycells]);
    h_deact = act == 0; h_deact = reshape(h_deact,[1,Xcells*Ycells]);
    h_act = h_act&(y>=0); h_deact = h_deact&(y>=0);
    y1d = (-preChamber:(Ycells-1-preChamber))*sqrt(3/4)*15/1000;
    
    plot(x(h_act),y(h_act),'r.','linewidth',3)
    hold on
    plot(x(h_deact),y(h_deact),'b.','linewidth',3)
    plot(act_avg*0.2325,y1d,'r--','linewidth',2)
    plot((1-act_avg)*0.2325,y1d,'b--','linewidth',2)
    axis([-0.05,0.2325+0.05,0,1.286])
    set(gca,'box','off','linewidth',2)
    title('NF-\kappaB activation')
    
    hfig=figure(fig_count); fig_count=fig_count+1;
    set(hfig,'Position',[200,200,250,550])
    avgNFkB=zeros(size(NFkB,1),(size(NFkB,2)-1)/Xcells);
    for i = 1:Xcells
        centerNFkB = i:Xcells:(((Ycells-1)*Xcells)+i);
        avgNFkB = avgNFkB + NFkB(:,1+centerNFkB);
    end
    avgNFkB=avgNFkB/Xcells; time=NFkB(200:300,1)*60;
    imagesc(y1d,time,avgNFkB(200:300,:))
    set(gca,'FontSize',15,'LineWidth',3)
    xlabel('Position[mm]','Fontsize',15);
    ylabel('Time[min]','Fontsize',15);
    title('NF-\kappaB [N_{tot}]')
    
    hfig=figure(fig_count); fig_count=fig_count+1;
    set(hfig,'Position',[200,200,250,550])
    COLORMAP(:,2)=[1 0 0];COLORMAP(:,3)=[0.75 0 0.25];COLORMAP(:,4)=[0.5 0 0.5];COLORMAP(:,5)=[0.25 0 0.75];COLORMAP(:,1)=[0 0 1];
    
    for j=1:max(peakNum)
        h1 = find(peakNum>=j);
        plot(y(h1),cellfun(@(v) v(j), peakTime(1,h1))*60,'.','LineWidth',2,'COLOR',COLORMAP(:,mod(j,5)+1))
        hold on
    end
    axis([0,1.3,0,300])
    set(gca,'FontSize',15,'LineWidth',3)
    xlabel('Position[mm]','Fontsize',15)
    ylabel('Time[min]','Fontsize',15)
    title('Peak time and position','Fontsize',12)
    set(gca,'YDir','reverse');  
    
end
end

% Function which loads data:
function [Data,TNF,NFkB,Gene1,Gene2,Gene3,TNFtotal]=loadData(handletime,handleS1,ens)

str1 = ['DATA_time_',num2str(handletime),'_S_',handleS1,'_',ens,'.txt'];
str2 = ['TNF_time_',num2str(handletime),'_S_',handleS1,'_',ens,'.txt'];
str3 = ['NFkB_time_',num2str(handletime),'_S_',handleS1,'_',ens,'.txt'];
str4 = ['Gene1_time_',num2str(handletime),'_S_',handleS1,'_',ens,'.txt'];
str5 = ['Gene2_time_',num2str(handletime),'_S_',handleS1,'_',ens,'.txt'];
str6 = ['Gene3_time_',num2str(handletime),'_S_',handleS1,'_',ens,'.txt'];
str7 = ['TotTNF_time_',num2str(handletime),'_S_',handleS1,'_',ens,'.txt'];

% Specify path if necessary. Remove old paths to avoid loading wrong files.
path1 = 'DATA\\';
addpath(path1);

Data = load(str1);
TNF = load(str2);
NFkB = load(str3);
Gene1 = load(str4);
Gene2 = load(str5);
Gene3 = load(str6);
TNFtotal = load(str7);
end

function [avgNF,stdNF] = AvgMaxNFkB(NFkB,Xcells,OnlyAct,threshold)
% NFkB data: Each row t, x1y1,x2y1,x3y1, ... , x10y1, x1y2, ... ,

NFkB = NFkB(:,2:end);   % Remove time entries.
maxNFkB = max(NFkB(201:end,:));    % Find the maximum NFkB activation in time after activation. 201: t= 0.
maxNFkB = reshape(maxNFkB,[Xcells,size(NFkB,2)/Xcells]);
if OnlyAct == 0
    avgNF = mean(maxNFkB);
    stdNF = sqrt(mean((maxNFkB-ones(Xcells,1)*avgNF).^2));
else
    boolAct = maxNFkB > threshold;
    number = sum(boolAct); number(number==0) = 1;
    avgNF = sum(maxNFkB.*boolAct)./number;
    avgNFmat = ones(Xcells,1)*avgNF;
    stdNF = sqrt(sum((maxNFkB.*boolAct-avgNFmat.*boolAct).^2)./number);
end

end

% Makes a movie of the chamber activation.
function movie(x,y,Data,TNF,NFkB,handletime,handleS1)

movstr = ['Chamber',num2str(handletime),'hrS1_',handleS1];
mov = VideoWriter(movstr);
mov.FrameRate = 7.5;
scrsz = get(0,'ScreenSize');
open(mov);
Xcells = Data(1); Ycells = Data(2);

hFig = figure(1);
set(hFig, 'Position', [200 100 500 630])

for j = 200 : length(TNF(:,1))
    
    subaxis('MarginTop',.15,'MarginLeft',.15,'MarginRight',.15)
    suptitle(['Time = ',sprintf('%0.01f',TNF(j,1)),'hr'])
    %%%%%%%%%%%%%
    %%%  TNF  %%%
    %%%%%%%%%%%%%
    subaxis(1,2,1)
    scatter(x,-y, 50, TNF(j,2:end)', 'filled', 'o')
    set(gca, 'box', 'off')
    set(gca,'XTick',[0,0.2])
    hcol=colorbar; caxis([0 10]);
    set(hcol, 'YTick', [0,2,4,6,8,10])
    set(hcol, 'YTickLabel', [0,2,4,6,8,10])
    set(hcol,'FontSize',15,'LineWidth',2)
    set(gca,'FontSize',15,'LineWidth',3)
    
    title(['TNF'])
    axis([0,x(end),-y(end),-y(1)+.05])
    axis equal
    
    hold off
    
    
    %%%%%%%%%%%%%
    %%% NF-kB %%%
    %%%%%%%%%%%%%
    subaxis(1,2,2)
    scatter(x,-y, 50, NFkB(j,2:end)', 'filled', 'o')
    set(gca,'XTick',[0,0.2])
    set(gca,'YTick',[2,5])
    set(gca, 'box', 'off')
    hcol=colorbar();
    caxis([0 .4]);
    set(hcol, 'YTick', [0,0.1,0.2,0.3,0.4,0.5])
    set(hcol, 'YTickLabel', [0,0.1,0.2,0.3,0.4,0.5])
    set(hcol,'FontSize',15,'LineWidth',2)
    set(gca,'FontSize',15,'LineWidth',3)
    
    title(['Nuclear NF-\kappaB'])
    axis([0,x(end),-y(end),-y(1)+.05])
    axis equal
    frame = getframe(gcf);
    writeVideo(mov,frame);
    %     pause();
end
close(mov);
end

% Plot traces of NFkB and TNF in individual cells as function of time. The
% inputs xcell and ycell are arrays of cells which should be compared.
function plot_cell(TNF,NFkB,xcell,ycell,Xcells) % xcell [1,Xcells], ycell [1,Ycells];
time = NFkB(:,1);
NFkB = NFkB(:,2:end); TNF = TNF(:,2:end);% Remove time entries.
TNFtrace = TNF(:,xcell + (ycell-1)*Xcells);
NFkBtrace = NFkB(:,xcell + (ycell-1)*Xcells);

subplot(2,1,1)
plot(time,TNFtrace,'LineWidth',2)
axis([-2,10,0,5])
set(gca,'box','off','linewidth',2)
title(['TNF-profile in distance ',num2str((ycell(1)-1)*sqrt(3/4)*15/1000,2),' mm from valve'])
ylabel('TNF [C_{TNF}]')

subplot(2,1,2)
plot(time,NFkBtrace,'LineWidth',2)
axis([-2,10,0,1])
set(gca,'box','off','linewidth',2)
title('NF\kappaB')
xlabel('Time [hr]');    ylabel('NF-\kappaB [N_{tot}]')
end

% Calculate timing of first peak
function [distance,peak_time,oscillations,act,act_avg,act_std,osc_avg,osc_std,osc_med,osc_25,osc_75,per_avg,per_std,per2_avg,per2_std,peakNum,peakPos,peakTime,area_avg,area_std,area_med,area_25,area_75]=FirstPeak(NFkB,Xcells,Ycells,Threshold,preChamber,peakChange)
time = NFkB(:,1);
NFkB = NFkB(:,2:end);

peak_time = zeros(1,Xcells*Ycells); peakNum = zeros(1,Xcells*Ycells); peakPos = 1:(Xcells*Ycells);
oscillations = zeros(1,Xcells*Ycells);
period = zeros(1,Xcells*Ycells);
period2 = zeros(1,Xcells*Ycells);
area = zeros(1,Xcells*Ycells);

distance = ((-preChamber:(Ycells-1-preChamber))*sqrt(3/4)*15/1000)'*ones(1,Xcells);
distance = reshape(distance',[1,Xcells*Ycells]);

for i = 1:Xcells*Ycells
    
    [val,pos] = findpeaks(NFkB(:,i),'MINPEAKDISTANCE',10,'MINPEAKHEIGHT',Threshold);
    [peaks_new_val,peaks_new_pos]=RealPeaks(pos,NFkB(:,i),time,0,peakChange);
    
    if isempty(peaks_new_val)==0
        peak_index = time(peaks_new_pos)>0;
        peakNum(i) = length(peak_index);
    else
        peak_index = [];
        peakNum(i) = 0;
    end
    
    if isempty(peak_index)
        peak_time(i)=time(end);
        peakTime{i}=0;
    else
        handle = find(peak_index);
        if isempty(handle)==0
            peak_time(i)=time(peaks_new_pos(handle(1)));
            peakTime{i}=time(peaks_new_pos(handle));
            oscillations(i)=length(handle);
            if length(handle)==2
                period(i)=time(peaks_new_pos(handle(2)))-time(peaks_new_pos(handle(1)));
                period2(i)=0;
            elseif length(handle)>=3
                period(i)=time(peaks_new_pos(handle(2)))-time(peaks_new_pos(handle(1)));
                period2(i)=time(peaks_new_pos(handle(3)))-time(peaks_new_pos(handle(2)));
            else
                period(i)=0;
                period2(i)=0;
            end
        else
            peak_time(i) = time(end);
            oscillations(i)=0;
            period(i)=0;
            period2(i)=0;
        end
    end
    
    
    % Find local maxima
    [peak_val,peak_pos] = findpeaks(NFkB(:,i),'MINPEAKDISTANCE',10,'MINPEAKHEIGHT',Threshold);
    [peaks_new_val,peaks_new_pos]=RealPeaks(peak_pos,NFkB(:,i),time,0,peakChange);
    
    if length(peaks_new_pos)==1
        [min_val,min_pos]=min(NFkB(peaks_new_pos:floor(end*3/4),i));
        thandle = and((time>=0),time<=time(peaks_new_pos+min_pos-1));
        area(i) = trapz(NFkB(thandle,i))-0.0436*(time(peaks_new_pos(1)+min_pos-1)-0);
    elseif length(peaks_new_pos)>1;
        [min_val,min_pos]=min(NFkB(peaks_new_pos(1):floor(end*3/4),i));
        thandle = and((time>=0),time<=time(peaks_new_pos(1)+min_pos-1));
        area(i) = trapz(NFkB(thandle,i))-0.0436*(time(peaks_new_pos(1)+min_pos-1)-0);
    elseif isempty(peaks_new_pos)==1;
        area(i) = 0;
    end
    
end

% Calculate activation probability and standard deviation:
act = oscillations > 0;
act = reshape(act,[Xcells,Ycells]);
act_avg = mean(act);
act_std = sqrt(act_avg.*(1-act_avg)/Xcells); % Standard deviation of the probability from binomial distributed data.

% Calculate the positions of the hexagonal grid.
x = zeros(1,Xcells*Ycells);
y = zeros(1,Xcells*Ycells);
b = (0:(Xcells-1))*15/1000;

% Create x and y position for the cells.
for i = 0:(Ycells-1)
    if (mod(i,2) == 0)
        x(i*Xcells+1:(i+1)*Xcells) = b + .5*15/1000;
    else
        x(i*Xcells+1:(i+1)*Xcells) = b;
    end
    y(i*Xcells+1:(i+1)*Xcells) = (i-preChamber)*sqrt(3/4)*15/1000;
end

% Calculate area of the first peak and standard deviation:
area = reshape(area,[Xcells,Ycells]);
area_avg = mean(area);
area_std = sqrt(mean((area-ones(Xcells,1)*area_avg).^2)); % Standard deviation of the probability from binomial distributed data.

area_med = median(area);
area_25  = quantile(area,.25);
area_75  = quantile(area,.75);

% Calculate mean and deviation of number of oscillations:
oscillations = reshape(oscillations,[Xcells,Ycells]);
osc_avg = mean(oscillations);
osc_std = sqrt(mean((oscillations-ones(Xcells,1)*osc_avg).^2));

osc_med = median(oscillations);
osc_25  = quantile(oscillations,.25);
osc_75  = quantile(oscillations,.75);

% Calculate mean and deviation of period between first two oscillations:
period = reshape(period,[Xcells,Ycells]);
number = sum(period>0);
number(number==0)=1;
per_avg = sum(period)./number;
boolPer = period > 0; per_avgmat = ones(Xcells,1)*per_avg;
per_std = sum((period.*boolPer-per_avgmat.*boolPer).^2)./number;

period2 = reshape(period2,[Xcells,Ycells]);
number = sum(period2>0);
number(number==0)=1;
per2_avg = sum(period2)./number;
boolPer = period2 > 0; per_avgmat = ones(Xcells,1)*per2_avg;
per2_std = sum((period2.*boolPer-per_avgmat.*boolPer).^2)./number;

end

