% Written by Minjun Son, PhD. in 2020 for Son et al 2021.
% This code 'HeatMap&Traces_Dose' is copyright 2020, Minjun Son.
% It should be considered 'freeware'- and may be distributed freely in its original form when properly attributed.

%% Heatmap for 100 ng/ml TNF source dose
clc
clear

warning('off','all')

% importing 100 ng/ml TNF experiment
load('FixedSource_TNF100_DC3');
% imaging started at 23:25:21
% first feeding started at 00:03:50
% second feeding started at 01:53:53

osD = [38, 148];
osI = [0, 0, 0, 0, .7, .7, .7, .7, 1.4, 1.4, 1.4, 1.4, 2.1, 2.1, 2.1, 2.1, 2.8, 2.8, 2.8, 2.8, 3.5, 3.5, 3.5, 3.5]; % in min
tra100 = [];
ind_pos100 = []; % This will store the index of the particular cell chosen and the cell position

thr = 5; AG = [0 90];
xlim = [-12, 270]; ylim = [.1, 2.7]; zlim = [-.5 5];


xm = 80; xx = 0:xm; lolo=.45;
withoffset = 1; % 1 means yes, 0 means no

iit = 45;
range = -12:6:iit*6; % in min
figure(1); clf;
% subplot(1,3,1);
hold on
for aa=1:12
    if ~isempty(R{aa,1})
        temp = R{aa,2}(:,2:3);
        if contains(R{aa,1}, '1-')
            os = -osD(1) + osI(aa);
        else
            os = -osD(2) + osI(aa);
        end
        for bb=1:size(temp,1)
            x = (temp{bb,1}(:,1)-1)*interval + os;
            y = mean(temp{bb,2}(:,4))./1000;
            if withoffset == 1 
                adj =  xm - ((lolo*xm^2-y)./lolo).^(1/2) -18; %-randi(18);
            elseif withoffset == 0
                adj = -18;
            end
            
            if x(1) <= range(1)+adj && x(end) >= range(end)+adj+7
                
                qwe = temp{bb,1}(:,2);
                qwe = smooth(qwe,5, 'lowess');
                [~,I] = min(abs(x - (range(1)+adj)));
                z = qwe(I:I+length(range)-1)';
                z = z - mean(z(1:3));
                
                tra100 = [tra100; z];
                ind_pos100 = [ind_pos100; aa, bb, y];
                
                z(z > thr) = thr;
                y = y*ones(1,length(range));
                surface([range+adj;range+adj],[y;y],[z;z],[z;z],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',1);
            end
        end
    end
end
hold off
set(gca, 'TickLength', [.01, .01], 'LineWidth', 1, 'FontSize', 12, 'TickDir', 'out');
set(gca, 'XLim', [min(range) max(range)]+adj, 'YLim', ylim, 'ZLim', zlim, 'View', AG, 'YDir', 'reverse');
xlabel('Time (min)', 'FontSize', 12); ylabel('Distance (mm)', 'FontSize', 12); % zlabel('NF-kB');
title('100 ng/ml', 'FontSize', 14);

%% Traces for 100 ng/ml TNF source dose
ylim1 = [-.2 3];
yrange = [
    200,  400;
    1000,  1200;
    2000,  2200;
    ];
yrange = yrange*.001;

load DFC_Dye_8h
[nr,nc] = size(I);
X = .7*2*((1:nc)-11); % Time
Y = 5*((1:nr)-1); % position (px in 20X)
Y = Y*.001;
I = I./1000;

aa = 1; % 1 for first section in the ragen, 3 for last section in the range
z = find(ind_pos100(:,3) >= yrange(aa,1) & ind_pos100(:,3) <= yrange(aa,2));
temp = tra100(z,:);
z = randsample(length(z),10);
tempm = mean(temp);
figure(2);clf;
hold on
plot(range+adj, temp(z,:), '-', 'color', [.3 .3 .3 .5], 'linewidth', 1);
plot(min(range):max(range), I(round((yrange(aa,1)+yrange(aa,2))*100), (1:length(min(range):max(range)))), '-', 'color', [0 0 1 1], 'linewidth', 1.5);
plot(range+adj, tempm, '-', 'color', [1 0 0 1], 'linewidth', 2);
hold off
set(gca, 'XLim', [min(range), max(range)]+adj, 'YLim', ylim1);
set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
xlabel('Time (min)', 'FontSize', 14); ylabel('NF-\kappaB', 'FontSize', 12);

yyaxis right 
ylabel('Fluor. Dye', 'FontSize', 14);
set(gca,'ycolor','b', 'ylim', [0 1]);

%% Heatmap for 30 ng/ml TNF source dose
load('FixedSource_TNF30_DC3');
% imaging started at 23:01:42
% first feeding started at 23:40:41
% second feeding started at 01:29:57
osD = [39, 148];
osI = [0, 0, 0, 0, .7, .7, .7, .7, 1.4, 1.4, 1.4, 1.4, 3.5, 3.5, 3.5, 3.5, 2.8, 2.8, 2.8, 2.8, 2.1, 2.1, 2.1, 2.1]; % in min
tra30 = [];
ind_pos30 = []; % This will store the index of the particular cell chosen and the cell position

osD = [39, 148];
osI = [0, 0, 0, 0, .7, .7, .7, .7, 1.4, 1.4, 1.4, 1.4, 3.5, 3.5, 3.5, 3.5, 2.8, 2.8, 2.8, 2.8, 2.1, 2.1, 2.1, 2.1]; % in min
tra30 = [];
ind_pos30 = []; % This will store the index of the particular cell chosen and the cell position

thr = 5; AG = [0 90];
xlim = [-12, 240]; ylim = [.1, 2.7]; zlim = [-.5 5];

xm = 80; xx = 0:xm; lolo=.45;
withoffset = 0; % 1 means yes, 0 means no

iit = 45;
range = -12:6:iit*6; % in min
figure(3); clf;
hold on
for aa=1:12
    if ~isempty(R{aa,1})
        temp = R{aa,2}(:,2:3);
        if contains(R{aa,1}, '1-')
            os = -osD(1) + osI(aa);
        else
            os = -osD(2) + osI(aa);
        end
        for bb=1:size(temp,1)
            x = (temp{bb,1}(:,1)-1)*interval + os;
            y = mean(temp{bb,2}(:,4))./1000;
            if withoffset == 1 
                adj =  xm - ((lolo*xm^2-y)./lolo).^(1/2) -18; %-randi(18);
            elseif withoffset == 0
                adj = -18;
            end
            
            if x(1) <= range(1)+adj && x(end) >= range(end)+adj+7
                
                qwe = temp{bb,1}(:,2);
                qwe = smooth(qwe,5, 'lowess');
                [~,I] = min(abs(x - (range(1)+adj)));
                z = qwe(I:I+length(range)-1)';
                z = z - mean(z(1:3));
                
                tra30 = [tra30; z];
                ind_pos30 = [ind_pos30; aa, bb, y];
                
                z(z > thr) = thr;
                y = y*ones(1,length(range));
                surface([range+adj;range+adj],[y;y],[z;z],[z;z],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',1.3);
                if y(1) > 1.401 && y(1) < 1.403
                    y = 1.42*ones(1,length(range));
                    surface([range+adj;range+adj],[y;y],[z;z],[z;z],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',1.3);
                end
            end
        end
    end
end

% plot3(xx, lolo*xm^2-lolo*(xm-xx).^2, thr*ones(1,length(xx)), 'r', 'linewidth', 2);
hold off
set(gca, 'TickLength', [.01, .01], 'LineWidth', 1, 'FontSize', 12, 'TickDir', 'out');
set(gca, 'XLim', [min(range) max(range)]+adj, 'YLim', ylim, 'ZLim', zlim, 'View', AG, 'YDir', 'reverse');
xlabel('Time (min)', 'FontSize', 14); ylabel('Distance (mm)', 'FontSize', 14); % zlabel('NF-kB');
title('30 ng/ml', 'FontSize', 14);

%% Traces for 30 ng/ml TNF source dose
ylim1 = [-.2 2];
yrange = [
    200,  400;
    1000,  1200;
    2000,  2200;
    ];
yrange = yrange*.001;

load DFC_Dye_8h
[nr,nc] = size(I);
X = .7*2*((1:nc)-11); % Time
Y = 5*((1:nr)-1); % position (px in 20X)
Y = Y*.001;
I = I./1000;

aa = 1; % 1 for first section in the ragen, 3 for last section in the range
z = find(ind_pos30(:,3) >= yrange(aa,1) & ind_pos30(:,3) <= yrange(aa,2));
temp = tra30(z,:);
z = randsample(length(z),10);
tempm = mean(temp);
figure(4);clf;
hold on
plot(range+adj, temp(z,:), '-', 'color', [.3 .3 .3 .5], 'linewidth', 1);
plot(min(range):max(range), .3*I(round((yrange(aa,1)+yrange(aa,2))*100), (1:length(min(range):max(range)))), '-', 'color', [0 0 1 1], 'linewidth', 1.5);
plot(range+adj, tempm, '-', 'color', [1 0 0 1], 'linewidth', 2);
hold off
set(gca, 'XLim', [min(range), max(range)]+adj, 'YLim', ylim1);
set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
xlabel('Time (min)', 'FontSize', 14); ylabel('NF-\kappaB', 'FontSize', 14);

yyaxis right 
ylabel('Fluor. Dye', 'FontSize', 14);
set(gca,'ycolor','b', 'ylim', [0 1]);

%% Heatmap for 10 ng/ml TNF source dose
% importing 10 ng/ml TNF experiment
load('FixedSource_TNF10_DC3');
% imaging started at 19:04:27
% first feeding started at 19:49:11
% second feeding started at 21:36:33
osD = [45, 152];
osI = [0, 0, 0, 0, .7, .7, .7, .7, 1.4, 1.4, 1.4, 1.4, 2.1, 2.1, 2.1, 2.1, 2.8, 2.8, 3.15, 3.15, 3.5, 3.5, 3.5, 3.5]; % in min
tra10 = [];
ind_pos10 = []; % This will store the index of the particular cell chosen and the cell position

osD = [39, 148];
osI = [0, 0, 0, 0, .7, .7, .7, .7, 1.4, 1.4, 1.4, 1.4, 3.5, 3.5, 3.5, 3.5, 2.8, 2.8, 2.8, 2.8, 2.1, 2.1, 2.1, 2.1]; % in min
tra30 = [];
ind_pos30 = []; % This will store the index of the particular cell chosen and the cell position

thr = 5; AG = [0 90];
xlim = [-12, 240]; ylim = [.1, 2.7]; zlim = [-.5 5];

xm = 80; xx = 0:xm; lolo=.45;
withoffset = 0; % 1 means yes, 0 means no

iit = 45;
range = -12:6:iit*6; % in min
figure(5); clf;
hold on
for aa=[1,2,3,4,9,10,11,12]
    if ~isempty(R{aa,1})
        temp = R{aa,2}(:,2:3);
        if contains(R{aa,1}, '1-')
            os = -osD(1) + osI(aa);
        else
            os = -osD(2) + osI(aa);
        end
        for bb=1:size(temp,1)
            x = (temp{bb,1}(:,1)-1)*interval + os;
            y = mean(temp{bb,2}(:,4))./1000;
            if withoffset == 1 
                adj =  xm - ((lolo*xm^2-y)./lolo).^(1/2) -18; %-randi(18);
            elseif withoffset == 0
                adj = -18;
            end
            
            if x(1) <= range(1)+adj && x(end) >= range(end)+adj+7
                
                qwe = temp{bb,1}(:,2);
                qwe = smooth(qwe,5, 'lowess');
                [~,I] = min(abs(x - (range(1)+adj)));
                z = qwe(I:I+length(range)-1)';
                z = z - mean(z(1:3));
                
                tra10 = [tra10; z];
                ind_pos10 = [ind_pos10; aa, bb, y];
                
                z(z > thr) = thr;
                y = y*ones(1,length(range));
                surface([range+adj;range+adj],[y;y],[z;z],[z;z],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',1.5);
            end
        end
    end
end

% plot3(xx, lolo*xm^2-lolo*(xm-xx).^2, thr*ones(1,length(xx)), 'r', 'linewidth', 2);
hold off
set(gca, 'TickLength', [.01, .01], 'LineWidth', 1, 'FontSize', 12, 'TickDir', 'out');
set(gca, 'XLim', [min(range) max(range)]+adj, 'YLim', ylim, 'ZLim', zlim, 'View', AG, 'YDir', 'reverse');
xlabel('Time (min)', 'FontSize', 14); ylabel('Distance (mm)', 'FontSize', 14);
title('10 ng/ml', 'FontSize', 14);

%% Traces for 10 ng/ml TNF source dose
ylim1 = [-.2 1.5];
yrange = [
    200,  400;
    1000,  1200;
    2000,  2200;
    ];
yrange = yrange*.001;

load DFC_Dye_8h
[nr,nc] = size(I);
X = .7*2*((1:nc)-11); % Time
Y = 5*((1:nr)-1); % position (px in 20X)
Y = Y*.001;
I = I./1000;

aa = 1; % 1 for first section in the ragen, 3 for last section in the range
z = find(ind_pos10(:,3) >= yrange(aa,1) & ind_pos10(:,3) <= yrange(aa,2));
temp = tra10(z,:);
z = randsample(length(z),10);
tempm = mean(temp);
figure(6);clf;
hold on
plot(range+adj, temp(z,:)+.0, '-', 'color', [.3 .3 .3 .5], 'linewidth', 1);
plot(min(range):max(range), .1*I(round((yrange(aa,1)+yrange(aa,2))*100), 0+(1:length(min(range):max(range))))+.03, '-', 'color', [0 0 1 1], 'linewidth', 1.5);
plot(range+adj, tempm+.0, '-', 'color', [1 0 0 1], 'linewidth', 2);
hold off
set(gca, 'XLim', [min(range), max(range)]+adj, 'YLim', ylim1);
set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
xlabel('Time (min)', 'FontSize', 14); ylabel('NF-\kappaB', 'FontSize', 14);

yyaxis right 
ylabel('Fluor. Dye', 'FontSize', 14);
set(gca,'ycolor','b', 'ylim', [0 1]);

