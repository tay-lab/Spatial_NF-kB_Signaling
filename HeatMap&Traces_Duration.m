% Written by Minjun Son, PhD. in 2020 for Son et al 2021.
% This code 'HeatMap&Traces_Duration' is copyright 2020, Minjun Son.
% It should be considered 'freeware'- and may be distributed freely in its original form when properly attributed.

%% Heatmap for 15min release duration with 100 ng/ml source TNF
clc
clear
warning('off','all')

load('DifferentDuration_TNF100_DC2'); % For 15min
AG = [0 90];
thr = 5;
thickness = 2.5;
adj = -42;
xlim = [-12, 240]; ylim = [.1, 2.7]; zlim = [-.5 5];
iit = 60;
range = -18:6:iit*6; % in min
tra100 = [];
ind_pos100 = []; % This will store the index of the particular cell chosen and the cell position

figure(1); clf;
hold on
for aa=[2,3,4]
    if ~isempty(R{aa,1})
        temp = R{aa,2}(:,2:3);
        for bb=1:size(temp,1)
            x = (temp{bb,1}(:,1)-1)*interval;
            y = round(mean(temp{bb,2}(:,4))./1000,2);
            
            if x(1)+adj <= range(1) && x(end)+adj >= range(end)
                
                qwe = temp{bb,1}(:,2);
                qwe = smooth(qwe,5, 'lowess');
                
                [~,I] = min(abs(x+adj));
                z = qwe((I-3):(I+length(range)-4))';
                z = z - mean(z(1:3));
                
                tra100 = [tra100; z];
                z(z > thr) = thr;
                ind_pos100 = [ind_pos100; aa, bb, y];
                y = y*ones(1,length(z));
                surface([range;range],[y;y],[z;z],[z;z],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',thickness);
            end
            
        end
    end
end

hold off
set(gca, 'TickLength', [.01, .01], 'LineWidth', 1, 'FontSize', 12, 'TickDir', 'out');
set(gca, 'XLim', [min(range) max(range)], 'YLim', ylim, 'ZLim', zlim, 'View', AG, 'YDir', 'reverse');
xlabel('Time (min)', 'FontSize', 14); ylabel('Distance (mm)', 'FontSize', 14); % zlabel('NF-kB');
title('15 min', 'FontSize', 14);

%% Traces for 15min release duration with 100 ng/ml source TNF
ylim1 = [-.2 3];

yrange = [
    100,  400;
    1000,  1200;
    1900,  2300;
    ];
yrange = yrange*.001;

load DFC_Dye_15min
[nr,nc] = size(I);
X = .7*2*((1:nc)-11); % Time
Y = 5*((1:nr)-1); % position (px in 20X)
Y = Y*.001;
I = I./1000;

aa = 1; % 1 for the first position, 3 for the last position
z = find(ind_pos100(:,3) >= yrange(aa,1) & ind_pos100(:,3) <= yrange(aa,2));
temp = tra100(z,:);
z = randsample(length(z),10);
tempm = mean(temp);
figure(2);clf;
hold on
plot(range, temp(z,:), '-', 'color', [.3 .3 .3 .5], 'linewidth', 1);
plot(1-12:361-12, I(round((yrange(aa,1)+yrange(aa,2))*100), (1:361)), '-', 'color', [0 0 1 1], 'linewidth', 1.5);
% plot(range, tempm+.1, '-', 'color', [1 0 0 1], 'linewidth', 2);
hold off
set(gca, 'XLim', [min(range), max(range)], 'YLim', ylim1);
set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
xlabel('Time (min)', 'FontSize', 14); ylabel('NF-\kappaB', 'FontSize', 14);

yyaxis right
ylabel('Fluor. Dye', 'FontSize', 14);
set(gca,'ycolor','b', 'ylim', [0 1]);

%% Heatmap for 30min release duration with 100 ng/ml source TNF
AG = [0 90];
thr = 5;
thickness = 2.5;
adj = -42;
xlim = [-12, 240]; ylim = [.1, 2.7]; zlim = [-.5 5];
iit = 60;
range = -18:6:iit*6; % in min
tra100 = [];
ind_pos100 = []; % This will store the index of the particular cell chosen and the cell position

figure(3); clf;
hold on
for aa=[5,6,7,8]
    if ~isempty(R{aa,1})
        temp = R{aa,2}(:,2:3);
        for bb=1:size(temp,1)
            x = (temp{bb,1}(:,1)-1)*interval;
            y = round(mean(temp{bb,2}(:,4))./1000,2);
            
            if x(1)+adj <= range(1) && x(end)+adj >= range(end)
                
                qwe = temp{bb,1}(:,2);
                qwe = smooth(qwe,5, 'lowess');
                [~,I] = min(abs(x+adj));
                z = qwe((I-3):(I+length(range)-4))';
                z = z - mean(z(1:3));
                
                tra100 = [tra100; z];
                z(z > thr) = thr;
                ind_pos100 = [ind_pos100; aa, bb, y];
                y = y*ones(1,length(z));
                surface([range;range],[y;y],[z;z],[z;z],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',thickness);
            end
            
        end
    end
end

hold off
set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
set(gca, 'XLim', [min(range) max(range)], 'YLim', ylim, 'ZLim', zlim, 'View', AG, 'YDir', 'reverse');
xlabel('Time (min)', 'FontSize', 14); ylabel('Distance (mm)', 'FontSize', 14); % zlabel('NF-kB');
title('30 min', 'FontSize', 14);


%% Traces for 30min release duration with 100 ng/ml source TNF
ylim1 = [-.2 3];
yrange = [
    200,  400;
    1000,  1200;
    1900,  2300;
    ];
yrange = yrange*.001;

load DFC_Dye_30min
[nr,nc] = size(I);
X = .7*2*((1:nc)-11); % Time
Y = 5*((1:nr)-1); % position (px in 20X)
Y = Y*.001;
I = I./1000;

aa = 1; % 1 for the first position, 3 for the last position
z = find(ind_pos100(:,3) >= yrange(aa,1) & ind_pos100(:,3) <= yrange(aa,2));
temp = tra100(z,:);
z = randsample(length(z),10);
tempm = mean(temp);
figure(4);clf;
hold on
plot(range, temp(z,:), '-', 'color', [.3 .3 .3 .5], 'linewidth', 1);
plot(1-12:361-12, I(round((yrange(aa,1)+yrange(aa,2))*100), (1:361)), '-', 'color', [0 0 1 1], 'linewidth', 1.5);
% plot(range, tempm+.1, '-', 'color', [1 0 0 1], 'linewidth', 2);
hold off
set(gca, 'XLim', [min(range), max(range)], 'YLim', ylim1);
set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
xlabel('Time (min)', 'FontSize', 14); ylabel('NF-\kappaB', 'FontSize', 14);

yyaxis right
ylabel('Fluor. Dye', 'FontSize', 14);
set(gca,'ycolor','b', 'ylim', [0 1]);


%% Heatmap for 60min release duration with 100 ng/ml source TNF
AG = [0 90];
thr = 5;
thickness = 2.5;
adj = -42;
xlim = [-12, 240]; ylim = [.1, 2.7]; zlim = [-.5 5];
iit = 60;
range = -18:6:iit*6; % in min
tra100 = [];
ind_pos100 = []; % This will store the index of the particular cell chosen and the cell position

figure(5); clf;
hold on
for aa=[9,10,11,12]
    if ~isempty(R{aa,1})
        temp = R{aa,2}(:,2:3);
        for bb=1:size(temp,1)
            x = (temp{bb,1}(:,1)-1)*interval;
            y = round(mean(temp{bb,2}(:,4))./1000,2);
            
            if x(1)+adj <= range(1) && x(end)+adj >= range(end)
                
                qwe = temp{bb,1}(:,2);
                qwe = smooth(qwe,5, 'lowess');
                [~,I] = min(abs(x+adj));
                z = qwe((I-3):(I+length(range)-4))';
                z = z - mean(z(1:3));
                tra100 = [tra100; z];
                z(z > thr) = thr;
                ind_pos100 = [ind_pos100; aa, bb, y];
                y = y*ones(1,length(z));
                surface([range;range],[y;y],[z;z],[z;z],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',thickness);
            end
            
        end
    end
end

hold off
set(gca, 'TickLength', [.01, .01], 'LineWidth', 1, 'FontSize', 12, 'TickDir', 'out');
set(gca, 'XLim', [min(range) max(range)], 'YLim', ylim, 'ZLim', zlim, 'View', AG, 'YDir', 'reverse');
xlabel('Time (min)', 'FontSize', 14); ylabel('Distance (mm)', 'FontSize', 14); % zlabel('NF-kB');
title('60 min', 'FontSize', 14);

%% Traces for 60min release duration with 100 ng/ml source TNF
ylim1 = [-.2 3];

yrange = [
    200,  400;
    1000,  1200;
    1900,  2300;
    ];
yrange = yrange*.001;

load DFC_Dye_1h
[nr,nc] = size(I);
X = .7*2*((1:nc)-11); % Time
Y = 5*((1:nr)-1); % position (px in 20X)
Y = Y*.001;
I = I./1000;

aa = 1; % 1 for the first position, 3 for the last position
z = find(ind_pos100(:,3) >= yrange(aa,1) & ind_pos100(:,3) <= yrange(aa,2));
temp = tra100(z,:);
z = randsample(length(z),10);
tempm = mean(temp);
figure(6);clf;
hold on
plot(range, temp(z,:), '-', 'color', [.3 .3 .3 .5], 'linewidth', 1);
plot(1-12:361-12, I(round((yrange(aa,1)+yrange(aa,2))*100), (1:361)), '-', 'color', [0 0 1 1], 'linewidth', 1.5);
% plot(range, tempm+.1, '-', 'color', [1 0 0 1], 'linewidth', 2);
hold off
set(gca, 'XLim', [min(range), max(range)], 'YLim', ylim1);
set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
xlabel('Time (min)', 'FontSize', 14); ylabel('NF-\kappaB', 'FontSize', 14);

yyaxis right
ylabel('Fluor. Dye', 'FontSize', 14);
set(gca,'ycolor','b', 'ylim', [0 1]);


