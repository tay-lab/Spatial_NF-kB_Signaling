clc
clear

warning('off','all')

% importing 100 ng/ml TNF experiment
load('20200903_CV_TNF100_DC3_All_J');
% imaging started at 23:25:21
% first feeding started at 00:03:50
% second feeding started at 01:53:53

backadj = 1; % 1 means yes, 0 means no
osD = [38, 148];
osI = [0, 0, 0, 0, .7, .7, .7, .7, 1.4, 1.4, 1.4, 1.4, 2.1, 2.1, 2.1, 2.1, 2.8, 2.8, 2.8, 2.8, 3.5, 3.5, 3.5, 3.5]; % in min
tra100 = [];
ind_pos100 = []; % This will store the index of the particular cell chosen and the cell position

thr = 5; AG = [0 90];
xlim = [0, 150]; ylim = [100, 2700]; zlim = [-.5 5];
ws = 120; ss = 120;

xm = 80; xx = 0:xm; lolo=.45;
iit = 25;
withoffset = 1; % 1 means yes, 0 means no

range = 0:6:iit*6; % in min
figure(10); clf;
subplot(1,3,1);
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
            y = mean(temp{bb,2}(:,4));
            if withoffset == 1 
                adj =  xm - ((lolo*xm^2-y)./lolo).^(1/2)-18; %-randi(18);
            elseif withoffset == 0
                adj = -18;
            end
            
            if x(1) <= range(1)+adj && x(end) >= range(end)+adj+7
                
                qwe = temp{bb,1}(:,2);
                qwe = smooth(qwe,5, 'lowess');
                [~,I] = min(abs(x - (range(1)+adj)));
                z = qwe(I:I+iit)';
%                 z = interp1(x',qwe',range+adj,'makima');
                if backadj == 1
                    z = msbackadj(range', z', 'WindowSize', ws, 'Stepsize', ss);
                    z = z';
                elseif backadj == 2
                    z = z - mean(z(1:3));
                end
                tra100 = [tra100; z];
                ind_pos100 = [ind_pos100; aa, bb, y];
                
                z(z > thr) = thr;
                y = y*ones(1,length(range));
                surface([range;range],[y;y],[z;z],[z;z],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',1);
            end
        end
    end
end
% plot3(xx, lolo*xm^2-lolo*(xm-xx).^2, thr*ones(1,length(xx)), 'r', 'linewidth', 2);
hold off
set(gca, 'XLim', xlim, 'YLim', ylim, 'ZLim', zlim, 'View', AG, 'YDir', 'reverse');
set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
xlabel('Time (min)', 'FontSize', 14); ylabel('Distance (\mum)', 'FontSize', 14); % zlabel('NF-kB');
title('100 ng/ml', 'FontSize', 14);

range = 0:6:iit*6;
% importing 30 ng/ml TNF experiment
load('20200909_CV_TNF30_DC3_All_J');
% imaging started at 23:01:42
% first feeding started at 23:40:41
% second feeding started at 01:29:57
osD = [39, 148];
osI = [0, 0, 0, 0, .7, .7, .7, .7, 1.4, 1.4, 1.4, 1.4, 3.5, 3.5, 3.5, 3.5, 2.8, 2.8, 2.8, 2.8, 2.1, 2.1, 2.1, 2.1]; % in min
tra30 = [];
ind_pos30 = []; % This will store the index of the particular cell chosen and the cell position
subplot(1,3,2);
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
            y = mean(temp{bb,2}(:,4));
            
            if withoffset == 1 
                adj =  xm - ((lolo*xm^2-y)./lolo).^(1/2) -18; %-randi(18);
            elseif withoffset == 0
                adj = -18;
            end
            
            if x(1) <= range(1)+adj && x(end) >= range(end)+adj+7
                
                qwe = temp{bb,1}(:,2);
                qwe = smooth(qwe,5, 'lowess');
                [~,I] = min(abs(x - (range(1)+adj)));
                z = qwe(I:I+iit)';
%                 z = interp1(x',qwe',range+adj,'makima');
                if backadj == 1
                    z = msbackadj(range', z', 'WindowSize', ws, 'Stepsize', ss);
                    z = z';
                elseif backadj == 2
                    z = z - mean(z(1:3));
                end
                tra30 = [tra30; z];
                ind_pos30 = [ind_pos30; aa, bb, y];
                
                z(z > thr) = thr;
                y = y*ones(1,length(range));
                surface([range;range],[y;y],[z;z],[z;z],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',1.2);
                
            end
        end
    end
end
% plot3(xx, lolo*xm^2-lolo*(xm-xx).^2, thr*ones(1,length(xx)), 'r', 'linewidth', 2);
hold off
set(gca, 'XLim', xlim, 'YLim', ylim, 'ZLim', zlim, 'View', AG, 'YDir', 'reverse');
set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
xlabel('Time (min)', 'FontSize', 14); % ylabel('Distance (\mum)', 'FontSize', 14); % zlabel('NF-kB');
title('30 ng/ml', 'FontSize', 14);

range = 0:6:iit*6;
% importing 10 ng/ml TNF experiment
load('20200906_CV_TNF10_DC3_All_J');
% imaging started at 19:04:27
% first feeding started at 19:49:11
% second feeding started at 21:36:33
osD = [45, 152];
osI = [0, 0, 0, 0, .7, .7, .7, .7, 1.4, 1.4, 1.4, 1.4, 2.1, 2.1, 2.1, 2.1, 2.8, 2.8, 3.15, 3.15, 3.5, 3.5, 3.5, 3.5]; % in min
tra10 = [];
ind_pos10 = []; % This will store the index of the particular cell chosen and the cell position
subplot(1,3,3);
hold on
for aa=[1,2,3,4,9,10,11,12] % [1,2,3,4,9,10,11,12,17,18,19,20,21,22]
    if ~isempty(R{aa,1})
        temp = R{aa,2}(:,2:3);
        if contains(R{aa,1}, '1-')
            os = -osD(1) + osI(aa);
        else
            os = -osD(2) + osI(aa);
        end
        for bb=1:size(temp,1)
            x = (temp{bb,1}(:,1)-1)*interval + os;
            y = mean(temp{bb,2}(:,4));
            
            if withoffset == 1 
                adj =  xm - ((lolo*xm^2-y)./lolo).^(1/2) -18; %-randi(18);
            elseif withoffset == 0
                adj = -18;
            end
            
            if x(1) <= range(1)+adj && x(end) >= range(end)+adj+7
                qwe = temp{bb,1}(:,2);
                qwe = smooth(qwe,5, 'lowess');
                [~,I] = min(abs(x - (range(1)+adj)));
                z = qwe(I:I+iit)';
%                 z = interp1(x',qwe',range+adj,'makima');
                if backadj == 1
                    z = msbackadj(range', z', 'WindowSize', ws, 'Stepsize', ss);
                    z = z';
                elseif backadj == 2
                    z = z - mean(z(1:3));
                end
                tra10 = [tra10; z];
                ind_pos10 = [ind_pos10; aa, bb, y];
                
                z(z > thr) = thr;
                y = y*ones(1,length(range));
                surface([range;range],[y;y],[z;z],[z;z],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',1.4);
                
            end
        end
    end
end
% plot3(xx, lolo*xm^2-lolo*(xm-xx).^2, thr*ones(1,length(xx)), 'r', 'linewidth', 2);
hold off
set(gca, 'XLim', xlim, 'YLim', ylim, 'ZLim', zlim, 'View', AG, 'YDir', 'reverse');
set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
xlabel('Time (min)', 'FontSize', 14); % ylabel('Distance (mm)', 'FontSize', 14); % zlabel('NF-kB');
title('10 ng/ml', 'FontSize', 14);

% tra100 = tra100(:,1:31);
% temp = randperm(length(tra100),1000);
% tra100 = tra100(temp,:);
% ind_pos100 = ind_pos100(temp,:);

%% Scatter Plot
% xrange1 = [2:1:13];
% xrange2 = [13:1:24];
xrange1 = [4:1:10]; % 8
alt_pos1 = 8;
xrange2 = [14:1:22]; % 17
alt_pos2 = 17;
AUCrange = 2;

pk100 = []; pk30 = []; pk10 = [];

for aa=1:size(tra100,1)
    [M,I] = max(tra100(aa,:));
    if ~isempty(find(xrange1 == I, 1)) || ~isempty(find(xrange2 == I, 1))
        [M1,I1] = max(tra100(aa,xrange1)); I1 = I1+xrange1(1)-1;
        [M2,I2] = max(tra100(aa,xrange2)); I2 = I2+xrange2(1)-1;
        pk100 = [pk100; M1, M2, sum(tra100(aa,I1+(-AUCrange:AUCrange))), sum(tra100(aa,I2+(-AUCrange:AUCrange))), ind_pos100(aa,3)];
    else
        pk100 = [pk100; tra100(aa,alt_pos1), tra100(aa,alt_pos2), sum(tra100(aa,alt_pos1+(-AUCrange:AUCrange))), sum(tra100(aa,alt_pos2+(-AUCrange:AUCrange))), ind_pos100(aa,3)];
    end
end
for aa=1:size(tra30,1)
    [M,I] = max(tra30(aa,:));
    if ~isempty(find(xrange1 == I, 1)) || ~isempty(find(xrange2 == I, 1))
        [M1,I1] = max(tra30(aa,xrange1)); I1 = I1+xrange1(1)-1;
        [M2,I2] = max(tra30(aa,xrange2)); I2 = I2+xrange2(1)-1;
        pk30 = [pk30; M1, M2, sum(tra30(aa,I1+(-AUCrange:AUCrange))), sum(tra30(aa,I2+(-AUCrange:AUCrange))), ind_pos30(aa,3)];
    else
        pk30 = [pk30; tra30(aa,alt_pos1), tra30(aa,alt_pos2), sum(tra30(aa,alt_pos1+(-AUCrange:AUCrange))), sum(tra30(aa,alt_pos2+(-AUCrange:AUCrange))), ind_pos30(aa,3)];
    end
end
for aa=1:size(tra10,1)
    [M,I] = max(tra10(aa,:));
    if ~isempty(find(xrange1 == I, 1)) || ~isempty(find(xrange2 == I, 1))
        [M1,I1] = max(tra10(aa,xrange1)); I1 = I1+xrange1(1)-1;
        [M2,I2] = max(tra10(aa,xrange2)); I2 = I2+xrange2(1)-1;
        pk10 = [pk10; M1, M2, sum(tra10(aa,I1+(-AUCrange:AUCrange))), sum(tra10(aa,I2+(-AUCrange:AUCrange))), ind_pos10(aa,3)];
    else
        pk10 = [pk10; tra10(aa,alt_pos1), tra10(aa,alt_pos2), sum(tra10(aa,alt_pos1+(-AUCrange:AUCrange))), sum(tra10(aa,alt_pos2+(-AUCrange:AUCrange))), ind_pos10(aa,3)];
    end
end

thr = .2;
z = find(pk100(:,4)./pk100(:,3) < thr); pk100(z,:) = [];
z = find(pk30(:,4)./pk30(:,3) < thr); pk30(z,:) = [];
z = find(pk10(:,4)./pk10(:,3) < thr); pk10(z,:) = [];

z1=randsample(1:length(pk100),1200);
z2=randsample(1:length(pk30),700);
z3=randsample(1:length(pk10),500);

sz=5; MFA=.5;
figure(1);clf;
% hold on
% s1 = scatter(pk100(z1,1), pk100(z1,2)./pk100(z1,1), sz, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.8, 0, 0], 'MarkerFaceAlpha', MFA);
% s2 = scatter(pk30(z2,1), pk30(z2,2)./pk30(z2,1), sz, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.8, 0, .8], 'MarkerFaceAlpha', MFA);
% s3 = scatter(pk10(z3,1), pk10(z3,2)./pk10(z3,1), sz, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0, 0, .8], 'MarkerFaceAlpha', MFA);
% hold off


figure(2);clf;
% hold on
% s4 = scatter(pk100(z1,3), pk100(z1,4)./pk100(z1,3), sz, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.8, 0, 0], 'MarkerFaceAlpha', MFA);
% s5 = scatter(pk30(z2,3), pk30(z2,4)./pk30(z2,3), sz, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.8, 0, .8], 'MarkerFaceAlpha', MFA);
% s6 = scatter(pk10(z3,3), pk10(z3,4)./pk10(z3,3), sz, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0, 0, .8], 'MarkerFaceAlpha', MFA);
% hold off

figure(3);clf;
hold on
data = log10([pk100(z1,3), pk100(z1,4)./pk100(z1,3)]);
s1 = scatter3(data(:,1),data(:,2), 0*ones(length(data(:,1)),1), sz, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.8, 0, 0], 'MarkerFaceAlpha', MFA);
C = cov(data);
M = mean(data);
z = [.2]'; % Vector of probabilities
[f,~] = Isolines(M,C,z);
f1 = fimplicit(f, 'linewidth', .5, 'color', [.8, 0, 0]);
hold off
grid()
set(gca,'XLim',[-1 1.3],'YLim', [-1 1],'ZLim', [0 1], 'ZTick', [], 'View', [-10 60]);
set(gca, 'TickLength', [.04, .04], 'LineWidth', 1.3, 'FontSize', 14);
% legend({'Data','Mean',label{:}});




pk100m = [];
pk30m = [];
pk10m = [];

tra = [];
% yrange = [
%     300,  500;
%     700,  900;
%     1100,  1300;
%     1500,  1700;
%     1900,  2100;
%     2300,  2500
%     ];
yrange = [
    200,  400;
    400,  700;
    700,  1100;
    1100,  1600;
    1600,  2200;
    2200,  2800;
%     2200,  2700;
%     2300,  2700;
%     2200,  2400
    ];

for aa=1:size(yrange,1)
    z = find(pk100(:,end) >= yrange(aa,1) & pk100(:,end) < yrange(aa,2));
    pk100m = [pk100m; median(pk100(z,1)), median(pk100(z,2)), mean(pk100(z,1)), mean(pk100(z,2)), median(pk100(z,3)), median(pk100(z,4)), mean(pk100(z,3)), mean(pk100(z,4))];
    z = find(pk30(:,end) >= yrange(aa,1) & pk30(:,end) < yrange(aa,2));
    pk30m = [pk30m; median(pk30(z,1)), median(pk30(z,2)), mean(pk30(z,1)), mean(pk30(z,2)), median(pk30(z,3)), median(pk30(z,4)), mean(pk30(z,3)), mean(pk30(z,4))];
    z = find(pk10(:,end) >= yrange(aa,1) & pk10(:,end) < yrange(aa,2));
    pk10m = [pk10m; median(pk10(z,1)), median(pk10(z,2)), mean(pk10(z,1)), mean(pk10(z,2)), median(pk10(z,3)), median(pk10(z,4)), mean(pk10(z,3)), mean(pk10(z,4))];
end

figure(1);
hold on
plot(pk100m(:,3), pk100m(:,4)./pk100m(:,3), '-^', 'color', [.8 0 0], 'linewidth', 1.5);
plot(pk30m(:,3), pk30m(:,4)./pk30m(:,3), '-^', 'color', [.8 0 .8], 'linewidth', 1.5);
plot(pk10m(:,3), pk10m(:,4)./pk10m(:,3), '-^', 'color', [0 0 .8], 'linewidth', 1.5);
hold off
set(gca, 'YScale', 'log', 'XScale', 'log', 'XLim', [0.08 4], 'YLim', [0.15 3]);
set(gca, 'TickLength', [.04, .04], 'LineWidth', 1.3, 'FontSize', 14);
xlabel('1st NF-\kappaB Peak Height', 'FontSize', 14); ylabel('Ratio (2nd/1st Peak Height)', 'FontSize', 14);




figure(2);
hold on
plot(pk100m(:,7), pk100m(:,8)./pk100m(:,7), '-^', 'color', [.8 0 0], 'linewidth', 1.5);
plot(pk30m(:,7), pk30m(:,8)./pk30m(:,7), '-^', 'color', [.8 0 .8], 'linewidth', 1.5);
plot(pk10m(:,7), pk10m(:,8)./pk10m(:,7), '-^', 'color', [0 0 .8], 'linewidth', 1.5);
hold off
set(gca, 'YScale', 'log', 'XScale', 'log', 'XLim', [0.3 15], 'YLim', [0.2 3]);
set(gca, 'TickLength', [.04, .04], 'LineWidth', 1.3, 'FontSize', 14);
xlabel('First NF-\kappaB Peak AUC', 'FontSize', 14); ylabel('Ratio (Second/First Peak AUC)', 'FontSize', 14);

%% IT Analysis (30 vs 100 ng/ml)
% getCC0 is the basic form, starting with equal probability for all samples.
% getCC1 is the one with global optimization.

yrange = [
    000,  400;
    400,  800;
    800,  1200;
    1200,  1600;
    1600,  2000;
    2000,  2400;
%    2400,  2800
    ];
% yrange = [
%     000,  300;
%     300,  600;
%     600,  900;
%     900,  1200;
%     1200,  1500;
%     1500,  1800;
%     1800,  2100;
%     2100,  2400;
%     % 2400,  2700
%     ];
dst = 5;
InfoScore = []; InfoScoreK = [];
thr = .2;

pk100a=pk100; pk30a=pk30; pk10a=pk10;
tra100a=tra100; tra30a=tra30; tra10a=tra10;
z = find(pk100a(:,1) < thr); pk100a(z,:) = []; tra100a(z,:) = [];
z = find(pk30a(:,1) < thr); pk30a(z,:) = []; tra30a(z,:) = [];
z = find(pk10a(:,1) < thr); pk10a(z,:) = []; tra10a(z,:) = [];
for aa=1:size(yrange,1)
    z = find(pk100a(:,end) >= yrange(aa,1) & pk100a(:,end) < yrange(aa,2));
    % temp1 = pk100(z,1:4);
    temp1 = [pk100a(z,2), pk100a(z,2)./pk100a(z,1)];
    temp11 = tra100a(z,:);
    for bb=1:size(yrange,1)
        z = find(pk30a(:,end) >= yrange(bb,1) & pk30a(:,end) < yrange(bb,2));
        % temp2 = pk30(z,1:4);
        temp2 = [pk30a(z,2), pk30a(z,2)./pk30a(z,1)];
        temp22 = tra30a(z,:);
        [I,Q,~]=getCC0({temp1',temp2'},dst,0);
        InfoScore(aa,bb) = I;
        
        tra = [temp11; temp22];
        Y = [zeros(size(temp11,1),1); ones(size(temp22,1),1)];
        t = templateTree('MaxNumSplits',20);
        Mdl = fitcensemble(tra,Y, 'Method','AdaBoostM1','KFold',10,'Learners',t); % , 'ScoreTransform','ismax'
        K = kfoldLoss(Mdl);
        InfoScoreK(aa,bb) = K;
    end
end


% varlabel = {'0~0.3','0.3~0.6','0.6~0.9','0.9~1.2','1.2~1.5','1.5~1.8','1.8~2.1','2.1~2.4'};
varlabel = {'0~0.4','0.4~0.8','0.8~1.2','1.2~1.6','1.6~2.0','2.0~2.4'};
z = find(InfoScore < .0001);
InfoScore(z) = 0;
figure(3);clf;
h = heatmap(varlabel,varlabel,round(InfoScore,2),'ColorLimits',[0 1]); colormap('autumn');
% h = heatmap(InfoScore); colormap('autumn');
set(gca, 'FontSize', 12);
xlabel('Distance (mm) in 30 ng/ml Sample'); ylabel('Distance (mm) in 100 ng/ml Sample');

figure(4);clf;
h = heatmap(varlabel,varlabel,round(1-InfoScoreK,2),'ColorLimits',[.5 1]); colormap('autumn');
% h = heatmap(1-InfoScoreK); colormap('autumn');
set(gca, 'FontSize', 12);
xlabel('Distance (mm) in 30 ng/ml Sample'); ylabel('Distance (mm) in 100 ng/ml Sample');

%% IT Analysis (10 vs 30 ng/ml)
% getCC0 is the basic form, starting with equal probability for all samples.
% getCC1 is the one with global optimization.

yrange = [
    000,  400;
    400,  800;
    800,  1200;
    1200,  1600;
    1600,  2000;
    2000,  2400;
%    2400,  2800
    ];
% yrange = [
%     000,  300;
%     300,  600;
%     600,  900;
%     900,  1200;
%     1200,  1500;
%     1500,  1800;
%     1800,  2100;
%     2100,  2400;
%     % 2400,  2700
%     ];
dst = 5;
InfoScore = []; InfoScoreK = [];
thr = .2;

pk100a=pk100; pk30a=pk30; pk10a=pk10;
tra100a=tra100; tra30a=tra30; tra10a=tra10;
z = find(pk100a(:,1) < thr); pk100a(z,:) = []; tra100a(z,:) = [];
z = find(pk30a(:,1) < thr); pk30a(z,:) = []; tra30a(z,:) = [];
z = find(pk10a(:,1) < thr); pk10a(z,:) = []; tra10a(z,:) = [];
for aa=1:size(yrange,1)
    z = find(pk30a(:,end) >= yrange(aa,1) & pk30a(:,end) < yrange(aa,2));
    % temp1 = pk100(z,1:4);
    temp1 = [pk30a(z,2), pk30a(z,2)./pk30a(z,1)];
    temp11 = tra30a(z,:);
    for bb=1:size(yrange,1)
        z = find(pk10a(:,end) >= yrange(bb,1) & pk10a(:,end) < yrange(bb,2));
        % temp2 = pk30(z,1:4);
        temp2 = [pk10a(z,2), pk10a(z,2)./pk10a(z,1)];
        temp22 = tra10a(z,:);
        [I,Q,~]=getCC0({temp1',temp2'},dst,0);
        InfoScore(aa,bb) = I;
        
        tra = [temp11; temp22];
        Y = [zeros(size(temp11,1),1); ones(size(temp22,1),1)];
        t = templateTree('MaxNumSplits',20);
        Mdl = fitcensemble(tra,Y, 'Method','AdaBoostM1','KFold',10,'Learners',t); % , 'ScoreTransform','ismax'
        K = kfoldLoss(Mdl);
        InfoScoreK(aa,bb) = K;
    end
end


% varlabel = {'0~0.3','0.3~0.6','0.6~0.9','0.9~1.2','1.2~1.5','1.5~1.8','1.8~2.1','2.1~2.4'};
varlabel = {'0~0.4','0.4~0.8','0.8~1.2','1.2~1.6','1.6~2.0','2.0~2.4'};
z = find(InfoScore < .0001);
InfoScore(z) = 0;
figure(5);clf;
h = heatmap(varlabel,varlabel,round(InfoScore,2),'ColorLimits',[0 1]); colormap('autumn');
% h = heatmap(InfoScore); colormap('autumn');
set(gca, 'FontSize', 12);
xlabel('Distance (mm) in 10 ng/ml Sample'); ylabel('Distance (mm) in 30 ng/ml Sample');

figure(6);clf;
h = heatmap(varlabel,varlabel,round(1-InfoScoreK,2),'ColorLimits',[.5 1]); colormap('autumn');
% h = heatmap(1-InfoScoreK); colormap('autumn');
set(gca, 'FontSize', 12);
xlabel('Distance (mm) in 10 ng/ml Sample'); ylabel('Distance (mm) in 30 ng/ml Sample');

%% From https://www.mathworks.com/matlabcentral/fileexchange/76676-isolines-quantiles-of-bivariate-normal-distribution

function[f,p] = Isolines(M,C,z)

% M: Mean vector, 1x2
% C: Covariance matrix, 2x2
% z: a nx2 matrix of points or a nx1 vector of probabilities

% Output:
% (a) case z is a vector of probabilities
% f: nx1 matrix, symbolic functions of the ellipses (isolines) corresponding to
% probabilities in z. I.e. draws 'outside' the ellipse area have
% probability <= z to occurr;
% p: p = z;

% (b) case z is a matrix of points
% f: nx1 matrix, simbolic function of the ellipse passing through each point in z
% (though each row)
% p: probability corresponding to each isoline

switch size(z,2)
    case 1
        type = 'Prob';
    case 2
        type = 'Point';
    otherwise
         error('Invalid input dimensions.')  
end

mx = M(1);
my = M(2);
S = inv(C);

a = S(1);
h = S(1,2);
b = S(2,2);

g = @(x,y) a*(x-mx).^2+b*(y-my).^2+2*h*(x-mx).*(y-my);

switch type
    case 'Prob'
        c = +2.*log(z);
        p = z;
    case 'Point'
        c = -g(z(:,1),z(:,2));
        p = exp(c/2);
end

syms x y
f = a*(x-mx).^2+b*(y-my).^2+2*h*(x-mx).*(y-my) + c == 0;

end