%% Importing the RAW data with offset
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
xlim = [0, 150]-18; ylim = [100, 2700]; zlim = [-.5 5];
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
                surface([range;range]-18,[y;y],[z;z],[z;z],...
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
                surface([range;range]-18,[y;y],[z;z],[z;z],...
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
%
% range = 0:6:iit*6;
% % importing 10 ng/ml TNF experiment
% load('20200906_CV_TNF10_DC3_All_J');
% % imaging started at 19:04:27
% % first feeding started at 19:49:11
% % second feeding started at 21:36:33
% osD = [45, 152];
% osI = [0, 0, 0, 0, .7, .7, .7, .7, 1.4, 1.4, 1.4, 1.4, 2.1, 2.1, 2.1, 2.1, 2.8, 2.8, 3.15, 3.15, 3.5, 3.5, 3.5, 3.5]; % in min
% tra10 = [];
% ind_pos10 = []; % This will store the index of the particular cell chosen and the cell position
% subplot(1,3,3);
% hold on
% for aa=[1,2,3,4,9,10,11,12] % [1,2,3,4,9,10,11,12,17,18,19,20,21,22]
%     if ~isempty(R{aa,1})
%         temp = R{aa,2}(:,2:3);
%         if contains(R{aa,1}, '1-')
%             os = -osD(1) + osI(aa);
%         else
%             os = -osD(2) + osI(aa);
%         end
%         for bb=1:size(temp,1)
%             x = (temp{bb,1}(:,1)-1)*interval + os;
%             y = mean(temp{bb,2}(:,4));
%
%             if withoffset == 1
%                 adj =  xm - ((lolo*xm^2-y)./lolo).^(1/2) -18; %-randi(18);
%             elseif withoffset == 0
%                 adj = -18;
%             end
%
%             if x(1) <= range(1)+adj && x(end) >= range(end)+adj+7
%                 qwe = temp{bb,1}(:,2);
%                 qwe = smooth(qwe,5, 'lowess');
%                 [~,I] = min(abs(x - (range(1)+adj)));
%                 z = qwe(I:I+iit)';
% %                 z = interp1(x',qwe',range+adj,'makima');
%                 if backadj == 1
%                     z = msbackadj(range', z', 'WindowSize', ws, 'Stepsize', ss);
%                     z = z';
%                 elseif backadj == 2
%                     z = z - mean(z(1:3));
%                 end
%                 tra10 = [tra10; z];
%                 ind_pos10 = [ind_pos10; aa, bb, y];
%
%                 z(z > thr) = thr;
%                 y = y*ones(1,length(range));
%                 surface([range;range]-18,[y;y],[z;z],[z;z],...
%                     'facecol','no',...
%                     'edgecol','interp',...
%                     'linew',1.4);
%
%             end
%         end
%     end
% end
% % plot3(xx, lolo*xm^2-lolo*(xm-xx).^2, thr*ones(1,length(xx)), 'r', 'linewidth', 2);
% hold off
% set(gca, 'XLim', xlim, 'YLim', ylim, 'ZLim', zlim, 'View', AG, 'YDir', 'reverse');
% set(gca, 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);
% xlabel('Time (min)', 'FontSize', 14); % ylabel('Distance (mm)', 'FontSize', 14); % zlabel('NF-kB');
% title('10 ng/ml', 'FontSize', 14);

% tra100 = tra100(:,1:31);
% temp = randperm(length(tra100),1000);
% tra100 = tra100(temp,:);
% ind_pos100 = ind_pos100(temp,:);

%% Linkage and Clustering for 100 ng/ml

m_c=10;

figure(1);clf;
% Z = linkage(tra100,'ward','euclidean');
Z = linkage(tra100,'weighted','correlation');
T = cluster(Z,'maxclust',m_c);
cutoff = median([Z(end-9,3) Z(end-8,3)]);
[dd, Ore1, Ore2] = dendrogram(Z,100,'ColorThreshold',cutoff);
set(dd, 'linewidth', 2);
% ylabel('Correlation Distance');
set(gca, 'fontsize', 12, 'XTick', []);

YY=[]; YY_ind=[];
[nr1,nc1]=size(tra100);
clusteringnum = m_c;
figure(2);clf;
counter=0;
for aa = [7, 9, 2, 1, 6, 5, 4, 3] % 1:clusteringnum
    counter=counter+1;
    subplot(2,4,counter)
    z = find(T == aa);
    temp = tra100(z,:);
    [nr,~]=size(temp);
    
    hold on
    if nr > 20
        toplot = randi(nr, 1,20);
    else
        toplot = 1:nr;
    end
    for bb=toplot
        plot(0:6:(nc1-1)*6, temp(bb,:), 'color', [.2+rand(1,3)*.8, .5], 'linewidth', .5);
    end
    plot(0:6:(nc1-1)*6, mean(temp), 'color', 'r', 'linewidth', 2);
    hold off
    
    set(gca, 'XLim', [0 (nc1-1)*6], 'Ylim', [0 3], 'fontsize', 11);
    %     title(['Group ' num2str(aa) ', ' num2str(round(nr./nr1*100,2)) '%']);
    title(['Group ' num2str(aa)]);
    %     if aa==1 || aa==3
    %         ylabel('NF-{\kappa}B Response');
    %     end
    %     if aa==3 || aa==4
    %         xlabel('Time (min)');
    %     end
    YY=[YY; ind_pos100(z,3)];
    YY_ind=[YY_ind; aa*ones(length(z),1)];
end

figure(3);clf;
% boxplot(YY, YY_ind);
ccc=[.7 .7 .7];
xx=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
LL=[7, 9, 2, 1, 6, 5, 4, 3];
hold on
for aa=1:8
    z=find(YY_ind(:,1)==LL(aa));
    Violin(YY(z), xx(aa), 'ViolinColor',ccc, 'BoxColor',[.8 0 0], 'ShowData', false, 'Width', .32, 'BoxWidth', 0.1, 'ShowMean', true); % This is actually showing median
end
hold off
set(gca,'xtick',[], 'yticklabel',[], 'TickLength', [.02, .02], 'LineWidth', .5, 'FontSize', 12, 'XLim', [0 9], 'YDir','reverse');


%% This is for example trace in Figure 1A

% boxplot(YY, YY_ind);
ccc=[.7 .7 .7];
xx=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
LL=[7, 9, 2, 1, 6, 5, 4, 3];

figure(4);clf;
aa=7;
z = find(T == aa);
temp = tra100(z,:);
[nr,~]=size(temp);
plot(0:6:(nc1-1)*6, temp(5,:), 'color', 'r', 'linewidth', 2);
set(gca, 'XLim', [0 (nc1-1)*6], 'Ylim', [-.5 3], 'fontsize', 11);



%% Linkage and Clustering for 30 ng/ml
figure(1);clf;
Z = linkage(tra30,'weighted','correlation');
T = cluster(Z,'maxclust',10);
cutoff = median([Z(end-9,3) Z(end-8,3)]);
[dd, T, outperm] = dendrogram(Z,'ColorThreshold',cutoff);
set(dd, 'linewidth', 1);
ylabel('Height');
set(gca, 'fontsize', 12);

YY=[]; YY_ind=[];
[nr1,nc1]=size(tra30);
clusteringnum = 10;
figure(2);clf;
for aa = 1:clusteringnum
    subplot(2,5,aa)
    z = find(T == aa);
    temp = tra30(z,:);
    [nr,~]=size(temp);
    
    hold on
    if nr > 20
        toplot = randi(nr, 1,20);
    else
        toplot = 1:nr;
    end
    for bb=toplot
        plot(0:6:(nc1-1)*6, temp(bb,:), 'color', [.2+rand(1,3)*.8, .5], 'linewidth', .5);
    end
    plot(0:6:(nc1-1)*6, mean(temp), 'color', 'r', 'linewidth', 2);
    hold off
    
    set(gca, 'XLim', [0 (nc1-1)*6], 'Ylim', [0 2], 'fontsize', 11);
    title(['Group ' num2str(aa) ', ' num2str(round(nr./nr1*100,2)) '%']);
    %     if aa==1 || aa==3
    %         ylabel('NF-{\kappa}B Response');
    %     end
    %     if aa==3 || aa==4
    %         xlabel('Time (min)');
    %     end
    YY=[YY; ind_pos30(z,3)];
    YY_ind=[YY_ind; aa*ones(length(z),1)];
end

figure(3);clf;
boxplot(YY, YY_ind);

%%
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

% thr = .1;
% z = find(pk100(:,1) < thr); pk100(z,:) = [];
% z = find(pk30(:,1) < thr); pk30(z,:) = [];
% z = find(pk10(:,1) < thr); pk10(z,:) = [];

z1=randsample(1:length(pk100),700);
z2=randsample(1:length(pk30),700);
z3=randsample(1:length(pk10),700);

figure(1);clf;
hold on
plot(pk100(z1,1), pk100(z1,2)./pk100(z1,1), '.', 'color', [.8 0 0 .5], 'MarkerSize', 7);
plot(pk30(z2,1), pk30(z2,2)./pk30(z2,1), '.', 'color', [0 .8 0 .5], 'MarkerSize', 7);
plot(pk10(z3,1), pk10(z3,2)./pk10(z3,1), '.', 'color', [0 0 .8 .5], 'MarkerSize', 7);
hold off
set(gca, 'YScale', 'log', 'XScale', 'log', 'XLim', [0.05 5], 'YLim', [0.05 5]);

figure(2);clf;
hold on
plot(pk100(z1,3), pk100(z1,4)./pk100(z1,3), '.', 'color', [.8 0 0 .5], 'MarkerSize', 7);
plot(pk30(z2,3), pk30(z2,4)./pk30(z2,3), '.', 'color', [0 .8 0 .5], 'MarkerSize', 7);
plot(pk10(z3,3), pk10(z3,4)./pk10(z3,3), '.', 'color', [0 0 .8 .5], 'MarkerSize', 7);
hold off
set(gca, 'YScale', 'log', 'XScale', 'log', 'XLim', [0.05 50], 'YLim', [0.05 5]);


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
    600,  800;
    1000,  1200;
    1400,  1600;
    1800,  2000;
    2200,  2400
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
plot(pk100m(:,3), pk100m(:,4)./pk100m(:,3), '-o', 'color', [.7 0 0 1], 'linewidth', 1.5);
plot(pk30m(:,3), pk30m(:,4)./pk30m(:,3), '-o', 'color', [0 .7 0 1], 'linewidth', 1.5);
plot(pk10m(:,3), pk10m(:,4)./pk10m(:,3), '-o', 'color', [0 0 .7 1], 'linewidth', 1.5);
hold off
figure(2);
hold on
plot(pk100m(:,7), pk100m(:,8)./pk100m(:,7), '-o', 'color', [.7 0 0 1], 'linewidth', 1.5);
plot(pk30m(:,7), pk30m(:,8)./pk30m(:,7), '-o', 'color', [0 .7 0 1], 'linewidth', 1.5);
plot(pk10m(:,7), pk10m(:,8)./pk10m(:,7), '-o', 'color', [0 0 .7 1], 'linewidth', 1.5);
hold off
