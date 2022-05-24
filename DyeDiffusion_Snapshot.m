clc
clear

imgadd = '';
[fname, pname] = uigetfile([imgadd '*.tif'], 'Pick a mat file to inspect.');
info = imfinfo([pname, fname]);
num_img = length(info);
I = {};
parfor frame_num = 1:num_img
    I{frame_num,1} = imread([pname, fname], frame_num);
end

%%

cmin = 3500; cmax = 9000;
Lpos = 395; Rpos = 630;
Tpos = 20;  Bpos = 1000;
% bg = [];
% for frame_num = 1:3
%     for aa=pos1:pos2
%         bg(aa,frame_num) = median(I{frame_num,1}(aa,Cen-thr:Cen+thr));
%     end
% end

% X = .7*2*(0:1:(num_img-4)); % Time
% Y = -10:5:5*(pos2-pos1+1-3); % position (px in 20X)

aa=10+40; % 10 is 0 min. Each frame is 1 min
J = imrotate(flipud(I{aa,1}),-0,'bilinear','crop');
J = J(Tpos:Bpos,Lpos:Rpos);
[nr,nc] = size(J);
X = 0.75*5*(0:1:nc-1);
Y = 0.75*5*(0:1:nr-1)-500;

figure(1); clf;
imagesc(X,Y,J, [cmin cmax]); axis image;
ylabel('Distance (\mum)', 'FontSize', 14);
set(gca,'xtick',[], 'YDir','reverse', 'TickLength', [.02, .02], 'LineWidth', 1, 'FontSize', 12);