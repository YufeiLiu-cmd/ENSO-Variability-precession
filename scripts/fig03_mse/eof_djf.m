%%%%%% 平均春和冬季
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end
clearvars -except data_root; clc; close all;
num =1 ;
for yyyy = 0:15:345
name = fullfile(data_root, sprintf('MSE_ave_month_%03d.mat', yyyy));
load(name, 'lat', 'lon', 'MSE_ave_monthly');
for month =1:12
    monthly(:,:,month)=mean(MSE_ave_monthly(:,:,month:12:end),3,'omitnan');
end
MSE_monthly_exp24(:,:,:,num)=monthly;
num=num+1
end
% [a,b,c,d] =size(MSE_monthly_exp24);
% ave = mean (MSE_monthly_exp24,4,'omitnan');
% summer = mean(ave(:,:,6:8),3,'omitnan');
% winter = mean(ave(:,:,[1,2,12]),3,'omitnan');
% save('season.mat','summer','winter','lat','lon');
%% 
winter_ipwp =squeeze(mean (MSE_monthly_exp24(:,:,[1,2,12],:),3,'omitnan'));
lonneed(1:10,:)=lon(183:192,:)-360;
lonneed(11:192,:)=lon(1:182,:);
lon=lonneed;
% 调整陆地在一起 
winter_ipwpeed(1:10,:,:,:)=winter_ipwp(183:192,:,:,:);
winter_ipwpeed(11:192,:,:,:)=winter_ipwp(1:182,:,:,:);
winter_ipwp= winter_ipwpeed;
dylat=41;
dylon=300;
idx_lon=find(lon >= -20 & lon <=dylon);
idx_lat=find(lat >= -dylat & lat <=dylat);
winter_ipwp = winter_ipwp(idx_lon,idx_lat,:);
lon_need = lon(idx_lon);
lat_need = lat(idx_lat);
[lon2d,lat2d]=meshgrid(lon_need,lat_need);
winter_ipwp_90 = winter_ipwp(:,:,7);
winter_ipwp_270 = winter_ipwp(:,:,19);
% %% 
% for i = 1:24
%     % 打开更大的绘图窗口
%     figure('Position', [100, 100, 1800, 1200]);
% 
%     % 获取数据
%     lock = winter_ipwp(:, :, i);
% 
%     % 绘图设置
%     m_proj('equidistant cylindrical', 'lon', [-20 dylon], 'lat', [-dylat+1 dylat-1]);
%     [c, h1] = m_contourf(lon2d, lat2d, lock', 250:1:350);
%     set(h1, 'LineColor', 'none');
%     hold on;
% 
%     m_coast('color', [0.7 0.7 0.7], 'linewidth', 2.5);
%     m_grid('box', 'on', 'linest', ':', 'linecolor', [0.7 0.7 0.7], ...
%         'linewidth', 0.1, 'tickdir', 'in', 'backcolor', 'w', ...
%         'ytick', -dylat+1:10:dylat-1, 'xtick', -20:40:dylon, ...
%         'YaxisLocation', 'right', 'Fontsize', 16);
% 
%     % 添加等高线
%     h3 = m_contour(lon2d, lat2d, lock', 250:5:350, 'LineWidth', 1.2, 'LineStyle', '-', 'Color', 'k');
%     clabel(h3, 'FontSize', 8, 'Color', 'k');
% 
%     % 加载配色方案
%     colorcode = load(fullfile(data_root, 'colormaps', 'hotcolr_19lev (3).rgb'));
%     colormap(gca, colorcode ./ 256);
% 
%     % 设置 colorbar
%     c = colorbar('southoutside');
%     caxis([300 350]);
%     ticks = 300:10:350;
%     set(c, 'Ticks', ticks, 'FontSize', 18);
%     set(get(c, 'Label'), 'String', 'MSE ANO', 'FontSize', 18);
% 
%     % 动态标题 & 文件名
%     phase_deg = (i-1)*15;
%     name = sprintf('DJF_MSE_%03ddeg', phase_deg);
%     title(sprintf('DJF - %03d° Precession Phase', phase_deg), 'FontSize', 20);
% 
%     % 设置背景
%     set(gcf, 'Color', [0.95 0.95 0.95]);
%     set(gca, 'BoxStyle', 'full', 'Box', 'on');
% 
%     % 保存图片
%     saveas(gcf, [name '.png']);
% 
%     % 或使用更高清版本：
%     % exportgraphics(gcf, [name '.png'], 'Resolution', 300);
% end

%% 做eof
var=winter_ipwp;
%% %%%%%%%%%%%%%%%%%%%%%%%%%EOF part: do not change %%%%%%%%%%%%%%%%%%%%%%
% please provide the var as input variable. The first two dimension is
% space, the third dimension is time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MON=1;YEAR=0;
EOF_hyang
%%%%%%%%%%%%%%%%%%%%EOF part: do not change %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this program will calculate the first three contriution terms 
% eof3d is 
%%%%%%%%%%%%%%%%%%%%EOF part: do not change %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eof_result = squeeze(eof3d(1,:,:));
    time_serise = pc(1,:);
%% 绘图
% 绘图部分
figure
width = 2000;
height = 600;
set(gcf, 'Position', [100, 100, width, height]);

%------------------ 子图1：时间序列 -------------------
subplot(1,10, [1 5]); % 占前两列
% 设置图像的宽度和高度
%  %设置图像的宽度和高度
subplotPos2 = get(gca, 'Position'); % 获取当前subplot的位置
subplotPos2(1) = 0.05; % 修改高度
subplotPos2(2) = 0.2; % 修改高度
subplotPos2(3) = 0.3; % 修改高度
subplotPos2(4) = 0.642; % 修改高度
set(gca, 'Position', subplotPos2);
x1 = 0:15:345;
xlim([0 345]);
xticks(0:15:345);
xticklabels({ ...
    '000°', '',  '030°','', ...
    '060°', '', '090°', '', ...
    '120°','', '150°', '', ...
    '180°', '',  '210°', '', ...
    '240°',  '', '270°', '', ...
    '300°', '', '330°',  ''
});
grid on;
% 左轴蓝色 
yyaxis left
y = time_serise;
%y(25) = y(1); % 闭合循环
plot(x1,y, 'b-', 'LineWidth', 2.5);
ylabel('Time series', 'FontSize', 16);
ylim([-0.4, 0.4]);
yticks(-0.4:0.1:0.4);


% 右轴红色
% yyaxis right
% late = 2;
% if month - late > 0
%     y2 = srad_exp24(:, month - late);
% else
%     y2 = srad_exp24(:, month + (12 - late));
% end
%y2 = srad_exp24(:,month+11);
%y2(25) = y2(1); % 闭合循环
% 超前2个月
% y2new(1:22)=y2(3:24);
% y2new(23:24)=y2(1:2);
% plot(x1, y2, 'r-', 'LineWidth', 2.5);
% ylabel('Solar radiation', 'FontSize', 16);
% ylim([300,390]);
% yticks(300:30:390);
% set(gca, 'FontSize',16, 'TickDir', 'in');
% title(sprintf('Month%02d.png', month));

%------------------ 子图2：EOF 空间图 -------------------
subplot(1, 10, [4 10]); % 占后四列

m_proj('equidistant cylindrical', 'lon', [-18 dylon], 'lat', [-dylat dylat]);
[c, h1] = m_contourf(lon2d, lat2d, eof_result', -30:0.05:30);
set(h1, 'LineColor', 'none');
hold on;
m_coast('color','[0.7 0.7 0.7]', 'linewidth', 2.5);
m_grid('box','on', 'linest', ':', 'linecolor',[0.7 0.7 0.7], ...
    'linewidth', 0.1, 'tickdir','in', 'backcolor','w', ...
    'ytick', -dylat:10:dylat, 'xtick', -20:40:dylon, ...
    'YaxisLocation','right', 'Fontsize', 16);
% h3 = m_contour(lon2d, lat2d,eof_result', 5:5:55, 'LineWidth', 1.2, 'LineStyle', '-', 'Color', 'k');
% h4 = m_contour(lon2d, lat2d,eof_result', [0], 'LineWidth', 2.5, 'LineStyle', '-', 'Color', 'k');
% h5 = m_contour(lon2d, lat2d,eof_result', -55:5:-5, 'LineWidth', 1.2, 'LineStyle', '--', 'Color', 'k');
% clabel(h3,  'FontSize', 8, 'Color', 'k');
% clabel(h4,  'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
% clabel(h5,  'FontSize', 8, 'Color', 'k');
% 配色和colorbar设置
colorcode = load(fullfile(data_root, 'colormaps', 'NCV_blu_red.rgb'));
colormap(gca, colorcode ./ 256);
c = colorbar('eastoutside');
caxis([-10 10]);
ticks = -10:10:10;
set(c, 'Ticks', ticks, 'FontSize', 18);
set(get(c, 'Label'), 'String', 'MSE(EOF)', 'FontSize', 18);
c.Position(1) = c.Position(1) + 0.03;  % 增加 X 方向偏移
% 背景统一
set(gcf, 'Color', [0.95 0.95 0.95]);
set(gca,'BoxStyle','full','Box','on')
set(gca,'DataAspectRatio',[10 5 1])
%  %设置图像的宽度和高度
subplotPos1 = get(gca, 'Position'); % 获取当前subplot的位置
%subplotPos(3) = 0.8; % 修改高度
subplotPos1(2) = 0.2; % 修改高度
set(gca, 'Position', subplotPos1);
subplotPos1 = get(gca, 'Position'); % 获取当前subplot的位置
%subplotPos(3) = 0.8; % 修改高度
subplotPos1(2) = 0.1100; % 修改高度
set(gca, 'Position', subplotPos1);
grid on;


% % 保存图像
print(gcf, sprintf('Precession_EOF1_MSE_WINTER.png', month), '-dpng');
title('winter:61.4%');
stop;

load(fullfile(data_root, 'aprc_sea.mat'),'winter');
winter_aprc =winter(idx_lon,idx_lat,:);
% h5 = m_contour(lon2d, lat2d, squeeze(winter_aprc(:,:,1)'), [8 8], ...
%     'LineWidth', 2, 'LineStyle', '-', 'Color', 'b');
% h6 = m_contour(lon2d, lat2d, squeeze(winter_aprc(:,:,2)'), [8 8], ...
%     'LineWidth', 2, 'LineStyle', '-', 'Color', 'r');
h6 = m_contour(lon2d, lat2d, winter_ipwp_90', [330.5 330.5], ...
    'LineWidth', 2, 'LineStyle', '-', 'Color', 'b');
h7 = m_contour(lon2d, lat2d, winter_ipwp_270',[330.5 330.5], ...
    'LineWidth', 2, 'LineStyle', '-', 'Color', 'r');
save('mse_eof_winter.mat','lon2d','lat2d','winter_ipwp_90','winter_ipwp_270','time_serise','eof_result','winter_aprc');
