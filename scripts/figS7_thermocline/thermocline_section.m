clc;clear;close all;
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end
folder = fullfile(data_root, 'data105and285thetao');
%%  循环遍历文件
num=1;
for yyyy = 105:180:285
filename = sprintf('thetao%03d.nc',yyyy);  % 生成文件名
dir= fullfile(folder, filename);
file_info = ncinfo(dir);  %获取ts文件信息
%  读数据
lon=ncread(dir,'lon');   %-180到180
lat=ncread(dir,'lat');
depth=ncread(dir,'depth');
thetao=ncread(dir,'thetao');
time=ncread(dir,'time');
% 调整为太平洋为中心
for i=1:size(lon,1);
if lon(i)<0;
lon(i)=lon(i)+360;
end
end
lon_1(1:180,1)=lon(181:360,1);
lon_1(181:360,1)=lon(1:180,1);
lon=lon_1;
clear lon_1;
t(1:180,:,:,:)=thetao(181:360,:,:,:);
t(181:360,:,:,:)=thetao(1:180,:,:,:);
thetao_need=t;
thetao_need(thetao_need <=-99) = NaN;
clear t;
% 变为20E-380
lon_need(1:341)=lon(20:360);
lon_need(342:360)=lon(1:19)+360;
lon=lon_need;
t(1:341,:,:,:)=thetao_need(20:360,:,:,:);
t(342:360,:,:,:)=thetao_need(1:19,:,:,:);
thetao=t;
clear thetao_need
clear t

% 选区域
idx_lat = find(lat >= -5 & lat <=5);
idx_lon = find(lon >= 40 & lon <= 280);
idx_depth=find(depth>=0 & depth<=230);
thetao_need=thetao(idx_lon,idx_lat,idx_depth,:);
[a,b,c,d]=size(thetao_need)
lon_need=lon(idx_lon);
depth_need=depth(idx_depth);
thetao_data=squeeze(mean(thetao_need,2,'omitnan'));
thermocline_exp24(:,:,:,num)=thetao_data;
num=num+1
end
%keyboard;  %% 储存两个极端相位 的 温跃层 
[a,b,c,d]=size(thermocline_exp24);
% 求均值
for i=1:2
    thetao_ave(:,:,i) = mean(squeeze(thermocline_exp24(:,:,:,i)),3,"omitnan");
end

% 求差异
    diff = squeeze(thetao_ave(:,:,1))-squeeze(thetao_ave(:,:,2));

% t检验
thetao105=reshape(squeeze(thermocline_exp24(:,:,:,1)),[240*16,1200]);
thetao285=reshape(squeeze(thermocline_exp24(:,:,:,2)),[240*16,1200]);

% 双样本t检验
[h,sig,ci]=ttest2(thetao105',thetao285',0.05);
h=reshape(h,[240,16]);
h(h==0)=nan;

save("thermocline.mat",'h','thetao_ave','lon_need','depth_need','diff');

