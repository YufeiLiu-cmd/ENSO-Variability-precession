
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end
clearvars -except data_root; clc; close all;

% % 文件夹路径
folder = fullfile(data_root, 'Precession24exp_0058_mean');
%  循环遍历文件
num=1
for yyyy = 0:15:345
   filename = sprintf('%03d_echam.nc', yyyy);  % 生成文件名
   dir= fullfile(folder, filename);
     file_info = ncinfo(dir);  % 获取文件信息
     lon = ncread(dir, 'lon');
     lat = ncread(dir, 'lat');
     tsw = ncread(dir, 'tsw'); % 读取变量 tsw 的数据
     aprc = ncread(dir, 'aprc'); % 
     aprc = aprc*86400;
     aprc_land = aprc;
     t=ncread(dir, 'time');

% 去掉陆地 
% file2=fullfile(data_root, 'orbzeroslm.nc');
slm=ncread(dir,'slm');
[m,n,t]=size(slm);
for i=1:m
    for j=1:n
        if slm(i,j,1)>0.5
            tsw(i,j,:)=nan;
        end
    end
end
clear i j ;
for q=1:m
    for w=1:n
        if slm(q,w,1)<0.5
            aprc_land(q,w,:)=nan;
        end
    end
end
clear q w;

tsw_exp24(:,:,:,num)=tsw;
aprc_land24(:,:,:,num)=aprc_land;
aprc_exp24(:,:,:,num)=aprc;
num=num+1;
end
ano_aprc = mean(squeeze(aprc_land24(:,:,:,8)),3)-mean(squeeze(aprc_land24(:,:,:,12)),3);

mean_tsw=mean(mean(tsw_exp24,4,'omitnan'),3,'omitnan');
mean_aprc_land24=mean(mean(aprc_land24,4,'omitnan'),3,'omitnan');
mean_aprc=mean(mean(aprc_exp24,4,'omitnan'),3,'omitnan');

save(fullfile(data_root, 'mean.mat'),'mean_aprc','mean_aprc_land24','mean_tsw','lon','lat','ano_aprc');
