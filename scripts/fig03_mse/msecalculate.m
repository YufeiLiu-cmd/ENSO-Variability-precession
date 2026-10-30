clc;clear;close all;
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end
folder1=fullfile(data_root, 'MSEgeopoth');
folder2=fullfile(data_root, 'aps');
num=1
for yyyy = 0:15:345
name1 = sprintf('%03dforMSE.nc', yyyy);  % 生成文件名
name2 = sprintf('%03dapsforMSE.nc', yyyy);  % 生成文件名
dir1= fullfile(folder1, name1);
file_info = ncinfo(dir1);  %获取ts文件信息
dir2=fullfile(folder2, name2);
lon = ncread(dir1, 'lon');
lat = ncread(dir1, 'lat');
aps = ncread(dir2, 'aps');  % 单位：Pa
st = ncread(dir1, 'st');  
q = ncread(dir1, 'q');
Z = ncread(dir1, 'geopoth');
plev = ncread(dir1, 'plev'); 
%% calculate MSE  S= Cp.st+g.Z+Lv.q
Cp=1005 ;  %J/kg.K
g=9.8;
Lv=2.26*10^6 ; %J/kg
MSE_4d = Cp*st+g*Z+Lv*q;  %J/kg
[a,b,c,d]=size(MSE_4d);

%% monthly
% for month = 1:12
% aps_monthly(:,:,month) = mean(aps(:,:,month:12:d),3,'omitnan');
% MSE_monthly(:,:,:,month) = mean(MSE_4d(:,:,:,month:12:d), 4, 'omitnan');
% end
% [a1,b1,c1,d1]= size(MSE_monthly);

%% total level ave
% 取每层中间值 
for i=1:16
    next=i+1;
    plev_mid(i)=plev(i)-(plev(i)-plev(next))/2;
end
plev_mid(17)=500;
for month = 1:d
    MSE_var=squeeze(MSE_4d(:,:,:,month));
    aps_var=squeeze(aps(:,:,month));
    % 第一层
    p_diff=aps_var-plev_mid(1);
    MSE_total(:,:,1)=p_diff.*MSE_var(:,:,1);
    % 第2-17层
    for i=1:16
      MSE_total(:,:,i+1)=(plev_mid(i)-plev_mid(i+1))*MSE_var(:,:,i+1); 
    end
  MSE_ave=sum(MSE_total,3,'omitnan')./(aps_var);
  MSE_ave_monthly(:,:,month)=MSE_ave/1000;  % 换算成KJ/kg
end
% [x,y]=meshgrid(lon,lat);
% figure(1)
% pcolor(x,y,squeeze(MSE_ave_monthly(:,:,2)'));shading interp; colorbar;
% figure(2)
% pcolor(x,y,squeeze(MSE_ave_monthly(:,:,8)'));shading interp; colorbar;
% stop
num=num+1
save_name = sprintf('MSE_ave_month_%03d.mat', yyyy);
save(save_name, 'lat', 'lon', 'MSE_ave_monthly');
end

