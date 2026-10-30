clearvars -except data_root; clc; close all;
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end

%% 读数据文件
num = 1;

for yyyy = 0:15:345

    dir = fullfile(data_root, sprintf('%03d', yyyy), 'fesom_regular_thetao_y3300-3399.nc');
    lon=ncread(dir,'lon');
    lat=ncread(dir,'lat');
    time=ncread(dir,'time');
    t=length(time);
    thetao = ncread(dir,'thetao');
    thetao = squeeze(thetao);

 
% stop  % removed in the public copy; it halted the script before the centroid loop





%% 读数据
clear lon lat lat2d lon2d ap lat1 lon1 idx_lon idx_lat
    
    ap=ncread(dir,'aprc');
    ap=ap*86400  ;%% 转换单位


    %  选区域并计算月份 
  idx_lon = find(lon >=30 & lon <= 240);
  idx_lat= find(lat >= -30 & lat <= 30);
  ap1=ap(idx_lon,idx_lat,:);
  [m,n,t]=size(ap1);
  
    %   求每个时间步长的质心

for t =1:1200
  var=squeeze(ap1(:,:,t));
  [lat1,lon1] = meshgrid(lat(idx_lat),lon(idx_lon));%网格
  % 选点
  point = [];
lon_needpoint=lon(idx_lon);
lat_needpoint=lat(idx_lat);
[m, n] = size(var);
for j = 1:n
   for i = 1:m
       if var(i, j) >= 4.93
            % 将满足条件的点的经纬度添加到 point 变量中
            point(end+1, 1) = lon_needpoint(i);
            point(end, 2) = lat_needpoint(j);
            point(end,3)=var(i, j);
        end
    end
end
%% 求质心
%经度坐标
centroid_lon = sum(point(:, 1) .* point(:, 3)) / sum(point(:, 3));
%纬度坐标
%centroid_lat = sum(point(:, 2) .* point(:, 3))/sum(point(:, 3));
barycenter(t,1)=centroid_lon;
%barycenter(t,2)=centroid_lat;

end
barycenter24(num,:) = barycenter;
num=num+1;
end
% keyboard;  % removed in the public copy; interactive breakpoint
barycentersolo=mean(barycenter24,2);
error24 = 1.96*(std(barycenter24, 0, 2)./ sqrt(1200));

x = 0:15:345; % 假设是 0 到 345 的等间隔横坐标
figure;
errorbar(x, barycentersolo, error24, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
barycenter=barycentersolo;
save(['exp24aprczhixin.mat'],'barycenter','error24');

