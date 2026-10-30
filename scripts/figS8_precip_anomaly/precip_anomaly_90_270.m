clc; clear; close all;
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end
fl=fullfile(data_root, 'aprc24')
exp =1
    for yyyy =90:180:270
name=sprintf('%03daprc.nc',yyyy);
dir=fullfile(fl,name);
file_info = ncinfo(dir);  %获取文件信息
var= ncread(dir, 'aprc')*86400; %读取变量 aprc的数据
lon=ncread(dir, 'lon');
lat=ncread(dir, 'lat');
t=ncread(dir,'time');
[m,n,t]=size(var);
% 陆地放在一起
lonneed(1:10,:)=lon(183:192,:)-360;
lonneed(11:192,:)=lon(1:182,:);
lon=lonneed;
vareed(1:10,:,:,:)=var(183:192,:,:,:);
vareed(11:192,:,:,:)=var(1:182,:,:,:);
var=vareed;
dylat=41;
dylon=300;
idx_lon=find(lon >= -20 & lon <=dylon);
idx_lat=find(lat >= -dylat & lat <=dylat);
lon_need = lon(idx_lon);
lat_need = lat(idx_lat);
[lon2d,lat2d]=meshgrid(lon_need,lat_need);
varneed= var(idx_lon,idx_lat,:,:);
[a,b,c,d] = size(varneed);
% ---- summer: June to August ----
for month = 6:8
    summer(:,:,:,month-5) = varneed(:,:,month:12:end);
end

[a,b,c,d] = size(summer);
summer_t = mean(summer, 4, 'omitnan'); 
summer_ave = mean(summer_t, 3, 'omitnan');

% ---- winter: December, January, February ----
num = 1;
for month = [12, 1, 2]
    winter(:,:,:,num) = varneed(:,:,month:12:end);
    num = num + 1;
end

[a,b,c,d] = size(winter);
winter_t = mean(winter, 4, 'omitnan'); 
winter_ave = mean(winter_t, 3, 'omitnan');

% ---- annual ----
annual_ave = mean(varneed,3,'omitnan');
annual_t = varneed;


% ---- save ----
winter_ave_exp(:,:,exp)  =winter_ave;
winter_t_exp(:,:,:,exp)  =winter_t;
summer_ave_exp(:,:,exp)  =summer_ave;
summer_t_exp(:,:,:,exp)  =summer_t;
annual_ave_exp(:,:,exp)  =annual_ave;
annual_t_exp(:,:,:,exp)  =annual_t;
exp =exp+1
    end

ano_summer = summer_ave_exp(:,:,1)-summer_ave_exp(:,:,2);
ano_winter = winter_ave_exp(:,:,1)-winter_ave_exp(:,:,2);
ano_annual = annual_ave_exp(:,:,1)-annual_ave_exp(:,:,2);

% pcolor(ano_summer');shading interp;colorbar;
% pcolor(ano_winter');shading interp;colorbar;
% pcolor(ano_annual');shading interp;colorbar;

% 双尾t 检验 summer
summer90 = reshape(summer_t_exp(:,:,:,1), [171*44, 100]);
summer270 = reshape(summer_t_exp(:,:,:,2), [171*44, 100]);
[h, sig, ci] = ttest2(summer90', summer270', 0.05);
hsummer = reshape(h, [171, 44]);
hsummer(h == 0) = nan;
% 双尾t 检验 winter
winter90 = reshape(winter_t_exp(:,:,:,1), [171*44, 100]);
winter270 = reshape(winter_t_exp(:,:,:,2), [171*44, 100]);
[h, sig, ci] = ttest2(winter90', winter270', 0.05);
hwinter = reshape(h, [171, 44]);
hwinter(h == 0) = nan;

% 双尾t 检验 annual
annual90 = reshape(annual_t_exp(:,:,:,1), [171*44, 1200]);
annual270 = reshape(annual_t_exp(:,:,:,2), [171*44, 1200]);
[h, sig, ci] = ttest2(annual90', annual270', 0.05);
hannual = reshape(h, [171, 44]);
hannual(h == 0) = nan;

save('anoaprc.mat','lon2d','lat2d','ano_summer','ano_winter','ano_annual','hannual','hsummer','hwinter');
