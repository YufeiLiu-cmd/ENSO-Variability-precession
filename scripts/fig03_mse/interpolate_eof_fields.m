clearvars -except data_root; clc; close all;
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end
load(fullfile(data_root, 'mse_eof_winter.mat'),'lon2d','lat2d','winter_ipwp_90','winter_ipwp_270','time_serise','eof_result','winter_aprc');
% 原始数据
lon = double(lon2d(1,:));           % 确保是 double
lat = double(lat2d(:,1));          % 确保是 double
[x, y] = meshgrid(lon,lat);    % 生成原始网格
v = double(eof_result);           % 注意 std24 是 160×24，要转置为 24×160 对应 meshgrid

% 插值网格
lon_intep = double(lon(1):0.1:lon(end));
lat_interp = double(lat(1):-0.1:lat(end));
[X, Y] = meshgrid(lon_intep,lat_interp);
V = interp2(x, y, v', X, Y);
lon2d = X;
lat2d = Y;
eof_result = V;
save(fullfile(data_root, 'mse_eof_winternew.mat'),'lon2d','lat2d','winter_ipwp_90','winter_ipwp_270','time_serise','eof_result','winter_aprc');
%%
clearvars -except data_root; clc; close all;
load(fullfile(data_root, 'mse_eof_winter.mat'),'lon2d','lat2d');
load(fullfile(data_root, 'mse_eof_summer.mat'),'time_serise','eof_result','summer_aprc');
% 原始数据
lon = double(lon2d(1,:));           % 确保是 double
lat = double(lat2d(:,1));          % 确保是 double
[x, y] = meshgrid(lon,lat);    % 生成原始网格
v = double(eof_result);           % 注意 std24 是 160×24，要转置为 24×160 对应 meshgrid

% 插值网格
lon_intep = double(lon(1):0.1:lon(end));
lat_interp = double(lat(1):-0.1:lat(end));
[X, Y] = meshgrid(lon_intep,lat_interp);
V = interp2(x, y, v', X, Y);
lon2d = X;
lat2d = Y;
eof_result = V;
save(fullfile(data_root, 'mse_eof_summernew.mat'),'time_serise','eof_result','summer_aprc');

%%
clearvars -except data_root; clc; close all;
load(fullfile(data_root, 'mse_eof_winter.mat'),'lon2d','lat2d');
load(fullfile(data_root, 'mse_eof_annual.mat'),'time_serise','eof_result','annual_aprc');
% 原始数据
lon = double(lon2d(1,:));           % 确保是 double
lat = double(lat2d(:,1));          % 确保是 double
[x, y] = meshgrid(lon,lat);    % 生成原始网格
v = double(eof_result);           % 注意 std24 是 160×24，要转置为 24×160 对应 meshgrid

% 插值网格
lon_intep = double(lon(1):0.1:lon(end));
lat_interp = double(lat(1):-0.1:lat(end));
[X, Y] = meshgrid(lon_intep,lat_interp);
V = interp2(x, y, v', X, Y);
lon2d = X;
lat2d = Y;
eof_result = V;
save(fullfile(data_root, 'mse_eof_annualnew.mat'),'time_serise','eof_result','annual_aprc');