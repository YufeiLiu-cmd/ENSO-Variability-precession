clearvars -except data_root; clc; close all;
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end
% 文件夹路径
folder =fullfile(data_root, 'psi')
num=1;
    for yyyy = 105:180:285
name = sprintf('%03dpsi.nc', yyyy);  % 生成文件名
file= fullfile(folder, name);
file_info = ncinfo(file);  %获取ts文件信息
psi=ncread(file,'psitime');
lon=ncread(file,'lon');
plev=ncread(file,'plev')./100;

for i=1:size(lon,1);
if lon(i)<0;
lon(i)=lon(i)+360;
end
end
lon_1(1:96,1)=lon(97:192,1);
lon_1(97:192,1)=lon(1:96,1);
lon=lon_1;
clear lon_1;
tem(1:96,:,:)=psi(97:192,:,:);
tem(97:192,:,:)=psi(1:96,:,:);
psi=tem;

% 求平均
psi_ave = mean(psi,3,'omitnan');
psi_ri(:,:,:,num) = psi;
psi_ave2(:,:,num)=psi_ave;
num=num+1
    end


    % 选区域
    idx_lon = find(lon >= 40 & lon <= 280);
    lon_need = lon(idx_lon);
    psi_ave =psi_ave2(idx_lon,:,:);
    psi_ri =psi_ri (idx_lon,:,:,:);
    [a,b,c,d]=size(psi_ri)
  
    % ano 
    diff = squeeze(psi_ave(:,:,1))-squeeze(psi_ave(:,:,2));
    % t检验
psi105=reshape(squeeze(psi_ri(:,:,:,1)),[128*17,1200]);
psi285=reshape(squeeze(psi_ri(:,:,:,2)),[128*17,1200]);

% 双样本t检验
[h,sig,ci]=ttest2(psi105',psi285',0.05);
h=reshape(h,[128,17]);
h(h==0)=nan;

save("psi.mat",'lon_need','h',"diff",'psi_ave','plev');
