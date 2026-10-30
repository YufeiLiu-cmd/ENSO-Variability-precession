%% ssta计算
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end
clearvars -except data_root; clc; close all;
%%  读数据
load(fullfile(data_root, 'thetao_allexp1.mat'),'a1','lon','lat');
load(fullfile(data_root, 'thetao_allexp2.mat'),'a2');
thetao_allexp(:,:,:,1:12)=a1;
thetao_allexp(:,:,:,13:24)=a2;
whos thetao_allexp
for num=1:1:24
   thetao_need= squeeze(thetao_allexp(:,:,:,num));
   thetao_need(thetao_need <=-99) = NaN;
   var = thetao_need;

   idx_lon = find(lon >= 120 & lon <=280);
  
   idx_var = var(idx_lon,:,:);

   % 计算月均值
   for i=1:160 
       for j = 1:10
           for month = 1:12
              ave(i,j,month)= mean(idx_var(i,j,month:12:2400),3,'omitnan');
           end
       end
   end
   % 计算ssta
   for i=1:160 
       for j = 1:10
           for t = 1:2400
              ssta(i,j,t)= idx_var(i,j,t)-ave(i,j,mod(t-1, 12) + 1);
           end
       end
   end
   % 计算std （ssta）
   stdssta = std(ssta,0,3,'omitnan');
   avestd = mean(stdssta,2,'omitnan');
   std24(:,num)=avestd;
    num=num+1
end
lona=lon(idx_lon);
%save(['meansstd24.mat'],'std24','lona');
figure(1);
pcolor(std24');shading interp;colorbar;
std24(:,25)=std24(:,1);
%%%%%%%%%%%

% 原始数据
lona = double(lona);           % 确保是 double
y = double(0:15:360);          % 确保是 double
[x, y] = meshgrid(lona, y);    % 生成原始网格
v = double(std24');           % 注意 std24 是 160×24，要转置为 24×160 对应 meshgrid

% 插值网格
lona_intep = double(120.5:0.5:279.5);
y_interp = double(0:0.5:360);
[X, Y] = meshgrid(lona_intep, y_interp);
std24_intep = interp2(x, y, v, X, Y);
figure(2);
pcolor(std24_intep');shading interp;colorbar;
save(['meansstd24_intep.mat'],'std24_intep','lona_intep','y_interp');