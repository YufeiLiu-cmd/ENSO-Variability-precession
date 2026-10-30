%% compute_nino3_enso.m
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end
% 批量处理 24 组岁差实验 tsw 数据 (每个实验 1200 个月)。
% 流程:
%   1. 循环读取每个 *tsw.nc 文件
%   2. 提取 Nino3 区域 (lat[-5,5], lon[210,270]) 区域平均 -> 原始 Nino3 序列
%   3. 去季节循环 (减月气候态) -> SSTA
%   4. 去线性长期趋势 -> 最终 Nino3 距平序列
%   5. ENSO 振幅 = SSTA 序列标准差
%   6. 卡方检验 (标准差的 95% 置信区间)
%   7. 存储为 .mat 文件
%   8. 可视化 (振幅-相位图, 顶刊风格)

clearvars -except data_root; clc; close all;

%% ---- 参数设置 ----
baseDir = fullfile(data_root, 'tsw');
latMin = -5;  latMax = 5;
lonMin = 210; lonMax = 270;
nMonthPerYear = 12;
ciLevel = 0.95;                   % 置信水平

expIDs = (0:15:345)';
nExp = numel(expIDs);
nTime = 1200;

% 坐标只读一次
firstFile = fullfile(baseDir, sprintf('%03dtsw.nc', expIDs(1)));
lon = double(ncread(firstFile, 'lon'));
lat = double(ncread(firstFile, 'lat'));
latIdx = find(lat >= latMin & lat <= latMax);
lonIdx = find(lon >= lonMin & lon <= lonMax);

% 预分配结果
nino3_all = zeros(nTime, nExp);
ssta_all  = zeros(nTime, nExp);

%% ---- 循环读文件计算 Nino3 ----
for k = 1:nExp
    yyyy  = expIDs(k);
    fname = sprintf('%03dtsw.nc', yyyy);
    fpath = fullfile(baseDir, fname);
    tsw = double(ncread(fpath, 'tsw'));   % [lon x lat x time]

    % Nino3 区域平均 -> 列向量
    nino3_temp = zeros(nTime, 1);
    for t = 1:nTime
        nino3_temp(t) = mean(tsw(lonIdx, latIdx, t), 'all');
    end
    nino3_all(:, k) = nino3_temp;

    % 去季节循环 -> SSTA
    tsMat = reshape(nino3_temp, nMonthPerYear, []);   % [12 x 100]
    clim  = mean(tsMat, 2);
    anom  = tsMat - repmat(clim, 1, 100);
    ssta  = anom(:);

    % 去长期趋势
    ssta = detrend(ssta, 1);
    ssta_all(:, k) = ssta;

    fprintf('EXP %03d done\n', yyyy);
end

%% ---- 计算标准差 (ENSO 振幅) 与卡方 95% 置信区间 ----
enso_amp    = zeros(nExp, 1);     % 每组标准差
chi2_ci_std = zeros(nExp, 2);     % 标准差 95% 置信区间 [下限, 上限]

alpha = 1 - ciLevel;
for k = 1:nExp
    x  = ssta_all(:, k);
    n  = numel(x);
    s2 = var(x, 0);               % 样本方差 (ddof=1)
    df = n - 1;

    enso_amp(k) = sqrt(s2);

    % 卡方分布分位点 (优先用工具箱, 无则回退 gaminv)
    chiLow  = chi2invLocal(alpha/2, df);
    chiHigh = chi2invLocal(1 - alpha/2, df);

    % 标准差置信区间: sqrt(df*s2/chiHigh) ~ sqrt(df*s2/chiLow)
    chi2_ci_std(k, 1) = sqrt(df * s2 / chiHigh);
    chi2_ci_std(k, 2) = sqrt(df * s2 / chiLow);
end

%% ---- 存储 .mat ----
matPath = fullfile(baseDir, 'nino3_enso_chi2_results.mat');
save(matPath, 'expIDs', 'nino3_all', 'ssta_all', ...
    'enso_amp', 'chi2_ci_std', 'ciLevel', ...
    'latMin', 'latMax', 'lonMin', 'lonMax', '-v7');
fprintf('已保存: %s\n', matPath);

%% ==================== 顶刊风格可视化 ====================
% x 轴: 岁差相位 (000 -> 345, 间隔 15度)
% 左 y 轴: ENSO 振幅 (std of SSTA, 0.5~1.4);  右 y 轴: 海洋层结 (0.4~1.1)
errLow  = enso_amp - chi2_ci_std(:, 1);       % 下误差
errHigh = chi2_ci_std(:, 2) - enso_amp;       % 上误差

% 配色
col_main = [0.80 0.15 0.15];      % 红色主线 (ENSO 振幅)
col_ci   = [0.90 0.55 0.55];      % 浅红置信区间
col_fill = [0.97 0.88 0.88];      % 背景填充
col_str  = [0.13 0.40 0.62];      % 蓝色 (海洋层结)

% 字体基准 (整体较之前 +2 号)
FS_axis  = 13;
FS_label = 16;
FS_title = 17;
FS_leg   = 13;

fig = figure('Color', 'w', 'Units', 'centimeters', ...
             'Position', [3 3 22 15]);
ax = axes('Parent', fig);

%% ---- 左 y 轴: ENSO 振幅 ----
yyaxis(ax, 'left'); hold(ax, 'on');

% 1) 卡方 95% CI 半透明色带
xBand = [expIDs; flipud(expIDs)];
yBand = [chi2_ci_std(:, 1); flipud(chi2_ci_std(:, 2))];
hBand = fill(xBand, yBand, col_fill, 'EdgeColor', 'none', ...
     'FaceAlpha', 0.7, 'Parent', ax);

% 2) 红色主线 + 垂直误差棒 (卡方 95% CI)
hAmp = errorbar(ax, expIDs, enso_amp, errLow, errHigh, ...
    '-o', 'Color', col_main, 'LineWidth', 1.8, ...
    'MarkerSize', 6, 'MarkerFaceColor', col_main, ...
    'MarkerEdgeColor', 'w', 'CapSize', 6);

set(ax, 'YColor', col_main, ...
        'YLim', [0.5, 1.4], 'YTick', 0.5:0.1:1.4);
ylabel(ax, 'ENSO AMPLITUDE (K)', ...
       'FontSize', FS_label, 'FontWeight', 'bold', 'Color', col_main);

%% ---- 右 y 轴: 海洋层结 ----
yyaxis(ax, 'right');
hStr = [];
if exist('stratification', 'var') == 1
    % 若已定义 stratification (24x1), 则画蓝色实线
    hStr = plot(ax, expIDs, stratification(:), '-s', ...
        'Color', col_str, 'LineWidth', 1.8, ...
        'MarkerSize', 6, 'MarkerFaceColor', col_str, ...
        'MarkerEdgeColor', 'w');
end
set(ax, 'YColor', col_str, ...
        'YLim', [0.4, 1.1], 'YTick', 0.4:0.1:1.1);
ylabel(ax, 'OCEAN STRATIFICATION', ...
       'FontSize', FS_label, 'FontWeight', 'bold', 'Color', col_str);

%% ---- x 轴 / 全局样式 ----
set(ax, 'XLim', [-10, 355], 'XTick', 0:15:345, ...
        'FontName', 'Helvetica', 'FontSize', FS_axis, ...
        'LineWidth', 1.1, 'TickDir', 'in', ...   % 刻度线全部朝内
        'Box', 'on', 'Layer', 'top');
xtickangle(ax, 90);
% x 轴标签补零成 000, 015, ...
ax.XTickLabel = arrayfun(@(v) sprintf('%03d', v), (0:15:345)', 'uni', 0);

grid(ax, 'on');
ax.GridLineStyle = '-';
ax.GridColor = [0.85 0.85 0.85];
ax.GridAlpha = 0.6;

xlabel(ax, 'Precession phase (\circ)', ...
       'FontSize', FS_label, 'FontWeight', 'bold', 'Interpreter', 'tex');
title(ax, 'ENSO amplitude', ...
      'FontSize', FS_title, 'FontWeight', 'bold');

%% ---- 图注: 放在图内右下角 ----
if isempty(hStr)
    legItems  = [hAmp, hBand];
    legLabels = {sprintf('ENSO amplitude (Ni\tildeno3 SSTA std)'), ...
                 sprintf('\\chi^2 %d%% CI', round(ciLevel*100))};
else
    legItems  = [hAmp, hBand, hStr];
    legLabels = {sprintf('ENSO amplitude (Ni\tildeno3 SSTA std)'), ...
                 sprintf('\\chi^2 %d%% CI', round(ciLevel*100)), ...
                 'Ocean stratification'};
end
lg = legend(legItems, legLabels, 'Interpreter', 'tex', ...
    'FontSize', FS_leg, 'Location', 'southeast', 'Box', 'on');
lg.Color = [1 1 1];
lg.EdgeColor = [0.7 0.7 0.7];

hold(ax, 'off');

% ---- 导出高分辨率图 ----
pngPath = fullfile(baseDir, 'enso_amplitude_vs_precession.png');
exportgraphics(fig, pngPath, 'Resolution', 300);
pdfPath = fullfile(baseDir, 'enso_amplitude_vs_precession.pdf');
exportgraphics(fig, pdfPath, 'ContentType', 'vector');
fprintf('可视化完成: %s\n', pngPath);

%% ================= 局部函数 =================
function x = chi2invLocal(p, df)
% 卡方分布逆累计分布函数 (优先工具箱, 否则回退)
    if exist('chi2inv', 'file') == 2
        x = chi2inv(p, df); return;
    end
    if exist('gaminv', 'file') == 2
        x = gaminv(p, df/2, 2); return;   % 卡方(df) = Gamma(df/2, 2)
    end
    % 无工具箱: 基于 gammainc 的二分法
    lo = 0; hi = max(10*df, 100);
    while gammainc(hi/2, df/2) < p, hi = hi*2; end
    for it = 1:200
        mid = 0.5*(lo+hi);
        if gammainc(mid/2, df/2) < p, lo = mid; else, hi = mid; end
    end
    x = 0.5*(lo+hi);
end
