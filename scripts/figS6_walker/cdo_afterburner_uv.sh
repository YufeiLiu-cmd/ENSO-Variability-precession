#!/bin/bash
set -x  # 启用调试模式
# 循环生成文件夹名并执行命令

for i in $(seq -f "%03g" 105 180 285); do 
# echam6 manual: https://www.cen.uni-hamburg.de/en/icdc/data/atmosphere/docs-atmo/echam6-userguide.pdf
# echam6 afterbuner: https://code.mpimet.mpg.de/projects/afterburner/wiki

# ***** time extents *****
module load cdo
cat << EOF > namelist.after
  &SELECT
    CODE = 131,132
    TYPE = 70,
    FORMAT = 2
  &END
EOF

out="${ENSO_DATA_ROOT:-./data}/walker"
mkdir -p "$out"
cdo after "${i}.nc" "${out}/${i}.nc" < namelist.after

done
# svo  138 1/s  vorticity
# sd   155 1/s  divergence
# lsp  152 -    nat. logarithm of surface
# geosp129 m2/s2 surface geopotential
# st   130 K    temperature
# q    133 -    specific humidity
# aps  134 Pa   surface pressure
# aprc 143 kg/(m2s)   convective precipitation
# slm  172 -             land sea mask
# tsw  103 K    surface temperature of water      
# 156  geopotential height
# 151  mean sea level pressure
