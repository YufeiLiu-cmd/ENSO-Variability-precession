
;--------------------------------------------------------------------------------------
;循环读文件
dirout=getenv("ENSO_DATA_ROOT")+"/walker/"
name= (/"105","285"/)
dir=getenv("ENSO_DATA_ROOT")+"/walker/"
diraps=getenv("ENSO_DATA_ROOT")+"/walker/"

do k=0

a = 6.37122e06 ;  radius of the earth
PI=3.1415926
g=9.8
lat1=-10
lat2=10
lon1=-180
lon2=180
psitime = new((/1200, 17, 192/), "float")  ; 初始化

  ;读取数据
fo       =addfile(dir+name(k)+".nc","r")
u       = fo->u
v       = fo->v
t       = fo->time

f1       =addfile(diraps+"select"+name(k)+".nc", "r")
Slp     = f1->aps

do  i = 0,1199

U=u(i,::-1,::-1,:)
;U=dim_avg_n_Wrap(u(:,::-1,::-1,:),0) ;  var131是U,var132是V,var134是地面气压
U=lonFlip(U) ;将经度从0~360度变成-180~180
;printVarSummary(U)

V=v(i,::-1,::-1,:)
;V=dim_avg_n_Wrap(v(:,::-1,::-1,:),0) ;对第1维求平均，并把第二维第三维数据反向
V=lonFlip(V)

PS_tmp=Slp(i,:,:)
PS_tmp=lonFlip(PS_tmp)
PS=PS_tmp({lat1:lat2},{lon1:lon2}) ;选取lat1~lat2,lon1~lon2范围内的数据

DIV=uv2dvG_Wrap(U,V) ;计算散度
DIV=DIV(:,::-1,:)   ;将第二维反向
uvd = dv2uvG_Wrap(DIV)  ;计算发散风分量,https://www.ncl.ucar.edu/Document/Functions/Built-in/dv2uvG-1.shtml
uu=uvd(0,:,{lat1:lat2},{lon1:lon2}) ;发散风u分量
vv=uvd(1,:,{lat1:lat2},{lon1:lon2}) ;发散风v分量

U_final=uu(plev|:,lon|:,lat|:) ;发散风u分量
U_final!1="lat"  ;将lon维度改成lat
U_final!2="lon"  ;将lat维度改成lon
PS_final=PS(lon|:,lat|:) ;地面气压维度（lat,lon）互换
PS_final!0="lat" ;将lon维度改成lat
PS_final!1="lon" ;将lat维度改成lon
lat=PS_final&lat ;读取新的lat变量，实际上是原来的lon
p=U_final&plev   ;读取气压层变量
psi=zonal_mpsi_Wrap(U_final,lat*0,p,PS_final) ; multiply 100 because pres is in mb
psi=psi/1e9      ;计算流函数
psi@units="svp"

psi!1   = "lon"
;printVarSummary(psi)
;print("Dimensions of psi: " + dimsizes(psi))
psitime(i,:,:) = psi


delete(U)
delete(V)
delete(PS)
delete(DIV)
delete(uvd)
delete(uu)
delete(vv)
delete(p)
delete(U_final)
delete(PS_final)
delete(psi)
end do
psitime!0 = "time"
; print("Final psitime values: ")
; print(psitime)
f       = addfile(dirout+name(k)+"psi.nc","c")
f->psitime  = psitime

end do 


