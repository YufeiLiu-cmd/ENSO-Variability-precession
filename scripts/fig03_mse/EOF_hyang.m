%% %%%%%%%%%%%%%%%%%%%%%%%%%EOF part: do not change %%%%%%%%%%%%%%%%%%%%%%
data_root = getenv('ENSO_DATA_ROOT');
if isempty(data_root)
    data_root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end
% please provide the var as input variable. The first two dimension is
% space, the third dimension is time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n,t]=size(var);tmp=var;
if MON==1
    for k=1:12
    for i=1:m
        for j=1:n
               Monthlymean(i,j,k)=mean(squeeze(tmp(i,j,k:12:t)));
               var(i,j,k:12:t)=tmp(i,j,k:12:t)-Monthlymean(i,j,k);
        end
    end
    end
elseif YEAR==1
    var=zeros([m,n,t/12]);
    for i=1:t/12
        var(:,:,i)=mean(tmp(:,:,(i-1)*12+1:i*12),3);
    end
end
%% converct the 3d field to 2d
count=0;clear M
[m,n,t]=size(var);
for j=1:m
    for k=1:n
        if ~isnan(mean(var(j,k,:)))
             count=count+1;
             M(:,count)=var(j,k,:);
        end
    end
end
%%
if size(M,1)<size(M,2);
    M=M';shift=1;
else
   shift=0;
end
addpath(fullfile(data_root, 'external', 'pcatool')); % not bundled; see README
[eofs,pc,expvar] = caleof(M,3,1);
%% convect the 2d file to 3d
if shift==1
   tmp=pc;pc=eofs;eofs=tmp;
end
count=0;clear M
[m,n,t]=size(var);eof3d=nan([min(size(eofs)),m,n]);
for j=1:m
    for k=1:n
        if ~isnan(mean(var(j,k,:)))
             count=count+1;
             eof3d(:,j,k)=eofs(:,count);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%EOF part: do not change %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this program will calculate the first three contriution terms 
% eof3d is 
%%%%%%%%%%%%%%%%%%%%EOF part: do not change %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
