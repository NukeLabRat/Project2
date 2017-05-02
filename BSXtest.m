% G=gpuDevice;
Number=10;
CpuTime=zeros(1,Number);
x=rand(2000);
x2=rand(2000);

for i=1:Number
Time=tic;
Res=norm(x2-x);
CpuTime(i)=toc(Time);
end
MeanCpu=mean(CpuTime)
bsxTime=zeros(1,Number);
% g=gpuArray(x);
% g2=gpuArray(x2);
twos=ones(size(x))+1;
halves=ones(size(x))-.5;
Gres=zeros(size(x));
for i=1:Number
Time=tic;
ResBSX=bsxfun(@minus,x2,x);
NormBSX=norm(ResBSX,2);
% ResRoot=bsxfun(@power,ResSquared,halves);
bsxTime(i)=toc(Time);

end
MeanBSX=mean(bsxTime)