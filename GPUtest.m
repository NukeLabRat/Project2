G=gpuDevice;
CpuTime=zeros(1,30);
x=rand(1000);
x2=rand(1000);

for i=1:10
Time=tic;
Res=norm(x2-x);
CpuTime(i)=toc(Time);
end
GpuTime=zeros(1,10);
MeanCpu=mean(CpuTime)
g=gpuArray(x);
g2=gpuArray(x2);
Gres=zeros(size(g));
for i=1:10
Time=tic;
Res=norm(g2-g);
GpuTime(i)=toc(Time);

end
MeanGpu=mean(GpuTime)