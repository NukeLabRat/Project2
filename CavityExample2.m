%% Shubham Maurya
%% Indian Institute of Space Science and Technology
clear all
clc
Re=100; % Reynolds Number
L=1; % Length of edges
m=11; % Number of grid points
n=11;
dx=L/(m+1); % Step size in length
dy=L/(n+1);
dt=0.01; %time step size
T=5; % end time
N=T/dt; % No. of time steps
%% Boundary conditions %%
%% Left wall
u1=zeros(1,m);
v1=zeros(1,m);
s1=zeros(1,m);

%% Bottom wall
u2=zeros(1,m);
v2=zeros(1,m);
s2=zeros(1,m);

%% Right wall
u3=zeros(1,m);
v3=zeros(1,m);
s3=zeros(1,m);

%% Top wall
u4=ones(1,m);
v4=zeros(1,m);
s4=zeros(1,m);

%% Initial conditions %%
u0=zeros(n,m);
v0=zeros(n,m);
psi0=zeros(n,m);
psi=zeros(n,m);
f=zeros(n,m);
w1=zeros(n,m);
wleft=zeros(1,n);
wdown=zeros(1,m);
wright=zeros(1,n);
wtop=zeros(1,m);

%% Initial vorticity field (at t=0, entire fluid is at rest)%%
w=zeros(n,m);
for i=1:n
    for j=1:m
        if (i==1)
            if(j==1)
                w(i,j)=-(v0(i,j+1)-v1(i)+u2(j)-u0(i+1,j))/(2*dx);
            else if(j==n)
                    w(i,j)=-(v3(i)-v0(i,j-1)+u2(j)-u0(i+1,j))/(2*dx);
                else
                    w(i,j)=-(v0(i,j+1)-v0(i,j-1)+u2(j)-u0(i+1,j))/(2*dx);
                end
            end
        else if (i==m)
                if(j==1)
                    w(i,j)=-(v0(i,j+1)-v1(i)+u0(i-1,j)-u4(j))/(2*dx);
                else if (j==n)
                        w(i,j)=-(v3(i)-v0(i,j-1)+u0(i-1,j)-u4(j))/(2*dx);
                    else
                        w(i,j)=-(v0(i,j+1)-v0(i,j-1)+u0(i-1,j)-u4(j))/(2*dx);
                    end
                end
            else
                if (j==1)
                    w(i,j)=-(v0(i,j+1)-v1(i)+u0(i-1,j)-u0(i+1,j))/(2*dx);
                else if (j==n)
                        w(i,j)=-(v3(i)-v0(i,j-1)+u0(i-1,j)-u0(i+1,j))/(2*dx);
                    else
                        w(i,j)=-(v0(i,j+1)-v0(i,j-1)+u0(i-1,j)-u0(i+1,j))/(2*dx);
                    end
                end
            end
        end
    end
end

%% Initial stream function %%
res=1;
iter=0;
while(res>0.0001)
    for i=1:n
        for j=1:m
            if (i==1)
                if (j==1)
                    psi(i,j)=0.25*(psi0(i+1,j)+psi0(i,j+1)+s2(j)+s1(i)+w(i,j)*dx^2);
                else if (j==m)
                        psi(i,j)=0.25*(psi0(i+1,j)+psi0(i,j-1)+s2(j)+s3(i)+w(i,j)*dx^2);
                    else
                        psi(i,j)=0.25*(psi0(i+1,j)+psi0(i,j+1)+psi0(i,j-1)+s2(j)+w(i,j)*dx^2);
                    end
                end
            else if (i==n)
                    if (j==1)
                        psi(i,j)=0.25*(psi0(i-1,j)+s4(j)+psi0(i,j+1)+s1(i)+w(i,j)*dx^2);
                    else if (j==m)
                            psi(i,j)=0.25*(s3(i)+s4(j)+psi0(i,j-1)+psi0(i-1,j)+w(i,j)*dx^2);
                        else
                            psi(i,j)=0.25*(psi0(i-1,j)+s4(j)+psi0(i,j-1)+psi0(i,j+1)+w(i,j)*dx^2);
                        end
                    end
                else
                    if (j==1)
                        psi(i,j)=0.25*(psi0(i+1,j)+psi0(i,j+1)+psi0(i-1,j)+s1(i)+w(i,j)*dx^2);
                    else if (j==m)
                            psi(i,j)=0.25*(s3(i)+psi0(i+1,j)+psi0(i,j-1)+psi0(i-1,j)+w(i,j)*dx^2);
                        else
                            psi(i,j)=0.25*(psi0(i+1,j)+psi0(i,j+1)+psi0(i,j-1)+psi0(i-1,j)+w(i,j)*dx^2);
                        end
                    end
                end
            end
        end
    end
    h=0;
    res=zeros(1,m*n);
    for i=1:n
        for j=1:m
            k=j+h;
            if (i==1)
                if (j==1)
                    res(k)=psi(i,j)-0.25*(psi(i+1,j)+psi(i,j+1)+s2(j)+s1(i)+w(i,j)*dx^2);
                else if (j==m)
                        res(k)=psi(i,j)-0.25*(psi(i+1,j)+psi(i,j-1)+s2(j)+s3(i)+w(i,j)*dx^2);
                    else
                        res(k)=psi(i,j)-0.25*(psi(i+1,j)+psi(i,j+1)+psi(i,j-1)+s2(j)+w(i,j)*dx^2);
                    end
                end
            else if (i==n)
                    if (j==1)
                        res(k)=psi(i,j)-0.25*(psi(i-1,j)+s4(j)+psi(i,j+1)+s1(i)+w(i,j)*dx^2);
                    else if (j==m)
                            res(k)=psi(i,j)-0.25*(s3(i)+s4(j)+psi(i,j-1)+psi(i-1,j)+w(i,j)*dx^2);
                        else
                            res(k)=psi(i,j)-0.25*(psi(i-1,j)+s4(j)+psi(i,j-1)+psi(i,j+1)+w(i,j)*dx^2);
                        end
                    end
                else
                    if (j==1)
                        res(k)=psi(i,j)-0.25*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+s1(i)+w(i,j)*dx^2);
                    else if (j==m)
                            res(k)=psi(i,j)-0.25*(s3(i)+psi(i+1,j)+psi(i,j-1)+psi(i-1,j)+w(i,j)*dx^2);
                        else
                            res(k)=psi(i,j)-0.25*(psi(i+1,j)+psi(i,j+1)+psi(i,j-1)+psi(i-1,j)+w(i,j)*dx^2);
                        end
                    end
                end
            end
        end
        h=m*i;
    end
    res=rms(res);
    psi0=psi;
    iter=iter+1;
end

%% Iterations start

time_step=0;
figure(1)
while(time_step<N)
    %% Boundary conditions for vorticity %%
    for i=1:n
        wleft(i)=-2*v1(i)/dx+2*(s1(i)-psi(i,1))/dx^2;
        wright(i)=2*v3(i)/dx+2*(s3(i)-psi(i,m))/dx^2;
    end
    for i=1:m
        wdown(i)=2*u2(i)/dy+2*(s2(i)-psi(1,i))/dy^2;
        wtop(i)=-2*u4(i)/dy+2*(s4(i)-psi(n,i))/dy^2;
    end
    %% Vorticity equation solver %%
    for i=1:n
        for j=1:m
            if i==1
                if j==1
                    up=0.25*((psi(i+1,j)-s2(j))+abs((psi(i+1,j)-s2(j))))/dy;
                    um=0.25*((psi(i+1,j)-s2(j))-abs((psi(i+1,j)-s2(j))))/dy;
                    vp=0.25*((-psi(i,j+1)+s1(i))+abs((psi(i,j+1)-s1(i))))/dx;
                    vm=0.25*((-psi(i,j+1)+s1(i))-abs((psi(i,j+1)-s1(i))))/dx;
                    f(i,j)=(w(i+1,j)+wdown(j)-4*w(i,j)+w(i,j+1)+wleft(i))/(Re*dx^2)-((psi(i+1,j)-s2(j))*(w(i,j+1)-wleft(i))+(s1(i)-psi(i,j+1))*(w(i+1,j)-wdown(j)))/(4*dx^2)+0.5*(up+um)*(w(i,j+1)+wleft(i)-2*w(i,j))/dx+0.5*(vp+vm)*(w(i+1,j)+wdown(j)-2*w(i,j))/dy;
                else if j==m
                        up=0.25*((psi(i+1,j)-s2(j))+abs((psi(i+1,j)-s2(j))))/dy;
                        um=0.25*((psi(i+1,j)-s2(j))-abs((psi(i+1,j)-s2(j))))/dy;
                        vp=0.25*((-s3(i)+psi(i,j-1))+abs((s3(i)-psi(i,j-1))))/dx;
                        vm=0.25*((-s3(i)+psi(i,j-1))-abs((s3(i)-psi(i,j-1))))/dx;
                        f(i,j)=(w(i+1,j)+wdown(j)-4*w(i,j)+wright(i)+w(i,j-1))/(Re*dx^2)-((psi(i+1,j)-s2(j))*(wright(i)-w(i,j-1))+(psi(i,j-1)-s3(i))*(w(i+1,j)-wdown(j)))/(4*dx^2)+0.5*(up+um)*(wright(i)+w(i,j-1)-2*w(i,j))/dx+0.5*(vp+vm)*(w(i+1,j)+wdown(j)-2*w(i,j))/dy;
                    else
                        up=0.25*((psi(i+1,j)-s2(j))+abs((psi(i+1,j)-s2(j))))/dy;
                        um=0.25*((psi(i+1,j)-s2(j))-abs((psi(i+1,j)-s2(j))))/dy;
                        vp=0.25*((-psi(i,j+1)+psi(i,j-1))+abs((psi(i,j+1)-psi(i,j-1))))/dx;
                        vm=0.25*((-psi(i,j+1)+psi(i,j-1))-abs((psi(i,j+1)-psi(i,j-1))))/dx;
                        f(i,j)=(w(i+1,j)+wdown(j)-4*w(i,j)+w(i,j+1)+w(i,j-1))/(Re*dx^2)-((psi(i+1,j)-s2(j))*(w(i,j+1)-w(i,j-1))+(psi(i,j-1)-psi(i,j+1))*(w(i+1,j)-wdown(j)))/(4*dx^2)+0.5*(up+um)*(w(i,j+1)+w(i,j-1)-2*w(i,j))/dx+0.5*(vp+vm)*(w(i+1,j)+wdown(j)-2*w(i,j))/dy;
                    end
                end
            else if i==n
                    if (j==1)
                        up=0.25*((s4(j)-psi(i-1,j))+abs((s4(j)-psi(i-1,j))))/dy;
                        um=0.25*((s4(j)-psi(i-1,j))-abs((s4(j)-psi(i-1,j))))/dy;
                        vp=0.25*((-psi(i,j+1)+s1(i))+abs((psi(i,j+1)-s1(i))))/dx;
                        vm=0.25*((-psi(i,j+1)+s1(i))-abs((psi(i,j+1)-s1(i))))/dx;
                        f(i,j)=(wtop(j)+w(i-1,j)-4*w(i,j)+w(i,j+1)+wleft(i))/(Re*dx^2)-((s4(j)-psi(i-1,j))*(w(i,j+1)-wleft(i))+(s1(i)-psi(i,j+1))*(wtop(j)-w(i-1,j)))/(4*dx^2)+0.5*(up+um)*(w(i,j+1)+wleft(i)-2*w(i,j))/dx+0.5*(vp+vm)*(wtop(j)+w(i-1,j)-2*w(i,j))/dy;
                    else if j==m
                            up=0.25*((s4(j)-psi(i-1,j))+abs((s4(j)-psi(i-1,j))))/dy;
                            um=0.25*((s4(j)-psi(i-1,j))-abs((s4(j)-psi(i-1,j))))/dy;
                            vp=0.25*((-s3(i)+psi(i,j-1))+abs((s3(i)-psi(i,j-1))))/dx;
                            vm=0.25*((-s3(i)+psi(i,j-1))-abs((s3(i)-psi(i,j-1))))/dx;
                            f(i,j)=(wtop(j)+w(i-1,j)-4*w(i,j)+wright(i)+w(i,j-1))/(Re*dx^2)-((s4(j)-psi(i-1,j))*(wright(i)-w(i,j-1))+(psi(i,j-1)-s3(i))*(wtop(j)-w(i-1,j)))/(4*dx^2)+0.5*(up+um)*(wright(i)+w(i,j-1)-2*w(i,j))/dx+0.5*(vp+vm)*(wtop(j)+w(i-1,j)-2*w(i,j))/dy;
                        else
                            up=0.25*((s4(j)-psi(i-1,j))+abs((s4(j)-psi(i-1,j))))/dy;
                            um=0.25*((s4(j)-psi(i-1,j))-abs((s4(j)-psi(i-1,j))))/dy;
                            vp=0.25*((-psi(i,j+1)+psi(i,j-1))+abs((psi(i,j+1)-psi(i,j-1))))/dx;
                            vm=0.25*((-psi(i,j+1)+psi(i,j-1))-abs((psi(i,j+1)-psi(i,j-1))))/dx;
                            f(i,j)=(wtop(j)+w(i-1,j)-4*w(i,j)+w(i,j+1)+w(i,j-1))/(Re*dx^2)-((s4(j)-psi(i-1,j))*(w(i,j+1)-w(i,j-1))+(psi(i,j-1)-psi(i,j+1))*(wtop(j)-w(i-1,j)))/(4*dx^2)+0.5*(up+um)*(w(i,j+1)+w(i,j-1)-2*w(i,j))/dx+0.5*(vp+vm)*(wtop(j)+w(i-1,j)-2*w(i,j))/dy;
                        end
                    end
                else
                    if j==1
                        up=0.25*((psi(i+1,j)-psi(i-1,j))+abs((psi(i+1,j)-psi(i-1,j))))/dy;
                        um=0.25*((psi(i+1,j)-psi(i-1,j))-abs((psi(i+1,j)-psi(i-1,j))))/dy;
                        vp=0.25*((-psi(i,j+1)+s1(i))+abs((psi(i,j+1)-s1(i))))/dx;
                        vm=0.25*((-psi(i,j+1)+s1(i))-abs((psi(i,j+1)-s1(i))))/dx;
                        f(i,j)=(w(i+1,j)+w(i-1,j)-4*w(i,j)+w(i,j+1)+wleft(i))/(Re*dx^2)-((psi(i+1,j)-psi(i-1,j))*(w(i,j+1)-wleft(i))-(s1(i)+psi(i,j+1))*(w(i+1,j)-w(i-1,j)))/(4*dx^2)+0.5*(up+um)*(w(i,j+1)+wleft(i)-2*w(i,j))/dx+0.5*(vp+vm)*(w(i+1,j)+w(i-1,j)-2*w(i,j))/dy;
                    else if j==m
                            up=0.25*((psi(i+1,j)-psi(i-1,j))+abs((psi(i+1,j)-psi(i-1,j))))/dy;
                            um=0.25*((psi(i+1,j)-psi(i-1,j))-abs((psi(i+1,j)-psi(i-1,j))))/dy;
                            vp=0.25*((-s3(i)+psi(i,j-1))+abs((s3(i)-psi(i,j-1))))/dx;
                            vm=0.25*((-s3(i)+psi(i,j-1))-abs((s3(i)-psi(i,j-1))))/dx;
                            f(i,j)=(w(i+1,j)+w(i-1,j)-4*w(i,j)+wright(i)+w(i,j-1))/(Re*dx^2)-((psi(i+1,j)-psi(i-1,j))*(wright(i)-w(i,j-1))+(psi(i,j-1)-s3(i))*(w(i+1,j)-w(i-1,j)))/(4*dx^2)+0.5*(up+um)*(wright(i)+w(i,j-1)-2*w(i,j))/dx+0.5*(vp+vm)*(w(i+1,j)+w(i-1,j)-2*w(i,j))/dy;
                        else
                            up=0.25*((psi(i+1,j)-psi(i-1,j))+abs((psi(i+1,j)-psi(i-1,j))))/dy;
                            um=0.25*((psi(i+1,j)-psi(i-1,j))-abs((psi(i+1,j)-psi(i-1,j))))/dy;
                            vp=0.25*((-psi(i,j+1)+psi(i,j-1))+abs((psi(i,j+1)-psi(i,j-1))))/dx;
                            vm=0.25*((-psi(i,j+1)+psi(i,j-1))-abs((psi(i,j+1)-psi(i,j-1))))/dx;
                            f(i,j)=(w(i+1,j)+w(i-1,j)-4*w(i,j)+w(i,j+1)+w(i,j-1))/(Re*dx^2)-((psi(i+1,j)-psi(i-1,j))*(w(i,j+1)-w(i,j-1))+(psi(i,j-1)-psi(i,j+1))*(w(i+1,j)-w(i-1,j)))/(4*dx^2)+0.5*(up+um)*(w(i,j+1)+w(i,j-1)-2*w(i,j))/dx+0.5*(vp+vm)*(w(i+1,j)+w(i-1,j)-2*w(i,j))/dy;
                        end
                    end
                end
            end
        end
    end
    for i=1:n
        for j=1:m
            w1(i,j)=w(i,j)+dt*f(i,j);
        end
    end
    w=w1;
    %% Streamfunction calculation %%
    res=1;
    iter=0;
    while(res>0.0001)
        for i=1:n
            for j=1:m
                if (i==1)
                    if (j==1)
                        psi(i,j)=0.25*(psi0(i+1,j)+psi0(i,j+1)+s2(j)+s1(i)+w(i,j)*dx^2);
                    else if (j==m)
                            psi(i,j)=0.25*(psi0(i+1,j)+psi0(i,j-1)+s2(j)+s3(i)+w(i,j)*dx^2);
                        else
                            psi(i,j)=0.25*(psi0(i+1,j)+psi0(i,j+1)+psi0(i,j-1)+s2(j)+w(i,j)*dx^2);
                        end
                    end
                else if (i==n)
                        if (j==1)
                            psi(i,j)=0.25*(psi0(i-1,j)+s4(j)+psi0(i,j+1)+s1(i)+w(i,j)*dx^2);
                        else if (j==m)
                                psi(i,j)=0.25*(s3(i)+s4(j)+psi0(i,j-1)+psi0(i-1,j)+w(i,j)*dx^2);
                            else
                                psi(i,j)=0.25*(psi0(i-1,j)+s4(j)+psi0(i,j-1)+psi0(i,j+1)+w(i,j)*dx^2);
                            end
                        end
                    else
                        if (j==1)
                            psi(i,j)=0.25*(psi0(i+1,j)+psi0(i,j+1)+psi0(i-1,j)+s1(i)+w(i,j)*dx^2);
                        else if (j==m)
                                psi(i,j)=0.25*(s3(i)+psi0(i+1,j)+psi0(i,j-1)+psi0(i-1,j)+w(i,j)*dx^2);
                            else
                                psi(i,j)=0.25*(psi0(i+1,j)+psi0(i,j+1)+psi0(i,j-1)+psi0(i-1,j)+w(i,j)*dx^2);
                            end
                        end
                    end
                end
            end
        end
        h=0;
        res=zeros(1,m*n);
        for i=1:n
            for j=1:m
                k=j+h;
                if (i==1)
                    if (j==1)
                        res(k)=psi(i,j)-0.25*(psi(i+1,j)+psi(i,j+1)+s2(j)+s1(i)+w(i,j)*dx^2);
                    else if (j==m)
                            res(k)=psi(i,j)-0.25*(psi(i+1,j)+psi(i,j-1)+s2(j)+s3(i)+w(i,j)*dx^2);
                        else
                            res(k)=psi(i,j)-0.25*(psi(i+1,j)+psi(i,j+1)+psi(i,j-1)+s2(j)+w(i,j)*dx^2);
                        end
                    end
                else if (i==n)
                        if (j==1)
                            res(k)=psi(i,j)-0.25*(psi(i-1,j)+s4(j)+psi(i,j+1)+s1(i)+w(i,j)*dx^2);
                        else if (j==m)
                                res(k)=psi(i,j)-0.25*(s3(i)+s4(j)+psi(i,j-1)+psi(i-1,j)+w(i,j)*dx^2);
                            else
                                res(k)=psi(i,j)-0.25*(psi(i-1,j)+s4(j)+psi(i,j-1)+psi(i,j+1)+w(i,j)*dx^2);
                            end
                        end
                    else
                        if (j==1)
                            res(k)=psi(i,j)-0.25*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+s1(i)+w(i,j)*dx^2);
                        else if (j==m)
                                res(k)=psi(i,j)-0.25*(s3(i)+psi(i+1,j)+psi(i,j-1)+psi(i-1,j)+w(i,j)*dx^2);
                            else
                                res(k)=psi(i,j)-0.25*(psi(i+1,j)+psi(i,j+1)+psi(i,j-1)+psi(i-1,j)+w(i,j)*dx^2);
                            end
                        end
                    end
                end
            end
            h=m*i;
        end
        res=rms(res);
        psi0=psi;
        iter=iter+1;
    end
    figure(1)
    contour(psi0,80)
    hold on
    title(sprintf('Streamlines for Re=%d and Time =%d s',Re,time_step*dt+dt))
    hold off
    drawnow
    time_step=time_step+1;
end
%% Vorticity field plot
figure(2)
contour(w1,80)
hold on
title(sprintf('Vorticity for Re=%d and Time =%d s',Re,time_step*dt))
hold off

%% Velocity field plot
u=zeros(n,m);
for i=1:n
    for j=1:m
        if i==1
            u(i,j)=(psi(i+1,j)-s2(j))/dy;
            
        else if i==m
                u(i,j)=(s4(j)-psi(i-1,j))/dy;
            else
                
                u(i,j)=(psi(i+1,j)-psi(i-1,j))/dy;
            end
        end
    end
end
v=zeros(n,m);
for i=1:n
    for j=1:m
        if j==1
            v(i,j)=(-psi(i,j+1)+s1(i))/dx;
            
        else if j==m
                v(i,j)=(-s3(i)+psi(i,j-1))/dx;
            else
                
                v(i,j)=(-psi(i,j+1)+psi(i,j-1))/dx;
            end
        end
    end
end
[x,y] = meshgrid(dx:dx:1-dx,dy:dy:1-dy);
figure(3)
quiver(x,y,u,v);
hold on
title('Velocity vectors')