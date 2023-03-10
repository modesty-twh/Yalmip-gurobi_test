clear
close all
clc
dbstop if error
tic

%% paratmeter
ncus='100';
% model parameters
Re=0.785; % conversion coefficient of electricity
Rg=2.62; % conversion coefficient of gasoline
lam_d=0.05; % alpha1
lam_v=1.6; % alpha2
p_idle=2.5; % idle power
L=500; % very large number
epsilon=0.001; % very small number

% load case
load('info_factory');
load('info_job');
load('info_vehicle');
load('capacity');

vWeight=capacity; % weight of vehicle = its capacity

nj=numel(info_job);
nf=numel(info_factory);
nv=numel(info_vehicle);
nop=0;
for i4=1:nj
    no2=info_job(i4).nop;
    for i5=1:no2
        nop=nop+1;
    end
end


%% 定义变量
for i=1:nj
    no1=info_job(i).nop;
    for j=1:no1
        for f=1:nf
            nm=info_factory(f).nm;
            for m=1:nm
                if sum(info_job(i).fmte{j}{f}(:,1)==m,'all')
                    y{i}{j}{f}{m}=binvar(nop,1);
                else
                    y{i}{j}{f}{m}=zeros(nop,1);
                end
            end
        end
    end
end

for i=1:nj
    no=info_job(i).nop;
    c{i}=sdpvar(no,1,'full');
end

cf=sdpvar(1,1,'full');

a=sdpvar(nj,1,'full');
Q=intvar(nj,1,'full');
g=sdpvar(nj+nf,nj+nf,'full');
for i=1:nj+nf
    g(i,i,:)=0;
end

for i=1:nj
    no=info_job(i).nop;
    c{i}=sdpvar(no,1,'full');
end

si=sdpvar(nj,1,'full');

%% 目标函数
% f1
tem=0;
for i=1:nj
    tem=max(c{i}(end),tem);
end
obj1=tem;

% f1
E1=0;
E2=0;
for f=1:nf
    nm=info_factory(f).nm;
    for m=1:nm
        T=0;
        for i=1:nj
            no=info_job(i).nop;
            for j=1:no
                pos=info_job(i).fmte{j}{f}(:,1)==m;
                if sum(pos)>0
                    E1=E1+info_job(i).fmte{j}{f}(pos,2)*info_job(i).fmte{j}{f}(pos,3)*sum(y{i}{j}{f}{m});
                    T=T+info_job(i).fmte{j}{f}(pos,2)*sum(y{i}{j}{f}{m});
                end
            end
        end
        E2=E2+p_idle*(cf-T);
    end    
end
obj2=Re*(E1+E2);

obj=obj1+obj2/200;

%% 约束条件
constraint=[];

tem=0;
for f=1:nf    
    for m=1:nm
        for i=1:nj
            no1=info_job(i).nop;
            for j=1:no1
                tem=max(tem,c{i}(j)*sum(y{i}{j}{f}{m}));
            end
        end
    end    
end
constraint=[constraint,cf>=tem];

for f=1:nf
    nm=info_factory(f).nm;
    for m=1:nm
        for i=1:nj
            no1=info_job(i).nop;
            for j=1:no1
                for i1=1:nj
                    no2=info_job(i1).nop;
                    for j1=1:no2
                        if sum(info_job(i).fmte{j}{f}(:,1)==m,'all')&&sum(info_job(i1).fmte{j1}{f}(:,1)==m,'all')
                            pos=info_job(i1).fmte{j1}{f}(:,1)==m;
                            for p=1:nop-1
                                constraint=[constraint,...
                                    c{i1}(j1)-info_job(i1).fmte{j1}{f}(pos,2)*y{i}{j}{f}{m}(p)*y{i1}{j1}{f}{m}(p+1)+...
                                    (1-y{i}{j}{f}{m}(p)*y{i1}{j1}{f}{m}(p+1))*L>=c{i}(j)];
                            end
                        end
                    end
                end
            end
        end
    end
end

for i=1:nj
    no=info_job(i).nop;
    for j=1:no
        tem=0;
        for f=1:nf
            nm=info_factory(f).nm;
            for m=1:nm
                tem=tem+sum(y{i}{j}{f}{m}(:),'all');
            end
        end
        constraint=[constraint,tem==1];
    end
end

for i=1:nj
    for f=1:nf
        no=info_job(i).nop;
        for j=1:no
            for jj=j+1:no
                tem1=0;
                tem2=0;
                nm=info_factory(f).nm;
                for m=1:nm
                    tem1=tem1+sum(y{i}{j}{f}{m}(:),'all');
                    tem2=tem2+sum(y{i}{jj}{f}{m}(:),'all');
                end
                constraint=[constraint,tem1==tem2];
            end
        end
    end
end

for f=1:nf
    nm=info_factory(f).nm;
    for m=1:nm
        for p=1:nop-1
            tem1=0;
            for i=1:nj
                no=info_job(i).nop;
                for j=1:no
                    tem1=tem1+y{i}{j}{f}{m}(p);
                end
            end
            tem2=0;
            for i=1:nj
                no=info_job(i).nop;
                for j=1:no
                    tem2=tem2+y{i}{j}{f}{m}(p+1);
                end
            end
            constraint=[constraint,tem1>=tem2];
        end
    end
end

for f=1:nf
    nm=info_factory(f).nm;
    for m=1:nm
        for p=1:nop
            tem=0;
            for i=1:nj
                no=info_job(i).nop;
                for j=1:no
                    tem=tem+y{i}{j}{f}{m}(p);
                end
            end
            constraint=[constraint,tem<=1];
        end
    end
end

for f=1:nf
    nm=info_factory(f).nm;
    for m=1:nm
        for i=1:nj
            no=info_job(i).nop;
            for j=2:no
                pos=info_job(i).fmte{j}{f}(:,1)==m;
                if sum(pos,'all')>0
                    constraint=[constraint,...
                        c{i}(j)-info_job(i).fmte{j}{f}(pos,2)*sum(y{i}{j}{f}{m}(:),'all')+(1-sum(y{i}{j}{f}{m}(:),'all'))*L>=c{i}(j-1)];
                end
            end
        end
    end
end

for f=1:nf
    nm=info_factory(f).nm;
    for i=1:nj
        no=info_job(i).nop;
        for j=1:no            
            for m=1:nm
                pos=info_job(i).fmte{j}{f}(:,1)==m;
                if sum(pos,'all')>0
                    constraint=[constraint,...
                        c{i}(j)-info_job(i).fmte{j}{f}(pos,2)*sum(y{i}{j}{f}{m}(:),'all')+(1-sum(y{i}{j}{f}{m}(:),'all'))*L>=0];
                end
            end
        end
    end
end

for i=1:nj
    for v=1:nv
        constraint=[constraint,sum(x(i,:,v))==sum(z(i,:,v))];
    end
end

for i=1:nj
    constraint=[constraint,sum(x(i,:,:),'all')==1];
end

for i=1:nj
    constraint=[constraint,si(i)>=c{i}(end)];
end

for f=1:nf
    for i1=1:nj
        for i2=1:nj
            for v=1:nv
                constraint=[constraint,si(i1)*z(i1,i2,v)*x(i1,f,v)*x(i2,f,v)==si(i2)*z(i1,i2,v)*x(i1,f,v)*x(i2,f,v)];
            end
        end
    end
end

for f=1:nf
    for i1=1:nj
        for i2=1:nj
            for v=1:nv
                constraint=[constraint,Q(i2)+info_job(i2).demand*z(i1,i2,v)*x(i1,f,v)*(1-x(i2,f,v))<=capacity];
            end
        end
    end
end

for f=1:nf
    for i1=1:nj
        for i2=1:nj
            for v=1:nv
                constraint=[constraint,Q(i1)*z(i1,i2,v)*x(i1,f,v)*(1-x(i2,f,v))==0];
            end
        end
    end
end

for f=1:nf
    for i1=1:nj
        for i2=1:nj
            for v=1:nv
                constraint=[constraint,(Q(i1)-info_job(i2).demand-Q(i2))*z(i1,i2,v)*x(i1,f,v)*x(i2,f,v)==0];
            end
        end
    end
end

for f=1:nf
    for i=1:nj
        for v=1:nv
            constraint=[constraint,a(i)-norm(info_factory(f).xy-info_job(i).xy) *z(end,i,v)*x(i,f,v)+...
                (1-z(end,i,v)*x(i,f,v))*L>=max(si(i),norm(info_factory(info_vehicle(v)).xy-info_factory(f).xy))];
        end
    end
end

for f=1:nf
    for i1=1:nj
        for i2=1:nj
            for v=1:nv
                if i1~=i2
                constraint=[constraint,...
                    a(i2)-norm(info_job(i1).xy-info_job(i2).xy)*z(i1,i2,v)*x(i1,f,v)*x(i2,f,v)+...
                    (1-z(i1,i2,v)*x(i1,f,v)*x(i2,f,v))*L>=a(i1)];
                end
            end
        end
    end
end

for f=1:nf
    for i1=1:nj
        for i2=1:nj
            for v=1:nv
                if i1~=i2
                constraint=[constraint,...
                    a(i2)-(norm(info_job(i1).xy-info_factory(f).xy)+norm(info_factory(f).xy-info_job(i2).xy))*...
                    z(i1,i2,v)*x(i2,f,v)*(1-x(i1,f,v))+...
                    (1-z(i1,i2,v)*x(i2,f,v)*(1-x(i1,f,v)))*L>=a(i1)];
                end
            end
        end
    end
end

for i=1:nj
    for f=1:nf
        constraint=[constraint,a(i)-norm(info_factory(f).xy-info_job(i).xy)*sum(x(i,f,:))>=si(i)];
    end    
end

for i=1:nj
    tem=0;
    for i1=1:nj+1
        for v=1:nv
            if i~=i1
                tem=tem+z(i,i1,v);
            end
        end
    end
    constraint=[constraint,tem==1];
end

for i=1:nj
    tem=0;
    for i1=1:nj+1
        for v=1:nv
            if i~=i1
                tem=tem+z(i1,i,v);
            end
        end
    end
    constraint=[constraint,tem==1];
end

for v=1:nv
    constraint=[constraint,sum(z(end,:,v))==1];
end

for v=1:nv
    constraint=[constraint,sum(z(:,end,v))==1];
end

for i=1:nj
    for v=1:nv
        constraint=[constraint,sum(z(:,i,v))==sum(z(i,:,v))];
    end
end
toc

%% 求解
ops = sdpsettings('solver', 'Gurobi+', 'verbose', 1, 'debug', 1, 'gurobi.TimeLimit',7200);
ops.gurobi.NonConvex = 2;