% 修改目标函数
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

toc

%% 求解
ops = sdpsettings('solver', 'Gurobi+', 'verbose', 1, 'debug', 1, 'gurobi.TimeLimit',7200);
ops.gurobi.NonConvex = 2;