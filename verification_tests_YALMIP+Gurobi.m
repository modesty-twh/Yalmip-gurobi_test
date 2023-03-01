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

% read case name
filepath=pwd;           %保存当前工作目录
cd(['E:\博士科研\03 科研成果\01 论文管理\Ppaer-9\revise 1\code\[exp5] cplex\case_FJSP'])           %把当前工作目录切换到指定文件夹

dz = cd;%得到你文件所在的名字
files = cellstr(ls([dz '\capacity*'])); %得到文件下的mat数据的名字
nfiles = size(files,1);

cd(filepath)            %切回原工作目录

nCase=1;
nameCase=files{nCase};

filepath=pwd;           %保存当前工作目录
cd(['E:\博士科研\03 科研成果\01 论文管理\Ppaer-9\revise 1\code\[exp5] cplex\case_FJSP'])           %把当前工作目录切换到指定文件夹

load([ 'info_factory',nameCase(9:end)]);
load([ 'info_job',nameCase(9:end)]);
load([ 'info_vehicle',nameCase(9:end)]);
info_job=info_job([1:5]);

load(nameCase);

% load([ 'info_factory_1']);
% load([ 'info_job_3']);
% load([ 'info_vehicle',nameCase(9:end)]);

cd(filepath)            %切回原工作目录

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

% 中间变量约束1
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

% 同机器时间关系
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

% y自身约束1：每个工序一定要加工一次且仅一次
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

% y自身约束2：同工件的工序必须同工厂加工
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

% y自身约束3：同一机器上的前一位置不小于后一位置
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

% y自身约束4：同一机器上的同一位置不能有多个工序
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

% 同工件时间关系
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

% 开始时间大于等于0
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
% ops = sdpsettings('solver','cplex','verbose',1);
ops = sdpsettings('solver', 'Gurobi+', 'verbose', 1, 'debug', 1, 'gurobi.TimeLimit',7200);
% ops = sdpsettings('solver', 'Gurobi+', 'verbose', 1, 'debug', 1);
ops.gurobi.NonConvex = 2;
% ops.cplex.display='on';
% ops.cplex.timelimit=60;
% ops.cplex.mip.tolerances.mipgap=0.1;

% 诊断求解可行性
disp('开始求解')
diagnostics=optimize(constraint,obj,ops);
if diagnostics.problem==0
    disp('Feasible!')
elseif diagnostics.problem == 1
    disp('Infeasible')
    pause();
% else
%     disp('Timeout, Display the current optimal solution')
end

%% 读结果
recy=[];
for i=1:nj
    no1=info_job(i).nop;
    for j=1:no1
        for f=1:nf
            nm=info_factory(f).nm;
            for m=1:nm
                if sum(info_job(i).fmte{j}{f}(:,1)==m,'all')
                    [s,~]=find(value(y{i}{j}{f}{m}(:))==1);
                    if ~isempty(s)
                        recy=[recy;[i,j,f,m,s]];
                    end
                end
            end
        end
    end
end

recC=[];
for i=1:nj
    recC=[recC;{value(c{i})}];
end

%% plot gantt
PlotGantt(recC,recy,info_job,info_factory,value(obj))
