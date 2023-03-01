% �޸�Ŀ�꺯��
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
filepath=pwd;           %���浱ǰ����Ŀ¼
cd(['E:\��ʿ����\03 ���гɹ�\01 ���Ĺ���\Ppaer-9\revise 1\code\[exp5] cplex\case_FJSP'])           %�ѵ�ǰ����Ŀ¼�л���ָ���ļ���

dz = cd;%�õ����ļ����ڵ�����
files = cellstr(ls([dz '\capacity*'])); %�õ��ļ��µ�mat���ݵ�����
nfiles = size(files,1);

cd(filepath)            %�л�ԭ����Ŀ¼

nCase=1;
nameCase=files{nCase};

filepath=pwd;           %���浱ǰ����Ŀ¼
cd(['E:\��ʿ����\03 ���гɹ�\01 ���Ĺ���\Ppaer-9\revise 1\code\[exp5] cplex\case_FJSP'])           %�ѵ�ǰ����Ŀ¼�л���ָ���ļ���

load([ 'info_factory',nameCase(9:end)]);
load([ 'info_job',nameCase(9:end)]);
load([ 'info_vehicle',nameCase(9:end)]);
info_job=info_job([1:5]);

load(nameCase);

% load([ 'info_factory_1']);
% load([ 'info_job_3']);
% load([ 'info_vehicle',nameCase(9:end)]);

cd(filepath)            %�л�ԭ����Ŀ¼

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


%% �������

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


%% Ŀ�꺯��
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

%% Լ������
constraint=[];

% �м����Լ��1
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

% ͬ����ʱ���ϵ
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

% y����Լ��1��ÿ������һ��Ҫ�ӹ�һ���ҽ�һ��
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

% y����Լ��2��ͬ�����Ĺ������ͬ�����ӹ�
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

% y����Լ��3��ͬһ�����ϵ�ǰһλ�ò�С�ں�һλ��
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

% y����Լ��4��ͬһ�����ϵ�ͬһλ�ò����ж������
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

% ͬ����ʱ���ϵ
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

% ��ʼʱ����ڵ���0
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

%% ���
% ops = sdpsettings('solver','cplex','verbose',1);
ops = sdpsettings('solver', 'Gurobi+', 'verbose', 1, 'debug', 1, 'gurobi.TimeLimit',7200);
% ops = sdpsettings('solver', 'Gurobi+', 'verbose', 1, 'debug', 1);
ops.gurobi.NonConvex = 2;
% ops.cplex.display='on';
% ops.cplex.timelimit=60;
% ops.cplex.mip.tolerances.mipgap=0.1;

% �����������
disp('��ʼ���')
diagnostics=optimize(constraint,obj,ops);
if diagnostics.problem==0
    disp('Feasible!')
elseif diagnostics.problem == 1
    disp('Infeasible')
    pause();
% else
%     disp('Timeout, Display the current optimal solution')
end

%% �����
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
