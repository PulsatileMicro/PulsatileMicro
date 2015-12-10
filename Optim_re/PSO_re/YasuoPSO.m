function [xm,fv] = YasuoPSO(fitness,N,swarms,c1,c2,M,D,lb,ub)
phi = c1 + c2;
if phi <= 4
    disp('c1 与 c2 的 和 必 须 大 于 4 ！');
    xm = NaN;
    fv = NaN;
    return;
end
format long;

%------初始化种群的个体------------
save('AdapCoeff_Wall_res.txt','swarms','-ascii','-v6','-append');
for i=1:N
    x(i,:) = swarms(i,:);
    for j=1:D

%         x(i,j)=randn;  %随机初始化位置

        v(i,j)=randn;  %随机初始化速度

    end

end

%------先计算各个粒子的适应度，并初始化Pi和Pg----------------------

for i=1:N

    p(i)=fitness(x(i,:));

    y(i,:)=x(i,:);

end

pg = x(N,:);             %Pg为全局最优
pgbest = p(N);
for i=1:(N-1)

%     if fitness(x(i,:))<fitness(pg)
% 
%         pg=x(i,:);
% 
%     end
    
    if p(i) < pgbest
        pg=x(i,:);
        pgbest=p(i);
    end
end

%------进入主要循环，按照公式依次迭代------------
vb = (ub-lb).*0.3; %% vmax 设置为问题空间的20%
for t=1:M
    tic;
    for i=1:N
        ksi = 2 / abs(2 - phi - sqrt(phi^2 - 4*phi));
        v(i,:) = v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        v(i,:) = ksi*v(i,:);
        %%disp(['第',num2str(t),'-',num2str(i),'次迭代随机步长为',num2str(v(i,:))]);
       %% test v bounds
        v(i,:)=arrayfun(@v_bounds_mutaion,v(i,:),vb);
        disp(['第',num2str(t),'-',num2str(i),'次迭代约束步长为',num2str(v(i,:))]);
        %%test the lower and uper bounds
        x(i,:)=x(i,:)+v(i,:);
        x(i,:)=arrayfun(@x_bounds_mutaion,x(i,:),lb,ub);
        pbesttmp=fitness(x(i,:));
        if pbesttmp<p(i)
            p(i)=pbesttmp;
            y(i,:)=x(i,:);
        end
%         pgbesttmp=fitness(pg);
        if p(i)<pgbest
            pg=y(i,:);
            pgbest=p(i);
        end

    end
%     Pbest(t)=pgbest;
    elaseTime=toc;
    disp(['YSPSO第',num2str(t),'次迭代运行时间',num2str(elaseTime)]);
    record=[elaseTime,pgbest];
    save('AdapCoeff_Wall_res.txt','record','pg','-ascii','-v6','-append');
end
xm = pg';
fv = pgbest;
end

function y = v_bounds_mutaion(vint, vb)
    if vint<-vb
        y=-vb+2*rand()*vb;

    elseif vint>vb
        y=vb-2*rand()*vb;
    else
        y = vint;
    end
end

function y = x_bounds_mutaion(xint,lb,ub)
    if xint<lb
        y=lb+0.8*rand()*(ub-lb);

    elseif xint>ub
        y=ub-0.8*rand()*(ub-lb);
    else
        y = xint;
    end
end