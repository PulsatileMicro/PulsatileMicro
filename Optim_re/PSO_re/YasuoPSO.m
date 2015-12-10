function [xm,fv] = YasuoPSO(fitness,N,swarms,c1,c2,M,D,lb,ub)
phi = c1 + c2;
if phi <= 4
    disp('c1 �� c2 �� �� �� �� �� �� 4 ��');
    xm = NaN;
    fv = NaN;
    return;
end
format long;

%------��ʼ����Ⱥ�ĸ���------------
save('AdapCoeff_Wall_res.txt','swarms','-ascii','-v6','-append');
for i=1:N
    x(i,:) = swarms(i,:);
    for j=1:D

%         x(i,j)=randn;  %�����ʼ��λ��

        v(i,j)=randn;  %�����ʼ���ٶ�

    end

end

%------�ȼ���������ӵ���Ӧ�ȣ�����ʼ��Pi��Pg----------------------

for i=1:N

    p(i)=fitness(x(i,:));

    y(i,:)=x(i,:);

end

pg = x(N,:);             %PgΪȫ������
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

%------������Ҫѭ�������չ�ʽ���ε���------------
vb = (ub-lb).*0.3; %% vmax ����Ϊ����ռ��20%
for t=1:M
    tic;
    for i=1:N
        ksi = 2 / abs(2 - phi - sqrt(phi^2 - 4*phi));
        v(i,:) = v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        v(i,:) = ksi*v(i,:);
        %%disp(['��',num2str(t),'-',num2str(i),'�ε����������Ϊ',num2str(v(i,:))]);
       %% test v bounds
        v(i,:)=arrayfun(@v_bounds_mutaion,v(i,:),vb);
        disp(['��',num2str(t),'-',num2str(i),'�ε���Լ������Ϊ',num2str(v(i,:))]);
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
    disp(['YSPSO��',num2str(t),'�ε�������ʱ��',num2str(elaseTime)]);
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