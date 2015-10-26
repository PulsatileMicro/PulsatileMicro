function [xm,fv] = YasuoPSO(fitness,N,swarms,c1,c2,M,D)
phi = c1 + c2;
if phi <= 4
    disp('c1 �� c2 �� �� �� �� �� �� 4 ��');
    xm = NaN;
    fv = NaN;
    return;
end
format long;

%------��ʼ����Ⱥ�ĸ���------------

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
pgbesttmp = p(N);
for i=1:(N-1)

%     if fitness(x(i,:))<fitness(pg)
% 
%         pg=x(i,:);
% 
%     end
    
    if p(i) < pgbesttmp
        pg=x(i,:);
        pgbesttmp=p(i);
    end
end

%------������Ҫѭ�������չ�ʽ���ε���------------

for t=1:M
    tic;
    for i=1:N
        ksi = 2 / abs(2 - phi - sqrt(phi^2 - 4*phi));
        v(i,:) = v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        v(i,:) = ksi*v(i,:);

        x(i,:)=x(i,:)+v(i,:);
        pbest=fitness(x(i,:));
        if pbest<p(i)

            p(i)=pbest;

            y(i,:)=x(i,:);

        end
        pgbesttmp=fitness(pg);
        if p(i)<pgbesttmp

            pg=y(i,:);

        end

    end
    Pbest(t)=pgbesttmp;
    elaseTime=toc;
    disp(['YSPSO��',num2str(t),'�ε�������ʱ��',num2str(elaseTime)]);
    record=[elaseTime,pgbesttmp];
    save('AdapCoeff_Wall_res.txt','record','pg','-ascii','-v6','-append');
end
xm = pg';
fv = pgbesttmp;



