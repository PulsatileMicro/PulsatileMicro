function [xm,fv] = QuaPSO(fitness,N,swarms,M,D,lb,ub)
%QPSO Summary of this function goes here
%   Detailed explanation goes here
alpha = 0.5; %����ϵ��beita����
format long;

%------��ʼ����Ⱥ�ĸ���------------
save('AdapCoeff_Wall_res.txt','swarms','-ascii','-v6','-append');
for i=1:N
    x(i,:) = swarms(i,:);
    %     for j=1:D
    %
    % %         x(i,j)=randn;  %�����ʼ��λ��
    %
    % %         v(i,j)=randn;  %�����ʼ���ٶ�
    %
    %     end
    
end

%------�ȼ���������ӵ���Ӧ�ȣ�����ʼ��Pi��Pg----------------------

for i=1:N
    
    fp(i)=fitness(x(i,:));
    
    p(i,:)=x(i,:);
    
end

pg = x(N,:);             %PgΪȫ������
fpg = fp(N);
for i=1:N-1
    if fp(i) < fpg
        pg=x(i,:);
        fpg=fp(i);
    end
end

%------������Ҫѭ�������չ�ʽ���ε���------------
for t=1:M
    tic;
    mbest=mean(p);
    beita=(1-alpha)*(M-t+1)/M+alpha;
    for i=1:N
        %   x(i,:)=arrayfun(@x_bounds_mutaion,x(i,:),lb,ub);
        for d=1:D
            rng('shuffle');
            fy=rand;
            p(i,d)=fy*p(i,d)+(1-fy)*pg(d);
%             rng('shuffle');
            fy=rand;
            if rand < 0.5
                x(i,d) = p(i,d)+beita*abs(mbest(d) - x(i,d))*log(1/fy);
            else
                x(i,d) = p(i,d)-beita*abs(mbest(d) - x(i,d))*log(1/fy);
            end
%             if x(i,d)<lb(d)
%                 x(i,d)=lb(d)+0.8*rand*(ub(d)-lb(d));
%             elseif x(i,d)>ub(d)
%                 x(i,d)=ub(d)-0.8*rand*(ub(d)-lb(d));
%             end
        end
        fpgtemp = fitness(x(i,:));
        if fpgtemp < fp(i)
            p(i,:) = x(i,:);
            fp(i) = fpgtemp;
            if fpgtemp < fpg
                pg = x(i,:);
                fpg = fpgtemp;
            end
        end
        
        %         disp(['��',num2str(t),'-',num2str(i),'�ε���Լ������Ϊ',num2str(v(i,:))]);
        
    end
    %   Pbest(t)=fpg;
    elaseTime=toc;
    disp(['YSPSO��',num2str(t),'�ε�������ʱ��',num2str(elaseTime)]);
    record=[elaseTime,fpg];
    save('AdapCoeff_Wall_res.txt','record','pg','-ascii','-v6','-append');
end
xm = pg';
fv = fpg;
end

