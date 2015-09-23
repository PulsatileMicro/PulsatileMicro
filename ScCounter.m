% �����ź� ����->���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sc,Jcmid]=ScCounter(Norder,BJc,FromNeg,ToNeg,Len,Sm,Jo,Lref,Eju)
%% ������λת�� %%%%
umLen=Len*1e3;  %��λת��um

%% ��Ҫ������ֵ����
VesNum=length(Norder);
Jc=1e-10*ones(VesNum,1);
Jcmid=zeros(VesNum,1);

%% ���߽��������ֵ %%%%
Jcindex=find(BJc>0);
for i=1:length(Jcindex)
  Jc(Jcindex(i))=(BJc(Jcindex(i))+umLen(Jcindex(i))*Sm(Jcindex(i)))*exp(-umLen(Jcindex(i))/Lref(i));
  Jcmid(Jcindex(i))=BJc(Jcindex(i))*exp(-umLen(Jcindex(i))/Lref(i))+0.5*umLen(Jcindex(i))*Sm(Jcindex(i))*exp(-umLen(Jcindex(i))/Lref(i));
end

if Eju==0  %�������
  for i=1:length(Norder)
    j=Norder(i); %����Ѫ��˳���������
    %�жϻ��Ѫ�ܣ����������ģ�
    ConvergeIndex=find(FromNeg==ToNeg(j)); 
    if length(ConvergeIndex)==2
      Jc(j)=(Jc(ConvergeIndex(1))+Jc(ConvergeIndex(2))+umLen(j)*Sm(j))*exp(-umLen(j)/Lref(j));
      Jcmid(j)=(Jc(ConvergeIndex(1))+Jc(ConvergeIndex(2)))*exp(-0.5*umLen(j)/Lref(j))+0.5*umLen(j)*Sm(j)*exp(-0.5*umLen(j)/Lref(j));
    end
    
    BifurIndex=find(ToNeg==FromNeg(j));  %�жϷֲ�Ѫ�ܣ����������ģ�
    if length(BifurIndex)==2
      Jc(BifurIndex(1))=(Jc(j)/2+umLen(BifurIndex(1))*Sm(BifurIndex(1)))*exp(-umLen(BifurIndex(1))/Lref(BifurIndex(1)));
      Jcmid(BifurIndex(1))=(Jc(j)/2)*exp(-0.5*umLen(BifurIndex(1))/Lref(BifurIndex(1)))+0.5*umLen(BifurIndex(1))*Sm(BifurIndex(1))*exp(-0.5*umLen(BifurIndex(1))/Lref(BifurIndex(1)));
      Jc(BifurIndex(2))=(Jc(j)/2+umLen(BifurIndex(2))*Sm(BifurIndex(2)))*exp(-umLen(BifurIndex(2))/Lref(BifurIndex(2)));
      Jcmid(BifurIndex(2))=(Jc(j)/2)*exp(-0.5*umLen(BifurIndex(2))/Lref(BifurIndex(2)))+0.5*umLen(BifurIndex(2))*Sm(BifurIndex(2))*exp(-0.5*umLen(BifurIndex(2))/Lref(BifurIndex(2)));
    elseif (length(find(FromNeg==FromNeg(j)))==1 && length(BifurIndex)==1)  %����Ѫ��
      Jc(BifurIndex(1))=(Jc(j)+umLen(BifurIndex(1))*Sm(BifurIndex(1)))*exp(-umLen(BifurIndex(1))/Lref(BifurIndex(1)));
      Jcmid(BifurIndex(1))=(Jc(j))*exp(-0.5*umLen(BifurIndex(1))/Lref(BifurIndex(1)))+0.5*umLen(BifurIndex(1))*Sm(BifurIndex(1))*exp(-0.5*umLen(BifurIndex(1))/Lref(BifurIndex(1)));
    else
      continue;
    end
  end
else %�쳣���
  Jcmid=zeros(VesNum,1);
end

Sc=Jcmid./(Jcmid+Jo);

end