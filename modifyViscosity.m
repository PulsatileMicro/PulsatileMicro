%ճ�Ͷ�����
function y=modifyViscosity(Porder,FromNeg,ToNeg,DebugVisc,T,modifyType,Alpha,Eju)

switch num2str(modifyType)
  case '0'   %������
    y=DebugVisc;  %Pa.s
  case '1'   %��bifurcation��SOR����
    DebugVisc=modifyBifurcation(Porder,FromNeg,ToNeg,DebugVisc,T,Alpha,Eju);
    y=DebugVisc;  %Pa.s
  case '2'   %������Ѫ����SOR����
    DebugVisc=modifyAll(DebugVisc,T,Alpha,Eju);
    y=DebugVisc;  %Pa.s
  otherwise
    error('Wrong setting of "modifyType" !');
end

%�Ӻ�������bifurcationӦ��SOR��������
  function y=modifyBifurcation(Porder,FromNeg,ToNeg,DebugVisc,T,Alpha,Eju)
    count=1;  %�ֲ������ֵ
    if Eju==0  %֮ǰ������������
      for i=1:length(Porder)
        j=Porder(i);  %��Ѫ��˳���ҳ�diver bifur
        bifurIndex=find(FromNeg==ToNeg(j));
        if length(bifurIndex)==2
          BifurMatrix(count,1:3)=[j bifurIndex(1) bifurIndex(2)];
          count=count+1;
        end
      end
      
      for i=1:length(BifurMatrix(:,1))
        if T>2
          Visc2=(1-Alpha)*DebugVisc(BifurMatrix(i,2),T)+Alpha*DebugVisc(BifurMatrix(i,2),T-1);
          Visc3=(1-Alpha)*DebugVisc(BifurMatrix(i,3),T)+Alpha*DebugVisc(BifurMatrix(i,3),T-1);
          DebugVisc(BifurMatrix(i,2),T)=Visc2;
          DebugVisc(BifurMatrix(i,3),T)=Visc3;
        end
      end
      y=DebugVisc;
    else
      y=DebugVisc;
    end
  end

%�Ӻ�������All vessel segmentsӦ��SOR��������
  function y=modifyAll(DebugVisc,T,Alpha,Eju)
    if Eju==0
      DebugVisc(:,T)=(1-Alpha)*DebugVisc(:,T)+Alpha*DebugVisc(:,T-1);
      y=DebugVisc(:,T);  %Pa.s
    else
      y=DebugVisc(:,T);  %Pa.s
    end
  end

end