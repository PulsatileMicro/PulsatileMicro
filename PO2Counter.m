%氧分压计算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PO2,SO2in,SO2mid,SO2out]=PO2Counter(Porder,BSO2in,FromNeg,ToNeg,Hd,Len,MeanFlow,Eju)
%氧体积数=血流量*0.5*红细胞比容*氧体积分数
%%PO2计算
%% 生理参数单位转换 %%%%
MeanFlow=MeanFlow./60./1e6;   %血流单位转换cm3/s
MPO2=4*1e-11;  %平均氧耗
MPO=MPO2*1e4/60*Len/10;  % 氧耗单位转换O2cm3.cm3/s

%% 主要参数初值矩阵 %%%%
VelNum=length(Porder);   %血管总数
SO2in=zeros(VelNum,1);   %入氧体积分数
SO2mid=1e-10*ones(VelNum,1); %血管中点氧体积分数
SO2out=zeros(VelNum,1);  %出氧体积分数
JconMin=zeros(VelNum,1);  %入氧体积数
JconMmid=zeros(VelNum,1);  %血管中点氧体积数
JconMout=zeros(VelNum,1);  %出氧体积数
PO2=zeros(VelNum,1);  %氧分压

%% 入边界各参数赋值 %%%%
SO2index=find(BSO2in>0);
for i=1:length(SO2index)
  SO2in(SO2index(i))=BSO2in(SO2index(i));
  JconMin(SO2index(i))=0.5*MeanFlow(SO2index(i)).*Hd(SO2index(i)).*SO2in(SO2index(i));
  JconMmid(SO2index(i))=JconMin(SO2index(i))-MPO((SO2index(i)))/2;
  JconMout(SO2index(i))=JconMin(SO2index(i))-MPO((SO2index(i)));
  SO2mid(SO2index(i))=JconMmid(SO2index(i))/(0.5*MeanFlow(SO2index(i)).*Hd(SO2index(i)));
  SO2out(SO2index(i))=JconMout(SO2index(i))/(0.5*MeanFlow(SO2index(i)).*Hd(SO2index(i)));
  PO2(SO2index(i))=Hill(SO2mid(SO2index(i)));
  if SO2mid(SO2index(i))<0 %氧体积分数默认不存在负值
    SO2mid(SO2index(i))=0;
    SO2out(SO2index(i))=0;
    PO2(SO2index(i))=0;
  else
    %氧分压计算公式：希尔公式
    PO2(SO2index(i))=Hill(SO2mid(SO2index(i)));
  end
end

if Eju==0  %正常情况
  for i=1:length(Porder)
    j=Porder(i); %导入计算顺序矩阵
    
    ConvergeIndex=find(ToNeg==FromNeg(j));  %判断汇聚血管
    if length(ConvergeIndex)==2   %汇聚
      
      %氧体积数守恒
      JconMin(j)=MeanFlow(ConvergeIndex(1))*0.5*Hd(ConvergeIndex(1))*SO2out(ConvergeIndex(1))+MeanFlow(ConvergeIndex(2))*0.5*Hd(ConvergeIndex(2))*SO2out(ConvergeIndex(2));
      SO2in(j)=JconMin(j)/(MeanFlow(j)*0.5*Hd(j));
      JconMmid(j)=JconMin(j)-MPO(j)/2;
      JconMout(j)=JconMin(j)-MPO(j);
      SO2mid(j)=JconMmid(j)/(MeanFlow(j)*0.5*Hd(j));
      SO2out(j)=JconMout(j)/(MeanFlow(j)*0.5*Hd(j));
      if Hd(j)==0
        SO2in(j)=0;
        SO2out(j)=0;
        PO2(j)=0;
      elseif SO2mid(j)<0
        SO2mid(j)=0;
        SO2out(j)=0;
        PO2(j)=0;
      elseif SO2mid(j)>1
        PO2(j)=95;
      else
        PO2(j)=Hill(SO2mid(j));
      end
    end
    
    BifurIndex=find(FromNeg==ToNeg(j));  %判断分叉血管
    if length(BifurIndex)==2   %分叉
      
      %氧体积分数传入下级子血管（氧体积数守恒）
      SO2in(BifurIndex(1))=SO2out(j);
      SO2in(BifurIndex(2))=SO2out(j);
      JconMin(BifurIndex(1))=MeanFlow(BifurIndex(1))*0.5*Hd(BifurIndex(1))*SO2in(BifurIndex(1));
      JconMin(BifurIndex(2))=MeanFlow(BifurIndex(2))*0.5*Hd(BifurIndex(2))*SO2in(BifurIndex(2));
      JconMmid(BifurIndex(1))=JconMin(BifurIndex(1))-MPO(BifurIndex(1))/2;
      JconMmid(BifurIndex(2))=JconMin(BifurIndex(2))-MPO(BifurIndex(2))/2;
      JconMout(BifurIndex(1))=JconMin(BifurIndex(1))-MPO(BifurIndex(1));
      JconMout(BifurIndex(2))=JconMin(BifurIndex(2))-MPO(BifurIndex(2));
      SO2mid(BifurIndex(1))=JconMmid(BifurIndex(1))/(MeanFlow(BifurIndex(1))*0.5*Hd(BifurIndex(1)));
      SO2mid(BifurIndex(2))=JconMmid(BifurIndex(2))/(MeanFlow(BifurIndex(2))*0.5*Hd(BifurIndex(2)));
      SO2out(BifurIndex(1))=JconMout(BifurIndex(1))/(MeanFlow(BifurIndex(1))*0.5*Hd(BifurIndex(1)));
      SO2out(BifurIndex(2))=JconMout(BifurIndex(2))/(MeanFlow(BifurIndex(2))*0.5*Hd(BifurIndex(2)));
      if Hd(BifurIndex(1))==0
        SO2mid(BifurIndex(1))=0;
        SO2out(BifurIndex(1))=0;
      end
      if Hd(BifurIndex(2))==0
        SO2mid(BifurIndex(2))=0;
        SO2out(BifurIndex(2))=0;
      end
      PO2(BifurIndex(1))=Hill(SO2mid(BifurIndex(1)));
      PO2(BifurIndex(2))=Hill(SO2mid(BifurIndex(2)));
      if SO2mid(BifurIndex(1))<0
        SO2out(BifurIndex(1))=0;
        PO2(BifurIndex(1))=0;
      end
      if SO2mid(BifurIndex(2))<0
        SO2out(BifurIndex(2))=0;
        PO2(BifurIndex(2))=0;
      end
      
    elseif (length(find(ToNeg==ToNeg(j)))==1 && length(BifurIndex)==1)  %单连血管
      
      %氧体积分数传入（氧体积数守恒）
      SO2in(BifurIndex(1))=SO2out(j);
      JconMin(BifurIndex(1))=MeanFlow(BifurIndex(1))*0.5*Hd(BifurIndex(1))*SO2in(BifurIndex(1));
      JconMmid(BifurIndex(1))=JconMin(BifurIndex(1))-MPO(BifurIndex(1))/2;
      JconMout(BifurIndex(1))=JconMin(BifurIndex(1))-MPO(BifurIndex(1));
      SO2mid(BifurIndex(1))=JconMmid(BifurIndex(1))/(MeanFlow(BifurIndex(1))*0.5*Hd(BifurIndex(1)));
      SO2out(BifurIndex(1))=JconMout(BifurIndex(1))/(MeanFlow(BifurIndex(1))*0.5*Hd(BifurIndex(1)));
      PO2(BifurIndex(1))=Hill(SO2mid(BifurIndex(1)));
      if SO2mid(j)<0
        SO2mid(j)=0;
        SO2out(j)=0;
        PO2(j)=0;
      elseif SO2mid(j)>1
        PO2(j)=95;
      else
        PO2(j)=Hill(SO2mid(j));
      end
      
    else %分叉血管的子血管（已在分叉血管中被计算）
      continue;
    end
  end
else
  PO2=zeros(VelNum,1);
end

end