%% 读取不同管径、不同输入频率、不同粘滞度、不同杨氏模量情况下的衰减和传播速度情况，并绘图
% 导入Origin绘图
clear;clc;close all;
[num1,txt1,raw1]=xlsread('Attenuation in single vessel.xlsx','小血管频域分析');

%% 不同管径
All_PIP_D=zeros(6,6);   % 行：不同管径，不同粘滞度，不同杨氏模量；列：不同频率
All_PIU_D=zeros(6,6);
All_PWV_D=zeros(6,6);
for i=1:6
  for j=1:6
    All_PIP_D(j,i)=num1(1+9*(i-1),4+5*(j-1));
    All_PIU_D(j,i)=num1(2+9*(i-1),4+5*(j-1));
    All_PWV_D(j,i)=num1(5+9*(i-1),4+5*(j-1));
  end
end

%% 不同粘滞度
All_PIP_V=zeros(6,5);   % 行：不同管径，不同粘滞度，不同杨氏模量；列：不同频率
All_PIU_V=zeros(6,5);
All_PWV_V=zeros(6,5);
for i=1:5
  for j=1:6
    All_PIP_V(j,i)=num1(55+9*(i-1),4+5*(j-1));
    All_PIU_V(j,i)=num1(56+9*(i-1),4+5*(j-1));
    All_PWV_V(j,i)=num1(59+9*(i-1),4+5*(j-1));
  end
end

%% 不同杨氏模量
All_PIP_E=zeros(6,8);   % 行：不同管径，不同粘滞度，不同杨氏模量；列：不同频率
All_PIU_E=zeros(6,8);
All_PWV_E=zeros(6,8);
for i=1:8
  for j=1:6
    All_PIP_E(j,i)=num1(100+9*(i-1),4+5*(j-1));
    All_PIU_E(j,i)=num1(101+9*(i-1),4+5*(j-1));
%     All_PWV_E(j,i)=num1(104+9*(i-1),4+5*(j-1));
  end
end

[num2,txt2,raw2]=xlsread('Attenuation in single vessel.xlsx','弹性vs.粘弹性');
%% 不同血管壁模型，不同管径
All_PIP_ViscElastic=zeros(6,6);   % 行：不同频率  列：不同管径，不同粘滞度，不同杨氏模量
All_PIU_ViscElastic=zeros(6,6);
All_PWV_ViscElastic=zeros(6,6);
for i=1:6
  for j=1:6
    All_PIP_ViscElastic(j,i)=num2(1+9*(i-1),4+5*(j-1));
    All_PIU_ViscElastic(j,i)=num2(2+9*(i-1),4+5*(j-1));
  end
end
All_PIP_ViscElastic=All_PIP_ViscElastic*100;
All_PIU_ViscElastic=All_PIU_ViscElastic*100;
