%% Hemodynamic analysis
function [MeanP1 MeanP2 MeanQ MeanU MeanA MeanVisc tAll PAll1 PAll2 UAll QAll AAll ViscAll]...
  =HemoAnal4Network(NumHisPt,NetType,VesNum,t_plot,all_plot,ODESolver)
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP
PAll1=zeros(NumHisPt,length(all_plot),VesNum);
PAll2=PAll1;UAll=PAll1;QAll=PAll1;AAll=PAll1;ViscAll=PAll1;tAll=PAll1;VesID=PAll1;WfAll=PAll1;WbAll=PAll1;PforwAll=PAll1;PbackwAll=PAll1;UforwAll=PAll1;UbackwAll=PAll1;
MeanP1=zeros(VesNum,NumHisPt);
MeanP2=zeros(VesNum,NumHisPt);
MeanQ=zeros(VesNum,NumHisPt);
MeanU=zeros(VesNum,NumHisPt);
MeanA=zeros(VesNum,NumHisPt);
MeanVisc=zeros(VesNum,NumHisPt);
for j=1:VesNum    % j: Vessel number
  if VesNum == 1
    fileName = [NetType '.his'];
  else
    fileName = [NetType '_' int2str(j) '.his'];
  end
  [t P U Q A Visc VesID]=readHisFile(fileName,NumHisPt,ODESolver);
  if ODESolver==ONED_EXP || ODESolver==ONED_IMP
    tAll(:,:,j)=reshape(t,NumHisPt,length(t)/NumHisPt);
    PAll1(:,:,j)=reshape(P,NumHisPt,length(P)/NumHisPt);
    UAll(:,:,j)=reshape(U,NumHisPt,length(U)/NumHisPt);
    QAll(:,:,j)=reshape(Q,NumHisPt,length(Q)/NumHisPt);
    AAll(:,:,j)=reshape(A,NumHisPt,length(A)/NumHisPt);
    ViscAll(:,:,j)=reshape(Visc,NumHisPt,length(Visc)/NumHisPt);
    for i=1:NumHisPt
      MeanP1(j,i)=mean(PAll1(i,t_plot,j));
      MeanQ(j,i)=mean(QAll(i,t_plot,j));
      MeanU(j,i)=mean(UAll(i,t_plot,j));
      MeanA(j,i)=mean(AAll(i,t_plot,j));
      MeanVisc(j,i)=mean(ViscAll(i,t_plot,j));
    end
  else
    tAll(:,:,j)=reshape(t,NumHisPt,length(t)/NumHisPt);
    QAll(:,:,j)=reshape(P,NumHisPt,length(Q)/NumHisPt);   % P->Q
    PAll1(:,:,j)=reshape(U,NumHisPt,length(P)/NumHisPt);  % U->P1
    PAll2(:,:,j)=reshape(Q,NumHisPt,length(P)/NumHisPt);  % Q->P2
    for i=1:NumHisPt
      MeanP1(j,i)=mean(PAll1(i,t_plot,j));
      MeanP2(j,i)=mean(PAll2(i,t_plot,j));
      MeanQ(j,i)=mean(QAll(i,t_plot,j));
    end
  end
%   disp(j)
  
end