segType=zeros(VesNum,1);
%818
% segType([216;544;547;559;561;47;57;70;80;143;681;683])=[3;ones(11,1)];
%Egg_571
segType([167;168;143;144;560;561;562;212;580;581;583;410;615;449;514;452])=ones(16,1);
segType([213;544;506;507;583;394;384;385;386;396])=2*ones(10,1);
segType=autoVesTypeAdv(segType,FromNeg,ToNeg,Porder);
segType(424)=3;
segType(428)=2;
%% 血管网络画图 %%%%
% drawType:1 - 血管参数 2 - 血管类型
drawType=2;
% segType=SegType;
switch num2str(drawType)
    case '1'
        plotForGeneralNetwork(MeanP,VesType,drawType);
    case '2'
        plotForGeneralNetwork(segType,VesType,drawType);
    otherwise
        error('Wrong setting of "drawType" !');
end