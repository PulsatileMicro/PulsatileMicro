function [signal,p_time,p_press,u_time,u_press]=findpu(original,f,flag)
%以下内容是针对当输入的序列不是标准状态而设置，如实验提供的序列就是上下倒转的。所以应该设置flag为0，
%但如果输入序列是标准状态下的序列，即没有上下倒置，没有左右倒置，则只需将flag设置为非0实数便可。
% Developed by Wen,Haoxiang
if flag==0
  m=mean(original);
  original=original-m;
  original=-original;
  original=original+m;
end
y=original;
n=length(y);
n1=0.16*f;
% n1=0.12*f;
windows_width=round(n1/2-1);%windows_width是半个窗口的宽度，前半个窗口，后半个窗口都是windows_width
clear n1;
for j=1+windows_width:n-windows_width,
  movewindows_y(j)=0;
  for i=(-1)*windows_width:windows_width,
    movewindows_y(j)=movewindows_y(j)+y(i+j);
  end
end
clear i;clear j;
%
%以下程序用于对移动窗积分后的信号进行微分。
movewindows_y=movewindows_y/(2*windows_width+1);
n1=1:n-windows_width-1;
difference_y(n1)=movewindows_y(n1+1)-movewindows_y(n1);%difference_y就是移动窗积分再微分后的信号。
clear n1;
signal=y;%备份原信号，最后输出时返回的是只去除基线漂移，没有进行移动窗积分及移动窗积分以后操作的信号。
%以下程序是因为移动窗积分后，头39点数据和尾39点数据没用，裁剪掉。
y=y(1+windows_width:n-windows_width-1);
difference_y=difference_y(:,1+windows_width:n-windows_width-1);
%
%以下程序是将微分后的信号处理，若大于0就平方之，如果少于0，就屏蔽掉经过这个处理后的信号是y4，
n1=length(difference_y);
for i=1:n1,
  if difference_y(i)>0
    y4(i)=difference_y(i)*difference_y(i);
  else
    y4(i)=0;
  end
end
%
%以下程序是打印出微分后信号与原信号的波形，观察U波与处理后的信号在时间域中的关系。
n1=1:n-2*windows_width-1;
m=mean(y);
y2=y-m;
% figure;plot(n1,y4,'r',n1,y2,'g');
clear y2;clear m;clear n1; %这4变量，以及上面的程序都是为了打印两个信号的关系
%
%设定阈值。以y4数学期望的1/n为阈值。n推荐为6。太小会将U波前面的一个小重搏波
%误认为U波，n太大会漏检。
threshold=0;
for n1=1:n-2*windows_width-1,
  threshold=threshold+y4(n1);
end
threshold=threshold/n1/6;
%
%以下是打印出阈值与y4信号的关系。为设定，调整阈值提供依据。
% figure;
n1=1:n-2*windows_width-1;
% figure;plot(n1,threshold,'r',n1,y4,'g');
new_array_index=0;%new_array_index为U点压强，时间序列的下标，找到一个U点自增1。
search_front=0;%必要时可调整为round(0.1*f);
% search_front=round(0.1*f);
search_last=round(0.16*f*1.25);
%search_front为达到阈值时开始搜索点的上限；search_last为下限
%search_last没什么问题。search_front要小心调整。太小找到的是U点后面的点，
%太大会导致误认方U波重的小重搏波为U波。
for a=5:round(0.5*f*1.75),
  if ((y4(a)==0)&(y4(a-1)==0)&(y4(a-2)==0)&(y4(a-3)==0)&(y4(a-4)==0))
    %     if ((y4(a)==0)&(y4(a-1)==0)&(y4(a-2)==0))
    i=a;
    break;
  end
end
clear a;

while i<n1(n-2*windows_width-1)-search_last-2,%开始对y4搜索，目的是寻找大于阈值的点
  i=i+1;
  if y4(i)>threshold%如果找到y4大于阈值的，开始在此点的附近寻找原信号的最低点，即U波
    new_array_index=new_array_index+1;  %因为U波必在这附近，所以必定能寻找到。所以下标加1
    if i-search_front<0   %防止刚开始搜索时，原信号的下标小于0
      ij=2;
    else
      ij=i-search_front+1;%防止搜索到信号的未端时，下标大于信号的点数。
    end
    max_value=y(ij-1);
    time_value=ij-1;
    %搜索最低点，并记下该点的值及坐标
    for j=ij:i+search_last,
      if max_value>y(j)
        max_value=y(j);
        time_value=j;
      end
    end
    u_press(new_array_index)=max_value;%u_press用于记录该点的值
    u_time(new_array_index)=time_value;%u_time用于记录该点的下标值，即时间
    i=time_value+round(0.4*f*1.25);%因为寻找到U波后的一段时间内，不可能再出现U波，所以设置一段死区
  end
end
clear i;clear j;clear ij;clear search_front;clear search_last;clear time_value;
%

%以下程序用于搜索最高点，即P波。思想是两个U波夹一个P波，在两个P波间找到一个最高点。便是P波
for i=1:new_array_index-1,
  search_front=u_time(i);
  search_last=round((u_time(i+1)+u_time(i))/2);
  min_value=y(search_front);
  time_value=search_front;
  for j=search_front+1:search_last,
    if min_value<y(j)
      min_value=y(j);
      time_value=j;
    end
    p_time(i)=time_value;
    p_press(i)=min_value;
  end
end
%
%以下程序是因为一开始的数据裁头裁尾了即把移动窗半个窗口的长度给裁掉了。所以现在把裁掉的数据补上。
p_time=p_time+windows_width;
u_time=u_time+windows_width;
n=length(p_press);
% clear u_time;
% clear u_press;
% for i=1:n,
%     search_last=p_time(i)-round(0.01*f);
%     if p_time(i)-round(0.16*f)<=6
%         search_front=6;
%     else
%         search_front=p_time(i)-round(0.16*f*2);
%     end
%     for j=search_last:-1:search_front,
%         if ((signal(j)<=signal(j-1))&(signal(j-1)<=signal(j-2))&(signal(j-2)<=signal(j-3))&(signal(j-3)<=signal(j-4)))%&(signal(j-4)<=signal(j-5)))
%             u_time(i)=j;
%             u_press(i)=signal(j);
%             break;
%         end
%     end
% end
%打印出原信号，最高点（P波），最低点（U波）标记出来。
% figure;
n1=1:length(signal);
% plot(n1,signal,'r',p_time,p_press,'*',u_time,u_press,'^');