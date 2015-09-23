function [signal,p_time,p_press,u_time,u_press]=findpu(original,f,flag)
%������������Ե���������в��Ǳ�׼״̬�����ã���ʵ���ṩ�����о������µ�ת�ġ�����Ӧ������flagΪ0��
%��������������Ǳ�׼״̬�µ����У���û�����µ��ã�û�����ҵ��ã���ֻ�轫flag����Ϊ��0ʵ����ɡ�
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
windows_width=round(n1/2-1);%windows_width�ǰ�����ڵĿ�ȣ�ǰ������ڣ��������ڶ���windows_width
clear n1;
for j=1+windows_width:n-windows_width,
  movewindows_y(j)=0;
  for i=(-1)*windows_width:windows_width,
    movewindows_y(j)=movewindows_y(j)+y(i+j);
  end
end
clear i;clear j;
%
%���³������ڶ��ƶ������ֺ���źŽ���΢�֡�
movewindows_y=movewindows_y/(2*windows_width+1);
n1=1:n-windows_width-1;
difference_y(n1)=movewindows_y(n1+1)-movewindows_y(n1);%difference_y�����ƶ���������΢�ֺ���źš�
clear n1;
signal=y;%����ԭ�źţ�������ʱ���ص���ֻȥ������Ư�ƣ�û�н����ƶ������ּ��ƶ��������Ժ�������źš�
%���³�������Ϊ�ƶ������ֺ�ͷ39�����ݺ�β39������û�ã��ü�����
y=y(1+windows_width:n-windows_width-1);
difference_y=difference_y(:,1+windows_width:n-windows_width-1);
%
%���³����ǽ�΢�ֺ���źŴ���������0��ƽ��֮���������0�������ε���������������ź���y4��
n1=length(difference_y);
for i=1:n1,
  if difference_y(i)>0
    y4(i)=difference_y(i)*difference_y(i);
  else
    y4(i)=0;
  end
end
%
%���³����Ǵ�ӡ��΢�ֺ��ź���ԭ�źŵĲ��Σ��۲�U���봦�����ź���ʱ�����еĹ�ϵ��
n1=1:n-2*windows_width-1;
m=mean(y);
y2=y-m;
% figure;plot(n1,y4,'r',n1,y2,'g');
clear y2;clear m;clear n1; %��4�������Լ�����ĳ�����Ϊ�˴�ӡ�����źŵĹ�ϵ
%
%�趨��ֵ����y4��ѧ������1/nΪ��ֵ��n�Ƽ�Ϊ6��̫С�ὫU��ǰ���һ��С�ز���
%����ΪU����n̫���©�졣
threshold=0;
for n1=1:n-2*windows_width-1,
  threshold=threshold+y4(n1);
end
threshold=threshold/n1/6;
%
%�����Ǵ�ӡ����ֵ��y4�źŵĹ�ϵ��Ϊ�趨��������ֵ�ṩ���ݡ�
% figure;
n1=1:n-2*windows_width-1;
% figure;plot(n1,threshold,'r',n1,y4,'g');
new_array_index=0;%new_array_indexΪU��ѹǿ��ʱ�����е��±꣬�ҵ�һ��U������1��
search_front=0;%��Ҫʱ�ɵ���Ϊround(0.1*f);
% search_front=round(0.1*f);
search_last=round(0.16*f*1.25);
%search_frontΪ�ﵽ��ֵʱ��ʼ����������ޣ�search_lastΪ����
%search_lastûʲô���⡣search_frontҪС�ĵ�����̫С�ҵ�����U�����ĵ㣬
%̫��ᵼ�����Ϸ�U���ص�С�ز���ΪU����
for a=5:round(0.5*f*1.75),
  if ((y4(a)==0)&(y4(a-1)==0)&(y4(a-2)==0)&(y4(a-3)==0)&(y4(a-4)==0))
    %     if ((y4(a)==0)&(y4(a-1)==0)&(y4(a-2)==0))
    i=a;
    break;
  end
end
clear a;

while i<n1(n-2*windows_width-1)-search_last-2,%��ʼ��y4������Ŀ����Ѱ�Ҵ�����ֵ�ĵ�
  i=i+1;
  if y4(i)>threshold%����ҵ�y4������ֵ�ģ���ʼ�ڴ˵�ĸ���Ѱ��ԭ�źŵ���͵㣬��U��
    new_array_index=new_array_index+1;  %��ΪU�������⸽�������Աض���Ѱ�ҵ��������±��1
    if i-search_front<0   %��ֹ�տ�ʼ����ʱ��ԭ�źŵ��±�С��0
      ij=2;
    else
      ij=i-search_front+1;%��ֹ�������źŵ�δ��ʱ���±�����źŵĵ�����
    end
    max_value=y(ij-1);
    time_value=ij-1;
    %������͵㣬�����¸õ��ֵ������
    for j=ij:i+search_last,
      if max_value>y(j)
        max_value=y(j);
        time_value=j;
      end
    end
    u_press(new_array_index)=max_value;%u_press���ڼ�¼�õ��ֵ
    u_time(new_array_index)=time_value;%u_time���ڼ�¼�õ���±�ֵ����ʱ��
    i=time_value+round(0.4*f*1.25);%��ΪѰ�ҵ�U�����һ��ʱ���ڣ��������ٳ���U������������һ������
  end
end
clear i;clear j;clear ij;clear search_front;clear search_last;clear time_value;
%

%���³�������������ߵ㣬��P����˼��������U����һ��P����������P�����ҵ�һ����ߵ㡣����P��
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
%���³�������Ϊһ��ʼ�����ݲ�ͷ��β�˼����ƶ���������ڵĳ��ȸ��õ��ˡ��������ڰѲõ������ݲ��ϡ�
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
%��ӡ��ԭ�źţ���ߵ㣨P��������͵㣨U������ǳ�����
% figure;
n1=1:length(signal);
% plot(n1,signal,'r',p_time,p_press,'*',u_time,u_press,'^');