function [newu_time,newu_press]=Get_new_u(y,u_time,u_press,p_time,p_press)
n=length(p_time);
for i=1:n,   
    diff=p_press(i)-u_press(i);
    min=u_press(i)+round(diff/4);
    max=u_press(i)+round(diff/2);
    ij=0;
    clear y1;
    clear x1;
    for j=u_time(i):p_time(i)
        if ((y(j)>=min)&(y(j)<=max))
            ij=ij+1;
            y1(ij)=y(j);
            x1(ij)=j;
        end
    end   
    a=polyfit(x1,y1,1);
    newu_time(i)=-a(2)/a(1);
    newu_press(i)=0;
end
figure;
n=1:length(y);
plot(n,y,'b',newu_time,newu_press,'*');

    
