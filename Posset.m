%%m�У�n�У�l��߾ࣻr�ұ߾ࣻb�ױ߾ࣻt�ϱ߾ࣻmm�оࣻnn�о�
%Output(l,4):l��������ʾsubplot�ı�ţ�ÿһ�е��ĸ�������subplot��position����
function Output=Posset(m,n,l,r,b,t,mm,nn)
Output=zeros(m*n,4);
Width=(1-l-r-(n-1)*nn)/n;
Height=(1-b-t-(m-1)*mm)/m;
count=0;
for i=1:m
    for j=1:n
        count=count+1;
        x=l+(j-1)*(Width+nn);
        y=1-t-Height-(i-1)*(Height+mm);
        Output(count,:)=[x y Width Height];
    end
end
        