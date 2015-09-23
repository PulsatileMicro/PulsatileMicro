%%m行；n列；l左边距；r右边距；b底边距；t上边距；mm行距；nn列距
%Output(l,4):l行数，表示subplot的编号；每一行的四个参数即subplot的position参数
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
        