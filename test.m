figure;
for i=1:5
  subplot(5,3,(i-1)*3+1);
  [AX,H1,H2]=plotyy(time,TypFlow_Init(:,i),time,NormTypFlow_Init(:,i));
  xlabel('Time (s)');
  set(get(AX(1),'Ylabel'),'String','Flow Rate (nL/min)')
  set(get(AX(2),'Ylabel'),'String','Normalized Flow Rate')
  
  subplot(5,3,(i-1)*3+2);
  [AX,H1,H2]=plotyy(time,TypFlow_C01(:,i),time,NormTypFlow_C01(:,i));
  xlabel('Time (s)');
  set(get(AX(1),'Ylabel'),'String','Flow Rate (nL/min)')
  set(get(AX(2),'Ylabel'),'String','Normalized Flow Rate')
  
  subplot(5,3,i*3);
  [AX,H1,H2]=plotyy(time,TypFlow_R01(:,i),time,NormTypFlow_R01(:,i));
  xlabel('Time (s)');
  set(get(AX(1),'Ylabel'),'String','Flow Rate (nL/min)')
  set(get(AX(2),'Ylabel'),'String','Normalized Flow Rate')
end