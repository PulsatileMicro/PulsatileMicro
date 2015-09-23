%% WIA
switch NetTypeName
  case 'SymNet'
    PpAll=zeros(NumHisPt,length(t_plot)-1,VesNum);
    PnAll=PpAll;UpAll=PpAll;UnAll=PpAll;
    cnt=0;
    for j=1:VesNum
      %     for j=[1 2 4 8 16 32 40 44 46]
      cnt=cnt+1;
      for i=1:NumHisPt
        WIA_pt=i;
        VesID=j;
        [A beta]=Eval_beta_A(2*sqrt(MeanA(j,i)/1e12/pi),VesParam(4,VesID),VesParam(3,VesID));
        [Pp,Pn,Up,Un]=WIA(PAll(WIA_pt,t_plot,VesID)*133,UAll(WIA_pt,t_plot,VesID)/1e3,AAll(WIA_pt,t_plot,VesID)/1e12,beta,VesParam(4,:),VesParam(3,:),VesParam(6,:));
        Pp=Pp/133;Pn=Pn/133;Up=Up*1e3;Un=Un*1e3;
        PpAll(WIA_pt,:,VesID)=Pp;
        PnAll(WIA_pt,:,VesID)=Pn;
        UpAll(WIA_pt,:,VesID)=Up;
        UnAll(WIA_pt,:,VesID)=Un;
        % 归一化血压
        maxP=max([PAll(WIA_pt,t_plot,VesID) max(PpAll(WIA_pt,:,VesID)) max(PnAll(WIA_pt,:,VesID))]);
        minP=min([PAll(WIA_pt,t_plot,VesID) min(PpAll(WIA_pt,:,VesID)) min(PnAll(WIA_pt,:,VesID))]);
        NormP=(PAll(WIA_pt,t_plot,VesID)-minP)./(maxP-minP);
        NormForwP=(PpAll(WIA_pt,:,VesID)-minP)./(maxP-minP);
        NormBackwP=(PnAll(WIA_pt,:,VesID)-minP)./(maxP-minP);
        % 归一化流速
        %         maxU=max([UAll(WIA_pt,t_plot,VesID) max(UpAll(WIA_pt,:,VesID)) max(UnAll(WIA_pt,:,VesID))]);
        %         minU=min([UAll(WIA_pt,t_plot,VesID) min(UpAll(WIA_pt,:,VesID)) min(UnAll(WIA_pt,:,VesID))]);
        maxU=max([UAll(WIA_pt,t_plot,VesID)]);
        minU=min([UAll(WIA_pt,t_plot,VesID)]);
        NormU=(UAll(WIA_pt,t_plot,VesID)-minU)./(maxU-minU);
        %         NormForwU=(UpAll(WIA_pt,:,VesID)-minU)./(maxU-minU);
        %         NormBackwU=(UnAll(WIA_pt,:,VesID)-minU)./(maxU-minU);
        %     subplot(3,VesNum,cnt+(i-1)*VesNum);hold on;
        subplot(2,NumHisPt,i);hold on;
        %     plot(tAll(1,t_plot,1),PAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
        %     plot(tAll(1,t_plot(1:end-1),1),PpAll(WIA_pt,:,VesID),'r--','LineWidth',2);
        %     plot(tAll(1,t_plot(1:end-1),1),PnAll(WIA_pt,:,VesID),'k:','LineWidth',2);
        plot(tAll(2,t_plot,2),NormP,'b','LineWidth',2);
        %         plot(tAll(2,t_plot(1:end-1),2),NormForwP,'r--','LineWidth',2);
        %         plot(tAll(2,t_plot(1:end-1),2),NormBackwP,'k:','LineWidth',2);
        xlabel('t(s)','FontName','Times New Roman','FontSize',20);
        ylabel('Norm. Pressure','FontName','Times New Roman','FontSize',20);
        set(gca,'FontName','Times New Roman','FontSize',20);
        %         legend('Total','Forward','Backward','FontSize',20);
        switch i
          case 1
            title('始端','FontName','Times New Roman','FontSize',20);
          case 2
            title('中点','FontName','Times New Roman','FontSize',20);
          case 3
            title('末端','FontName','Times New Roman','FontSize',20);
        end
        %     subplot(3,VesNum,j+(i+2)*VesNum);hold on;
        subplot(2,NumHisPt,i+NumHisPt);hold on;
        plot(tAll(2,t_plot,2),NormU,'b','LineWidth',2);
        %         plot(tAll(2,t_plot(1:end-1),2),NormForwU,'r--','LineWidth',2);
        %         plot(tAll(2,t_plot(1:end-1),2),NormBackwU,'k:','LineWidth',2);
      end
    end
  case 'Tree'
    PpAll=zeros(NumHisPt,length(t_plot)-1,VesNum);
    PnAll=PpAll;UpAll=PpAll;UnAll=PpAll;
    for j=1:VesNum
      %     for j=[1 2 4 8 16]
      figure;
      for i=1:NumHisPt
        WIA_pt=i;
        VesID=j;
        [A beta]=Eval_beta_A(2*sqrt(MeanA(j,i)/1e12/pi),VesParam(4,VesID),VesParam(3,VesID));
        [Pp,Pn,Up,Un]=WIA(PAll(WIA_pt,t_plot,VesID)*133,UAll(WIA_pt,t_plot,VesID)/1e3,AAll(WIA_pt,t_plot,VesID)/1e12,beta,VesParam(4,:),VesParam(3,:),VesParam(6,:));
        Pp=Pp/133;Pn=Pn/133;Up=Up*1e3;Un=Un*1e3;
        PpAll(WIA_pt,:,VesID)=Pp;
        PnAll(WIA_pt,:,VesID)=Pn;
        UpAll(WIA_pt,:,VesID)=Up;
        UnAll(WIA_pt,:,VesID)=Un;
        
        subplot(2,NumHisPt,i);hold on;
        plot(tAll(1,t_plot,1),PAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
        plot(tAll(1,t_plot(1:end-1),1),PpAll(WIA_pt,:,VesID),'r--','LineWidth',2);
        plot(tAll(1,t_plot(1:end-1),1),PnAll(WIA_pt,:,VesID),'k:','LineWidth',2);
        xlabel('t(s)','FontName','Times New Roman','FontSize',20);
        ylabel('Pressure(mmHg)','FontName','Times New Roman','FontSize',20);
        set(gca,'FontName','Times New Roman','FontSize',20);
        %         legend('Total','Forward','Backward','FontSize',20);
        switch i
          case 1
            title('始端','FontName','Times New Roman','FontSize',20);
          case 2
            title('中点','FontName','Times New Roman','FontSize',20);
          case 3
            title('末端','FontName','Times New Roman','FontSize',20);
        end
        subplot(2,NumHisPt,i+NumHisPt);hold on;
        plot(tAll(1,t_plot,1),UAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
        plot(tAll(1,t_plot(1:end-1),1),UpAll(WIA_pt,:,VesID),'r--','LineWidth',2);
        plot(tAll(1,t_plot(1:end-1),1),UnAll(WIA_pt,:,VesID),'k:','LineWidth',2);
        xlabel('t(s)','FontName','Times New Roman','FontSize',20);
        ylabel('Velocity(mm/s)','FontName','Times New Roman','FontSize',20);
        set(gca,'FontName','Times New Roman','FontSize',20);
      end
    end
  case 'Single'
    for i=1:NumHisPt
      %     for i=18:
      WIA_pt=i;
      VesID=1;
      [A beta]=Eval_beta_A(2*sqrt(MeanA(i)/1e12/pi),VesParam(4,VesID),VesParam(3,VesID));
      [Pp,Pn,Up,Un]=WIA(PAll(WIA_pt,t_plot,VesID)*133,UAll(WIA_pt,t_plot,VesID)/1e3,AAll(WIA_pt,t_plot,VesID)/1e12,beta,VesParam(4,:),VesParam(3,:),VesParam(7,:));
      Pp=Pp/133;Pn=Pn/133;Up=Up*1e3;Un=Un*1e3;
      PpAll(:,i)=Pp;PnAll(:,i)=Pn;UpAll(:,i)=Up;UnAll(:,i)=Un;
      % 归一化血压
      maxP=max([PAll(WIA_pt,t_plot,VesID) max(PpAll(:,i)) max(PnAll(:,i))]);
      minP=min([PAll(WIA_pt,t_plot,VesID) min(PpAll(:,i)) min(PnAll(:,i))]);
      NormP=(PAll(WIA_pt,t_plot,VesID)-minP)./(maxP-minP);
      NormForwP=(PpAll(:,i)-minP)./(maxP-minP);
      NormBackwP=(PnAll(:,i)-minP)./(maxP-minP);
      subplot(2,NumHisPt,i);hold on;
      plot(tAll(2,t_plot),PAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
      %       plot(tAll(2,t_plot(1)+indp(WIA_pt)),PAll(WIA_pt,t_plot(1)+indp(WIA_pt),VesID),'ro');
      plot(tAll(2,t_plot(1:end-1)),PpAll(:,i),'r--','LineWidth',2);
      plot(tAll(2,t_plot(1:end-1)),PnAll(:,i),'k:','LineWidth',2);
      %       plot(tAll(2,t_plot),NormP,'b','LineWidth',2);
      %       plot(tAll(2,t_plot(1:end-1)),NormForwP,'r--','LineWidth',2);
      %       plot(tAll(2,t_plot(1:end-1)),NormBackwP,'k:','LineWidth',2);
      %       ylim([0 1]);
      xlabel('t(s)','FontName','Times New Roman','FontSize',20);
      ylabel('Pressure(mmHg)','FontName','Times New Roman','FontSize',20);
      set(gca,'FontName','Times New Roman','FontSize',20);
      %   legend('Total','Forward','Backward','FontSize',20);
      switch i
        case 1
          title('始端','FontName','Times New Roman','FontSize',20);
        case 2
          title('中点','FontName','Times New Roman','FontSize',20);
        case 3
          title('末端','FontName','Times New Roman','FontSize',20);
      end
      
      % 归一化流速
      maxU=max([UAll(WIA_pt,t_plot,VesID) max(UpAll(:,i)) max(UnAll(:,i))]);
      minU=min([UAll(WIA_pt,t_plot,VesID) min(UpAll(:,i)) min(UnAll(:,i))]);
      NormU=(UAll(WIA_pt,t_plot,VesID)-minU)./(maxU-minU);
      NormForwU=(UpAll(:,i)-minU)./(maxU-minU);
      NormBackwU=(UnAll(:,i)-minU)./(maxU-minU);
      subplot(2,NumHisPt,NumHisPt+i);hold on;
      plot(tAll(2,t_plot),UAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
      %       plot(tAll(2,t_plot(1)+indu(WIA_pt)),UAll(WIA_pt,t_plot(1)+indu(WIA_pt),VesID),'ro');
      plot(tAll(2,t_plot(1:end-1)),UpAll(:,i),'r--','LineWidth',2);
      plot(tAll(2,t_plot(1:end-1)),UnAll(:,i),'k:','LineWidth',2);
      %       plot(tAll(2,t_plot),NormU,'b','LineWidth',2);
      %       plot(tAll(2,t_plot(1:end-1)),NormForwU,'r','LineWidth',2);
      %   NormBackwU=2*max(NormBackwU)-NormBackwU;
      %       plot(tAll(2,t_plot(1:end-1)),NormBackwU,'k','LineWidth',2);
      xlabel('t(s)','FontName','Times New Roman','FontSize',20);
      ylabel('Velocity(mm/s)','FontName','Times New Roman','FontSize',20);
      set(gca,'FontName','Times New Roman','FontSize',20);
    end
    %     PI_P=PI_P';
    %     PI_U=PI_U';
  case 'Junc'
    PpAll=zeros(NumHisPt,length(t_plot)-1,VesNum);
    PnAll=PpAll;UpAll=PpAll;UnAll=PpAll;
    figure;
    cnt=0;
    for j=1:1
      cnt=cnt+1;
      for i=1:NumHisPt
        WIA_pt=i;
        VesID=j;
        [A beta]=Eval_beta_A(2*sqrt(MeanA(j,i)/1e12/pi),VesParam(4,VesID),VesParam(3,VesID));
        [Pp,Pn,Up,Un]=WIA(PAll(WIA_pt,t_plot,VesID)*133,UAll(WIA_pt,t_plot,VesID)/1e3,AAll(WIA_pt,t_plot,VesID)/1e12,beta,VesParam(4,:),VesParam(3,:),VesParam(6,:));
        Pp=Pp/133;Pn=Pn/133;Up=Up*1e3;Un=Un*1e3;
        PpAll(WIA_pt,:,VesID)=Pp;
        PnAll(WIA_pt,:,VesID)=Pn;
        UpAll(WIA_pt,:,VesID)=Up;
        UnAll(WIA_pt,:,VesID)=Un;
        subplot(2,NumHisPt,i);hold on;
        plot(tAll(2,t_plot,2),PAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
        plot(tAll(2,t_plot(1:end-1),2),Pp,'r--','LineWidth',2);
        plot(tAll(2,t_plot(1:end-1),2),Pn,'k:','LineWidth',2);
        xlabel('t(s)','FontName','Times New Roman','FontSize',20);
        ylabel('Pressure(mmHg)','FontName','Times New Roman','FontSize',20);
        set(gca,'FontName','Times New Roman','FontSize',20);
        %         legend('Total','Forward','Backward');
        xlim([4 4.8]);
        switch i
          case 1
            title('始端','FontName','Times New Roman','FontSize',20);
          case 2
            title('中点','FontName','Times New Roman','FontSize',20);
          case 3
            title('末端','FontName','Times New Roman','FontSize',20);
        end
        subplot(2,NumHisPt,i+NumHisPt);hold on;
        plot(tAll(2,t_plot,2),UAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
        plot(tAll(2,t_plot(1:end-1),2),Up,'r--','LineWidth',2);
        plot(tAll(2,t_plot(1:end-1),2),Un,'k:','LineWidth',2);
        xlabel('t(s)','FontName','Times New Roman','FontSize',20);
        ylabel('Velocity(mm/s)','FontName','Times New Roman','FontSize',20);
        set(gca,'FontName','Times New Roman','FontSize',20);
        xlim([4 4.8]);
      end
    end
end