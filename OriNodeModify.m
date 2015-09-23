function [From,To]=OriNodeModify(Boundary,From,To)
From(:,2)=From(:,1);
To(:,2)=To(:,1);
for i=1:length(Boundary(:,1))
  if Boundary(i,3)>=0 && Boundary(i,2)~=0
    InputIndex=find(From(:,1)==Boundary(i,1), 1);
    if isempty(InputIndex)
      OutputIndex=find(To(:,1)==Boundary(i,1));
      From(OutputIndex,1)=To(OutputIndex,2);
      To(OutputIndex,1)=From(OutputIndex,2);
    end
  else
    OutputIndex=find(To(:,1)==Boundary(i,1), 1);
    if isempty(OutputIndex)
      InputIndex=find(From(:,1)==Boundary(i,1));
      From(InputIndex,1)=To(InputIndex,2);
      To(InputIndex,1)=From(InputIndex,2);
    end
  end
end
From(:,2)=[];
To(:,2)=[];