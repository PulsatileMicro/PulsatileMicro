function [ Output ] = Norm0_1( Input )
Output=(Input-min(Input))./(max(Input)-min(Input));
end

