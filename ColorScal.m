% Color scale
function InputRGB=ColorScal(ScaledInput,ColorOrder)
cmap=colormap(hot(ColorOrder));
InputRGB=zeros(length(ScaledInput),3);
MapValue=uint8((ScaledInput-min(ScaledInput))./(max(ScaledInput)-min(ScaledInput))*ColorOrder)+1;
MapValue(MapValue>ColorOrder)=ColorOrder;
InputRGB=cmap(MapValue,:);

