function [  ] = fPlotCodeCompBEM(field, BEM,BEMNoLoss,BEMNewTip,BEMShen,BEMXS,R)
if(R==0)
    R=1;
end
global legds;
if(isfield(BEMNoLoss,field) && isfield(BEMNoLoss,'r'))
    plot(BEMNoLoss.r/R,getfield(BEMNoLoss,field),'Color',[0.4 0.4 0.4],'LineWidth',2.2)
    legds{end+1}='BEMNoLoss';
end
if(~isempty(BEMNewTip))
    if(isfield(BEMNewTip,field) && isfield(BEMNewTip,'r'))
        plot(BEMNewTip.r/R,getfield(BEMNewTip,field),'-','Color','k')
        legds{end+1}='BEM New'; 
    end
end
if(isfield(BEM,field) && isfield(BEM,'r'))
    plot(BEM.r/R,getfield(BEM,field),'--','Color',fColrs(1))
    legds{end+1}='BEM';
end
if(isfield(BEMShen,field) && isfield(BEMShen,'r'))
    plot(BEMShen.r/R,getfield(BEMShen,field),':','Color',fColrs(1))
    legds{end+1}='BEM Shen';
end

if(isfield(BEMXS,field) && isfield(BEMXS,'r'))
    plot(BEMXS.r/R,getfield(BEMXS,field),':','Color',fColrs(2))
    legds{end+1}='BEM Xu Sankar';
end

if(R~=1)
    xlabel('$r/R$ [.]')
else
    xlabel('$r$ [m]')
end
end

