function [figuredata,figh]=plotVYBdata(figuredata,idx,alpha,setYScaletoLog,ligandname,colormat)

addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/flowanalysis
addpath ~/Documents/MATLAB/FACS/
addpath ~/Documents/MATLAB/barwitherr


for k=1:length(figuredata)
    fdn=fieldnames(figuredata(k));
    orilen=length(figuredata(k).samplenames);
    for j=1:length(fdn)
        if length(figuredata(k).(fdn{j}))==orilen
            vec=figuredata(k).(fdn{j});
%             fprintf('%s\n',fdn{j});
            figuredata(k).(fdn{j})=vec(idx);
        end
    end
    m=[figuredata(k).meannorm_minus' figuredata(k).meannorm_plus']+2;
s=[figuredata(k).signorm_minus' figuredata(k).signorm_plus'];

fold=figuredata(k).meanfold;
numsamps=[figuredata(k).numsamps' figuredata(k).numsamps'];
xlabelnames=figuredata(k).samplenames;

if setYScaletoLog
    err=[];
    err(:,:,1)=10.^(m)-10.^(m-s);
    err(:,:,2)=10.^(m+s)-10.^(m);
else
    err=[];
    err(:,:,1)=(10.^(m)-10.^(m-s))./sqrt(numsamps);
    err(:,:,2)=(10.^(m+s)-10.^(m))./sqrt(numsamps);
end
    
mm=10.^m;
% setfig(sprintf('%s %s',figuredata(k).condition,ligandname));clf
[h hErrorbar]=barwitherr(err,mm);


set(h,'linewidth',2)
set(h(1),'FaceColor',[0.6 0.6 0.6])
if nargin<6
    set(h(2),'FaceColor',[0.9 0.2 0.3])
else
    set(h(2),'FaceColor',colormat)
end

h(1).BaseLine.LineWidth=2;


L = get(gca,'XLim');
% xlabelnames=uniquenames(in);
if length(xlabelnames)>10
set(gca,'xticklabel',{'',xlabelnames{:},''})
NumTicks = length(xlabelnames)+2; 
    
elseif  length(xlabelnames)>5
set(gca,'xticklabel',{'',xlabelnames{:},''})
NumTicks = length(xlabelnames)+2; 
%     
% elseif length(xlabelnames)>5
% set(gca,'xticklabel',{'',xlabelnames{:},'',''})
% NumTicks = length(xlabelnames)+3; 

elseif length(xlabelnames)==5 % for zeatin
set(gca,'xticklabel',{'','',xlabelnames{(1)},'',xlabelnames{(2)},'',...
                      xlabelnames{(3)},'',xlabelnames{(4)},'',...
                      xlabelnames{(5)},'',''})
NumTicks = length(xlabelnames)+8; 


elseif length(xlabelnames)==4 % for zeatin
set(gca,'xticklabel',{'',xlabelnames{(1)},'',xlabelnames{(2)},'',...
                      xlabelnames{(3)},'',xlabelnames{(4)},''})
NumTicks = length(xlabelnames)+5; 

else
set(gca,'xticklabel',{xlabelnames{:}})
NumTicks = length(xlabelnames); 
end
set(gca,'TickLabelInterpreter','none')
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'fontsize',22)
set(gca,'linewidth',2)
YL=get(gca,'Ylim');


f=find(~isnan(fold));
foldlabel={};
for i=1:length(f)
    foldlabel{end+1}=sprintf('%0.1f',fold(f(i)));
end

if  ~setYScaletoLog
txtheight=10.^(m(:,2)')+9.5;
axis([0 length(xlabelnames)+1 0 max(YL)])
else 
txtheight=10.^(m(:,2)')*1.5;
axis([0 length(xlabelnames)+1 min(min(10.^m))*0.9 max(YL)*2.8])
end

txtwidth=(1:length(xlabelnames))-0.35;
text(txtwidth(f),txtheight(f),foldlabel,'fontsize',22)

pval=figuredata(k).pval((f));
hval=zeros(1,length(pval));

hval(pval<alpha & pval>alpha/10)=1;
hval(pval<(alpha/10) & pval>(alpha/100))=2;
hval(pval<(alpha/100) & pval>(alpha/1000))=3;
hval(pval<(alpha/1000))=4;


for j=1:length(hval)
    if ~setYScaletoLog
        astoffset=5;
       if hval(j)==1
           text(txtwidth(f(j)),txtheight(f(j))+astoffset,' *','fontsize',28)
       elseif hval(j)==2
           text(txtwidth(f(j))+0.05,txtheight(f(j))+astoffset,'**','fontsize',28)
       elseif hval(j)==3
           text(txtwidth(f(j))-0.1,txtheight(f(j))+astoffset,'***','fontsize',28)
       elseif hval(j)==4
           text(txtwidth(f(j))-0.2,txtheight(f(j))+astoffset,'****','fontsize',28)
       end
   
    else
        astoffset=1.3;
       if hval(j)==1
           text(txtwidth(f(j)),txtheight(f(j))*astoffset,' *','fontsize',28)
       elseif hval(j)==2
           text(txtwidth(f(j))+0.05,txtheight(f(j))*astoffset,'**','fontsize',28)
       elseif hval(j)==3
           text(txtwidth(f(j))-0.1,txtheight(f(j))*astoffset,'***','fontsize',28)
       elseif hval(j)==4
           text(txtwidth(f(j))-0.2,txtheight(f(j))*astoffset,'****','fontsize',28)
       end        
    end
end


legendnames={sprintf('0 mM %s',ligandname),sprintf('%s %s',figuredata(k).condition,ligandname)};
if setYScaletoLog
legend(legendnames,'location','northoutside')
else
    legend(legendnames,'location','northwest')
    legend(legendnames,'location','northoutside')

end

% legend('0 mM theophylline', '1 mM theophylline','location','southeast')
% legend('0 mM aciclovir', '5 mM aciclovir','location','southeast')
% legend('0 \muM noscapine', '500 \muM noscapine','location','northwest')
if setYScaletoLog
    set(gca,'YScale','log')
end
ylabel('GFP/mCherry')
end

end