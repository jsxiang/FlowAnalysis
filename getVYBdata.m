function [figuredata,figh]=getVYBdata(analyzedALL,expligandconc,ligandname,time,alpha,normsTRSVctl,setYScaletoLog,domedian,plotindividualsamples,switches,timecourse,timeassay);

addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/flowanalysis
addpath ~/Documents/robot/MatLab-utilities
addpath ~/Documents/MATLAB/FACS/

grouped=struct;

for k=1:length(analyzedALL)
    nameroots={};
    a=analyzedALL(k);
    for i=1:length(a.samplenames)
        re=regexpi(a.samplenames{i},'(\w*:?\w*)(-.*)','tokens');
        nameroots{end+1}=re{1}{1};

    end

    uniquenames=unique(nameroots);
    grouped(k).time=struct;

    for i=1:length(uniquenames)
        hasname=regexp(a.samplenames,strcat(uniquenames(i),'-'));
        hasnamei=find(~cellfun('isempty',hasname));
        grouped(k).time(i).nameroot=uniquenames{i};
        grouped(k).time(i).allnames=a.samplenames(hasnamei);

        if timecourse
            nameext={};
            for p=1:length(grouped(k).time(i).allnames)
                re=regexpi(grouped(k).time(i).allnames{p},'([A-Z]+[0-9]*.*)(_.*)','tokens');
                nameext{end+1}=re{1}{1};
            end
            grouped(k).time(i).nameext=nameext;
        end

        grouped(k).time(i).mu=a.mu(hasnamei);
        grouped(k).time(i).medF1=a.medFP1(hasnamei);
        grouped(k).time(i).sigma=a.sigma(hasnamei);
        grouped(k).time(i).medsigFP1=a.medsigFP1(hasnamei);
        grouped(k).time(i).omit=a.omit(hasnamei);        

        if switches
            n=regexpi(grouped(k).time(i).allnames,'(.*)(-lig)(.*)');
            nolig=~cellfun('isempty',n);
            w=regexpi(grouped(k).time(i).allnames,'(.*)(+lig)(.*)');
            wlig=~cellfun('isempty',w);
            grouped(k).time(i).noligmu=grouped(k).time(i).mu(nolig);
            grouped(k).time(i).wligmu=grouped(k).time(i).mu(wlig);
            grouped(k).time(i).noligmedF1=grouped(k).time(i).medF1(nolig);
            grouped(k).time(i).wligmedF1=grouped(k).time(i).medF1(wlig);
            grouped(k).time(i).nolig=nolig;
            grouped(k).time(i).wlig=wlig;
            grouped(k).time(i).fold=grouped(k).time(i).wligmu-grouped(k).time(i).noligmu;
            grouped(k).time(i).foldmedF1=grouped(k).time(i).wligmedF1-grouped(k).time(i).noligmedF1;
            
            nameswitch={};
            for p=1:length(grouped(k).time(i).allnames)
                re=regexpi(grouped(k).time(i).allnames{p},'(.*)([-|+]lig)(.*)','tokens');
                nameswitch{end+1}=re{1}{1};
            end
            grouped(k).time(i).nameswitch=nameswitch;
            
        end


        
    end


    for i=1:length(grouped(k).time)
        isbad=regexp(grouped(k).time(i).omit,'bad');
        isgood=cellfun('isempty',isbad);
        
        if switches
            goodmu=[grouped(k).time(i).mu(find(grouped(k).time(i).nolig .* isgood));
                    grouped(k).time(i).mu(find(grouped(k).time(i).wlig .* isgood))]';
            goodmedF1=[grouped(k).time(i).medF1(find(grouped(k).time(i).nolig .* isgood));
                    grouped(k).time(i).medF1(find(grouped(k).time(i).wlig .* isgood))]';
            grouped(k).time(i).goodmumean=mean(goodmu);
            grouped(k).time(i).goodmusigma=std(goodmu);
            grouped(k).time(i).goodmedF1mean=mean(goodmedF1);
            grouped(k).time(i).goodmedF1sigma=std(goodmedF1);
            try
            grouped(k).time(i).foldmean=grouped(k).time(i).goodmumean(2)-grouped(k).time(i).goodmumean(1);
            catch
                grouped(k).time(i).foldmean=nan;
            end
            
            try
                grouped(k).time(i).foldmedF1mean=grouped(k).time(i).goodmedF1mean(2)-grouped(k).time(i).goodmedF1mean(1);
            catch
                grouped(k).time(i).foldmedF1mean=nan;
            end
            goodsigma=[grouped(k).time(i).sigma(find(grouped(k).time(i).nolig .* isgood));
                    grouped(k).time(i).sigma(find(grouped(k).time(i).wlig .* isgood))]';
                
            goodmedsigFP1=[grouped(k).time(i).medsigFP1(find(grouped(k).time(i).nolig .* isgood));
                    grouped(k).time(i).medsigFP1(find(grouped(k).time(i).wlig .* isgood))]';
        else
            goodmu=grouped(k).time(i).mu(find(isgood));
            goodmedF1=grouped(k).time(i).medF1(find(isgood));
            goodsigma=grouped(k).time(i).sigma(find(isgood));
            goodmedsigFP1=grouped(k).time(i).medsigFP1(find(isgood));
            grouped(k).time(i).goodmumean=mean(goodmu);
            grouped(k).time(i).goodmusigma=std(goodmu);
            grouped(k).time(i).goodmedF1mean=mean(goodmedF1);
            grouped(k).time(i).goodmedF1sigma=std(goodmedF1);
        
        end

        if switches
            goodnames=grouped(k).time(i).nameswitch(find(isgood));
            goodnames=goodnames(1:(length(goodnames)/2));
        elseif timecourse
            goodnames=grouped(k).time(i).nameext(find(isgood));

        else
            goodnames=grouped(k).time(i).allnames(find(isgood));
        end

        if plotindividualsamples
        setfig(grouped(k).time(i).nameroot);
        if k==1
            clf
        end
        hold on
        size(goodmu);
        plot(goodmu,'-o','MarkerSize',20,'linewidth',4)

        L = get(gca,'XLim');
        xlabelnames=goodnames;
        if length(goodnames)>5
        set(gca,'xticklabel',{xlabelnames{:},''});
        NumTicks = length(xlabelnames); 
            
        elseif length(goodnames)>4
        set(gca,'xticklabel',{'',xlabelnames{:},''});
        NumTicks = length(xlabelnames)+1; 
        else
        set(gca,'xticklabel',{xlabelnames{:},''});
        NumTicks = length(xlabelnames); 
        end
        set(gca,'TickLabelInterpreter','none');
        set(gca,'XTickLabelRotation',45);
        set(gca,'XTick',linspace(L(1),L(2),NumTicks));
        title(grouped(k).time(i).nameroot)
        ylabel('GFP/mCherry')
        if k==length(analyzedALL)
        legend(time,'location','best')
        end
        set(gca,'fontsize',20);
        set(gca,'linewidth',2);
        ylim([-1.2 1.5])
        end
%         axis([0 length(goodnames) -1 1.5])
    end
    analyzedALL(k).grouped=grouped(k);
end


figh=[];
for k=1:length(analyzedALL)
    timeassay{k}
figh(k)=figure(k);clf

    
if switches
mus=zeros(length(grouped(k).time),2);
sigs=zeros(length(grouped(k).time),2);
fold=[];
    if domedian
            for i=1:length(grouped(k).time)
        try
            mus(i,:)=grouped(k).time(i).goodmedF1mean;
            sigs(i,:)=grouped(k).time(i).goodmedF1sigma;
        catch
            grouped(k).time(i).goodmedF1mean(numel(mus(i,:)))=0;
            mus(i,:)=grouped(k).time(i).goodmedF1mean;

            grouped(k).time(i).goodmedF1sigma(numel(sigs(i,:)))=0;
            sigs(i,:)=grouped(i).goodmedF1sigma;
        end
        fold(end+1)=grouped(k).time(i).foldmedF1mean;
            end
    else
    for i=1:length(grouped(k).time)
        try
            mus(i,:)=grouped(k).time(i).goodmumean;
            sigs(i,:)=grouped(k).time(i).goodmusigma;
        catch
            grouped(k).time(i).goodmumean(numel(mus(i,:)))=0;
            mus(i,:)=grouped(k).time(i).goodmumean;

            grouped(k).time(i).goodmusigma(numel(sigs(i,:)))=0;
            sigs(i,:)=grouped(i).goodmusigma;
        end
        fold(end+1)=grouped(k).time(i).foldmean;
    end
    end
elseif timecourse
    mus=zeros(length(analyzedALL(1).grouped),length(analyzedALL));
    sigs=zeros(length(analyzedALL(1).grouped),length(analyzedALL));
    for k=1:length(analyzedALL)
        if domedian
        for i=1:length(analyzedALL(k).grouped)
        mus(i,k)=analyzedALL(k).grouped(i).goodmedF1mean;
        sigs(i,k)=analyzedALL(k).grouped(i).goodmedF1sigma;
        end
        else
        for i=1:length(analyzedALL(k).grouped)
        mus(i,k)=analyzedALL(k).grouped(i).goodmumean;
        sigs(i,k)=analyzedALL(k).grouped(i).goodmusigma;
        end
        end
        
    end
else
    if domedian
    for i=1:length(grouped(k).time)
        mus(i)=grouped(k).time(i).goodmedF1mean;
        sigs(i)=grouped(k).time(i).goodmedF1sigma;
    end
    else
    for i=1:length(grouped(k).time)
        mus(i)=grouped(k).time(i).goodmumean;
        sigs(i)=grouped(k).time(i).goodmusigma;
    end
    end

end



% sort lowest basal first

[Y,I]=sortrows(mus,1);
xlabelnames=uniquenames(I);
fold=fold(I);


m=zeros(size(mus));
mm=m;
s=zeros(size(sigs));

if normsTRSVctl
    sci=find(~cellfun('isempty',regexp(xlabelnames,'sTRSVctl')));
strsvctlmus=repmat(Y(sci,:),length(Y(:,1)),1);
Y=Y-strsvctlmus+2;

else
strsvctlmus=zeros(length(Y(:,1)),2);
Y=Y-strsvctlmus+2;

end

fold=Y(:,2)-Y(:,1);

for i=1:length(Y(1,:))
    in=~(isnan(Y(:,i)) & isnan(sigs(:,i)));
    m(in,i)=Y(in,i);
    s(in,i)=sigs(in,i);
    mm(in,i)=10.^Y(in,i);
%     xlabelnames=uniquenames(in);
end

% test for significance
Hs=[];
Ps=[];

meannorm_minus=[];
signorm_minus=[];
meannorm_plus=[];
signorm_plus=[];
meanfold=[];
sigfold=[];

numsamps=[];

for i=1:length(grouped(k).time)
    
    [H,P]=ttest2(grouped(k).time(i).noligmu-strsvctlmus(1,1),grouped(k).time(i).wligmu-strsvctlmus(1,2),'alpha',alpha);
    Hs(end+1)=H;
    Ps(end+1)=P;
    
    meannorm_minus(i)=mean(grouped(k).time(i).noligmu-strsvctlmus(1,1));
    signorm_minus(i)=sqrt(std(grouped(k).time(i).noligmu).^2+std(strsvctlmus(1,1)).^2);
    meannorm_plus(i)=mean(grouped(k).time(i).wligmu-strsvctlmus(1,2));
    signorm_plus(i)=std(grouped(k).time(i).wligmu-strsvctlmus(1,2));
    signorm_plus(i)=sqrt(std(grouped(k).time(i).wligmu).^2+std(strsvctlmus(1,2)).^2);
    meanfold(i)=mean(10.^((grouped(k).time(i).wligmu-mean(strsvctlmus(1,2)))-(grouped(k).time(i).noligmu-mean(strsvctlmus(1,1)))));
    sigfold(i)=std(10.^((grouped(k).time(i).wligmu-mean(strsvctlmus(1,2)))-(grouped(k).time(i).noligmu-mean(strsvctlmus(1,1)))));
    numsamps(i)=length(grouped(k).time(i).noligmu);

end

m=[meannorm_minus(I)' meannorm_plus(I)']+2;
figuredata(k).condition=expligandconc{k};
figuredata(k).meannorm_minus=meannorm_minus(I);
figuredata(k).signorm_minus=signorm_minus(I);
figuredata(k).meannorm_plus=meannorm_plus(I);
figuredata(k).signorm_plus=signorm_plus(I);
% figuredata(k).meanfold=10.^fold;
figuredata(k).meanfold=meanfold(I);
figuredata(k).sigfold=sqrt(figuredata(k).signorm_plus.^2+figuredata(k).signorm_minus.^2);
figuredata(k).samplenames=uniquenames(I);
figuredata(k).numsamps=numsamps(I);



numsamps=[numsamps(I)' numsamps(I)'];

if setYScaletoLog
    
    s=[signorm_minus(I)' signorm_plus(I)'];
    err=[];
    err(:,:,1)=10.^(m)-10.^(m-s);
    err(:,:,2)=10.^(m+s)-10.^(m);
else
    s=[signorm_minus(I)' signorm_plus(I)'];
    err=[];
    err(:,:,1)=(10.^(m)-10.^(m-s))./sqrt(numsamps);
    err(:,:,2)=(10.^(m+s)-10.^(m))./sqrt(numsamps);
end
    
mm=10.^m;


[h hErrorbar]=barwitherr(err,mm);


set(h,'linewidth',2)
set(h(1),'FaceColor',[0.5 0.5 0.5])
set(h(2),'FaceColor',[0.9 0.2 0.3])
h(1).BaseLine.LineWidth=2;

hval=Hs(I);
pval=Ps(I);



L = get(gca,'XLim');
% xlabelnames=uniquenames(in);
xlabelnames=uniquenames;
if length(xlabelnames)>10
set(gca,'xticklabel',{'',xlabelnames{I},'',''})
NumTicks = length(xlabelnames)+3; 
    
elseif  length(xlabelnames)>8
set(gca,'xticklabel',{'',xlabelnames{I},''})
NumTicks = length(xlabelnames)+2; 
    
elseif length(xlabelnames)>5
set(gca,'xticklabel',{'',xlabelnames{I},''})
NumTicks = length(xlabelnames)+2; 

elseif length(xlabelnames)==5 % for zeatin
set(gca,'xticklabel',{'',xlabelnames{I(1)},'',xlabelnames{I(2)},'',...
                      xlabelnames{I(3)},'',xlabelnames{I(4)},'',...
                      xlabelnames{I(5)},''})
NumTicks = length(xlabelnames)+6; 


elseif length(xlabelnames)==4 % for zeatin
set(gca,'xticklabel',{'',xlabelnames{I(1)},'',xlabelnames{I(2)},'',...
                      xlabelnames{I(3)},'',xlabelnames{I(4)},''})
NumTicks = length(xlabelnames)+5; 

else
set(gca,'xticklabel',{xlabelnames{I}})
NumTicks = length(xlabelnames); 
end
set(gca,'TickLabelInterpreter','none')
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'fontsize',22)
set(gca,'linewidth',2)
% set(gca,'YScale','log')

YL=get(gca,'Ylim');
if normsTRSVctl
    ylabel('% norm. mCherry/BFP')
else
    ylabel('\mu')
    set(gca,'YScale','log')

end
% title('CSY1171')

% legend('0 mM 6R,S-folinic acid', '5 mM 6R,S-folinic acid','location','best')
if switches
f=find(~isnan(fold));
time={};
for i=1:length(f)
    time{end+1}=sprintf('%0.1f',10.^fold(f(i)));
end

if normsTRSVctl && ~setYScaletoLog
txtheight=10.^(m(:,2)')+7.5;
axis([0 length(uniquenames)+1 0 max(YL)])
elseif setYScaletoLog
txtheight=10.^(m(:,2)')*1.5;
axis([0 length(uniquenames)+1 min(min(10.^m))*0.9 max(YL)*2.8])
else
txtheight=10.^(m(:,2)')*1.5;
axis([0 length(uniquenames)+1 min(min(10.^m))*0.9 max(YL)*1.5])

end

txtwidth=(1:length(uniquenames))-0.35;
text(txtwidth(f),txtheight(f),time,'fontsize',18)
hval(pval<(alpha/10) & pval>(alpha/100))=2;
hval(pval<(alpha/100) & pval>(alpha/1000))=3;
hval(pval<(alpha/1000))=4;

figuredata(k).hval=hval;
figuredata(k).pval=pval;


for j=1:length(hval)
    if normsTRSVctl && ~setYScaletoLog
        astoffset=3;
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
end

legendnames={sprintf('0 mM %s',ligandname),sprintf('%s %s',expligandconc{k},ligandname)};
legend(legendnames,'location','southeast')
% legend('0 mM theophylline', '1 mM theophylline','location','southeast')
% legend('0 mM aciclovir', '5 mM aciclovir','location','southeast')
% legend('0 \muM noscapine', '500 \muM noscapine','location','northwest')
if setYScaletoLog
    set(gca,'YScale','log')
end
ylabel(sprintf('Flow RFU\n(mCherry/BFP)'))
end

end