
analyzed1=analyzeflow('inputfilename','input1mM.txt','plotfit2D',1,'thresh2D',[2.7 2.4],'omit',true,'muhisto',true,'returndata',true);
analyzed2=analyzeflow('inputfilename','input5mM.txt','plotfit2D',1,'thresh2D',[2.7 2.4],'omit',true,'muhisto',true,'returndata',true);
%%
analyzedflowplot=analyzeflow('inputfilename','theoAAAAA_flowplots.txt','plotfit2D',1,'dimplot',[1,1],'thresh2D',[2.7 2.4],'omit',true,'muhisto',true);

%%

analyzed(1)=analyzed1;
analyzed(2)=analyzed2;
%%
time={'1 mM','5 mM'};
expligandconc={'1 mM','5 mM'};
timeassay=time;
analyzedALL=struct;
for k=1:length(expligandconc)
    analyzedALL(k).samplenames=analyzed(k).samplenames;
    analyzedALL(k).omit=analyzed(k).omit;
    analyzedALL(k).mu=analyzed(k).mu;
    analyzedALL(k).medFP1=analyzed(k).medFP2;
    analyzedALL(k).sigma=analyzed(k).sigma;
    analyzedALL(k).medsigFP1=analyzed(k).medsigFP2;
end


%%
timecourse=true;
switches=true;
plotindividualsamples=false;
domedian=false;
ligandname='theophylline';
normsTRSVctl=true;
setYScaletoLog=true;
alpha=1e-2;

figuredata=getVYBdata(analyzedALL,expligandconc,ligandname,time,alpha,normsTRSVctl,setYScaletoLog,domedian,plotindividualsamples,switches,timecourse,timeassay);

%%
swn='TheoAAAAA';
combinedsw={};
combinedsc={};
for j=1:2
    if j==1
        nameh={strcat('(',swn,')(.*)(-lig)')};
        scname={strcat('(sTRSVctl)(.*)(-lig)')};
    else
        nameh={strcat('(',swn,')(.*)(\+lig)')};
        scname={strcat('(sTRSVctl)(.*)(\+lig)')};

    end
s=find(~cellfun('isempty',regexp(analyzed1.samplenames,nameh)));
        sc=find(~cellfun('isempty',regexp(analyzed1.samplenames,scname)));

s_mat=[];
for i=1:length(s)
    s_mat=[s_mat;analyzed1.data(s(i)).GFPmCherryratio];
    fprintf('num cells = %d\n',length(s_mat))
end
combinedsw{j}=s_mat;

sc_mat=[];
for i=1:length(sc)
    sc_mat=[sc_mat;analyzed1.data(sc(i)).GFPmCherryratio];
    fprintf('num cells = %d\n',length(sc_mat))
end
combinedsc{j}=sc_mat;

end


%% kolmogorov-smirnoff test of normality
bw=0.05;


figure(1);clf
hold on
falpha=0.2;

kstestpval=[];
ksteststat=[];
for j=1:length(combinedsw)
    s_mat=combinedsw{j};
    sc_mat=combinedsc{j};
setfig('GM histo')
h=histogram(log10(s_mat)-log10(mean(sc_mat))+2,'binwidth',bw);
h.Normalization='probability';
figh{j}=h;
% end

figure(1)
hold on
colormat2={[0.6 0.6 0.6],[0.7 0.2 0.2],colormat{1},colormat{1},[0.8 0.8 0.8],[0.8 0.8 0.8],[0.8 0.8 0.8]};
% for j=1:length(nameh)
    if l==1 && j==1
        lstyle='--';
%         falpha=0.4;
    else
        lstyle='-';
%         falpha=0.4;
    end
    falpha=falpha+0.1;
area(figh{j}.BinEdges(2:end)-figh{j}.BinWidth/2,figh{j}.Values,...
     'linewidth',2,'edgecolor',colormat2{j}/2,'linestyle',lstyle,'facecolor',colormat2{j}*(1-falpha/2),'facealpha',falpha)
    
    s_mat_stand=log10(s_mat)-log10(mean(sc_mat));
    
    [h,p,kss]=kstest(s_mat_stand,'alpha',0.01);
    kstestpval(end+1)=p;
    ksteststat(end+1)=kss;
end



xlabel('mCherry/BFP')
ylabel('Frequency')
set(gca,'Fontsize',26)
set(gca,'linewidth',2)
legend('0 mM Theophylline','1 mM Theophylline')



%%  END











