
analyzed1=analyzeflow('inputfilename','input1mM.txt','plotfit2D',0,'thresh2D',[2.7 2.4],'omit',true);
analyzed2=analyzeflow('inputfilename','input5mM.txt','plotfit2D',0,'thresh2D',[2.7 2.4],'omit',true);

%%%%
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





%%  END











