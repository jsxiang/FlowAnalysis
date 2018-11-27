function [analyzed]=analyzeflow(varargin)
% inputfilename references a file that contains 2 columns, the first 
% specifying sample name and 2nd specifying location of the file
% plotfit2D is true or false to indicate whether the plots showing the
% fits should be displayed
% thresh2D is a vector specifying the threshold fluor levels demarcating
% transduced/transfected cells 

addpath ~/Documents/MATLAB/FACS
addpath ~/Documents/MATLAB/fca_readfcs/
addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/barwitherr
defaults=struct('inputfilename','','plotfit2D',1,'thresh2D',[],'dimplot',[],...
                'cytometer','','numcells',[],'transduced',0,'omit',false,...
                'muhisto',false,'returndata',false);
args=processargs(defaults,varargin);
if isempty(args.plotfit2D)
    args.plotfit2D=true;
end

if isempty(args.thresh2D)
    args.thresh2D=[2.05 1.95];
end

if isempty(args.cytometer)
    args.cytometer='VYB';
end

if isempty(args.numcells)
    args.numcells=10000;
end

if isempty(args.transduced)
    args.transduced=false;
end

if isempty(args.returndata)
    args.returndata=false;
end

if isempty(args.omit)
    args.omit=false;
end


if isempty(args.muhisto)
    args.muhisto=false;
end


fid=fopen(args.inputfilename);
samples=textscan(fid,'%s');
fclose(fid);

samplenames={};
channel1={};
channel2={};
filenames={};
omitsample={};

if args.omit
    int=5;
else
    int=4;
end
for i=1:int:length(samples{1})
    samplenames{end+1}=samples{1}{i};
    channel1{end+1}=samples{1}{i+1};
    channel2{end+1}=samples{1}{i+2};
    filenames{end+1}=samples{1}{i+3};
    if args.omit
    try
        omitsample{end+1}=samples{1}{i+4};
    catch
        omitsample{end+1}='';
    end
    end
end

analyzed.samplenames=samplenames;
analyzed.files=filenames;
analyzed.channel1=channel1;
analyzed.channel2=channel2;
analyzed.muratio=struct;
if args.omit
    analyzed.omit=omitsample;
end
%%
for i=1:length(analyzed.samplenames)
%     filename=strcat('data/',filenames{i});
    fcsdat = fca_readfcs(analyzed.files{i});
%     try 
%         fcsdat=dat(1:args.numcells,:);
%     catch
%         fcsdat=dat;
%     end
%     clearvars dat
    if strcmp(args.cytometer,'VYB')
        
%         fprintf('%s\n',analyzed.files{i})
%         fprintf('size %d\n',size(fcsdat))

        fsch=fcsdat(:,5);
        fsca=fcsdat(:,4);
        ssca=fcsdat(:,2);
        ssch=fcsdat(:,3);
        viablecells=fsca>10^3.5 & fsca<10^4.6 & ssca>10.^4.6 & ssca<10.^5.05;
        singlets=(fsch>(0.65*fsca))&viablecells;

        setfig('plot viable');clf
        hold on
        plot(log10(fsca),log10(ssca),'.','color',[0.6 0.6 0.6])
        dscatter(log10(fsca(viablecells)),log10(ssca(viablecells)),'bins',[100,100]);
        dscatter(log10(fsca(viablecells)),log10(ssca(viablecells)),'plottype','contour','bins',[100,100]);
        xlabel('log_1_0(FSC-A)')
        ylabel('log_1_0(SSC-A)')
%         set(gca,'Yscale','log')
%         set(gca,'XScale','log')
        axis([2.5 5.4 4.2 5.4])
        set(gca,'fontsize',24)
        set(gca,'linewidth',2)
        
        p=sprintf('Fraction of parent = %0.2f',sum(viablecells)/length(viablecells));
        text(3.5,4.4,p,'fontsize',20);
        p=sprintf('Count in gate = %d',sum(viablecells));
        text(3.5,4.3,p,'fontsize',20);
        box on
        
        
        
        setfig('plot singlets');clf
        plot(log10(fsca(viablecells)),log10(fsch(viablecells)),'.','color',[0.6 0.6 0.6])
        xlabel('log_1_0(FSC-A)')
        ylabel('log_1_0(FSC-H)')
%         viablecells=fsca>10^3 & fsca<10^4.2; %only for Y42
%         set(gca,'Yscale','log')
%         set(gca,'XScale','log')
        axis([2.7 5.4 2.7 5.4])
        singlets=(fsch>(0.45*fsca))&viablecells;  %only for Y42 and T86 cdG
        
        hold on
%         plot(fsca(singlets),fsch(singlets),'r.')
        dscatter(log10(fsca(singlets)),log10(fsch(singlets)),'bins',[100,100])
        dscatter(log10(fsca(singlets)),log10(fsch(singlets)),'plottype','contour','bins',[100,100])
        
        p=sprintf('Fraction of parent = %0.2f',sum(singlets)/sum(viablecells));
        text(3.5,3.2,p,'fontsize',20);
        p=sprintf('Count in gate = %d',sum(singlets));
        text(3.5,3.0,p,'fontsize',20);
        set(gca,'fontsize',24)
        set(gca,'linewidth',2)
        box on
        hold off

        ssc=fcsdat(singlets,2);
        fsc=fcsdat(singlets,4);
        
        h=fcsdat(singlets,6);
        v=fcsdat(singlets,8);
        GFP=fcsdat(singlets,10);
        
        % remove saturated
        nonsat=(h<5.4e5)&(v<5.4e5)&(GFP<5.4e5);
        
        data(i).name=analyzed.samplenames{i};
        data(i).ssc=ssc(nonsat);
        data(i).fsc=fsc(nonsat);
        data(i).BFP=h(nonsat);
        data(i).mCherry=v(nonsat);
        data(i).GFP=GFP(nonsat);
        data(i).singlets=singlets;
        data(i).singlet_nonsat=singlets(find(nonsat));
    elseif strcmp(args.cytometer,'Carmen')
        ssc=fcsdat(:,13);
        fsc=fcsdat(:,11);
        h=fcsdat(:,19);
        v=fcsdat(:,21);
        GFP=fcsdat(:,15);
        data(i).name=analyzed.samplenames{i};
        data(i).ssc=ssc;
        data(i).fsc=fsc;
        data(i).BFP=h;
        data(i).mCherry=v;
        data(i).GFP=GFP;
    end
        
end

%% Fit mu and generate 2D plot

analyzed.mu=[];
analyzed.transduced=[];
analyzed.silenceVert=[];
analyzed.silenceHorz=[];
analyzed.twocolor=[];

ind=1:length(samplenames);
threshHorz=args.thresh2D(1)*ones(size(ind));
threshVert=args.thresh2D(2)*ones(size(ind));
% if args.plotfit2D
%     setfig('twocolorplot');clf
% end
for i=ind;

    if strcmp(analyzed.channel1{i},'mCherry')
        horzchannel=abs(data(ind(i)).mCherry);
    elseif strcmp(analyzed.channel1{i},'clover')
        horzchannel=abs(data(ind(i)).GFP);
    elseif strcmp(analyzed.channel1{i},'GFP')
        horzchannel=abs(data(ind(i)).GFP);
    elseif strcmp(analyzed.channel1{i},'BFP')
        horzchannel=abs(data(ind(i)).BFP);
    else
        horzchannel=abs(data(ind(i)).fsc);
    end

    if strcmp(analyzed.channel2{i},'mCherry')
        vertchannel=abs(data(ind(i)).mCherry);
    elseif strcmp(analyzed.channel2{i},'clover')
        vertchannel=abs(data(ind(i)).GFP);
    elseif strcmp(analyzed.channel2{i},'GFP')
        vertchannel=abs(data(ind(i)).GFP);
    elseif strcmp(analyzed.channel2{i},'BFP')
        vertchannel=abs(data(ind(i)).BFP);
    else
        vertchannel=abs(data(ind(i)).ssc);
    end
        
    
    % estimate mu based on transduced cells, defined as double color
    tdataindex=log10(horzchannel)>threshHorz(i)&log10(vertchannel)>threshVert(i);
    tdvert=vertchannel(tdataindex);
    tdhorz=horzchannel(tdataindex);
    [m,s]=normfit(log10(tdvert./tdhorz));
    analyzed.mu(i)=m;
    analyzed.sigma(i)=s;
    
    % estimate mu based on transduced cells, defined as all BFP
    tdataindexALL=log10(horzchannel)>threshHorz(i);
    tdvertALL=vertchannel(tdataindexALL);
    tdhorzALL=horzchannel(tdataindexALL);
    [mALL,sALL]=normfit(log10(tdvertALL./tdhorzALL));
    analyzed.muALL(i)=mALL;
    analyzed.sigmaALL(i)=sALL;
    
    % find median of the horizontal channel, everything above the horz
    % threshold 
    indALLhorz=log10(horzchannel)>threshHorz(i);
    ALLhorz=horzchannel(indALLhorz);
    analyzed.medFP1(i)=median(log10(ALLhorz));
    analyzed.medsigFP1(i)=std(log10(ALLhorz));
    
    analyzed.medFP2(i)=median(log10(tdvert));
    analyzed.medsigFP2(i)=std(log10(tdvert));
    
    transduced=(1-sum(log10(horzchannel)<threshHorz(i)&log10(vertchannel)<threshVert(i))/length(vertchannel))*100;
    silencedVert=sum(log10(horzchannel)<threshHorz(i)&log10(vertchannel)>threshVert(i))/length(vertchannel)*100;
    silencedHorz=sum(log10(horzchannel)>threshHorz(i)&log10(vertchannel)<threshVert(i))/length(vertchannel)*100;
    twocolor=sum(log10(horzchannel)>threshHorz(i)&log10(vertchannel)>threshVert(i))/length(vertchannel)*100;

    analyzed.transduced(i)=transduced;
    analyzed.silencedVert(i)=silencedVert;
    analyzed.silencedHorz(i)=silencedHorz;
    analyzed.twocolor(i)=twocolor;
    
    if args.plotfit2D
        setfig(samplenames{i});clf
        if isempty(args.dimplot)
            dimensionploty=ceil(sqrt(ind(end)));
            if dimensionploty*(dimensionploty-1)>ind(end)
                dimensionplotx=dimensionploty-1;
            else
                dimensionplotx=dimensionploty;
            end
        else
            dimensionplotx=args.dimplot(1);
            dimensionploty=args.dimplot(2);
        end
%         setfig('twocolorplot');
%         subplot(dimensionplotx,dimensionploty,i)
        hold on
%         if length(horzchannel)>args.numcells
%             plotdatx=log10(horzchannel(1:args.numcells));
%         else
            plotdatx=log10(horzchannel);
%         end
%         if length(vertchannel)>args.numcells
%             plotdaty=log10(vertchannel(1:args.numcells));
%         else
            plotdaty=log10(vertchannel);
%         end
        
        hold on
        plot(plotdatx,plotdaty,'.','color',[0.6 0.6 0.6])
        dscatter(plotdatx(tdataindexALL),plotdaty(tdataindexALL),'marker','o','BINS',[100 100])
        dscatter(plotdatx(tdataindexALL),plotdaty(tdataindexALL),'plottype','contour','BINS',[100 100])
        
        x=args.thresh2D(1)*ones(1,100);
        y=linspace(0, 5.4);
        
        yp=args.thresh2D(2)*ones(1,100);
        xp=linspace(0,5.4);
        
        if args.transduced
            td=sprintf('Transduced = %2.1f%%',transduced);
            text(1.1,4.9,td)
            plot(x,y,'k--')
            plot(xp,yp,'k--')
        end

%         z=linspace(args.thresh2D(1),5.4)+analyzed.muALL(i);
        z=linspace(args.thresh2D(1),5.4)+analyzed.muALL(i);
        plot(linspace(args.thresh2D(1),5.4),z,'k:','LineWidth',2)
%         mutxt=sprintf('\\mu = %1.3f',analyzed.muALL(i));
        mutxt=sprintf('\\mu = %1.2f',analyzed.muALL(i));
%         sigmatxt=sprintf('\\sigma = %1.2f',analyzed.sigmaALL(i));
        sigmatxt=sprintf('\\sigma = %1.2f',analyzed.sigmaALL(i));
        text(4.55,2.2,mutxt,'fontsize',20);
        text(4.6,1.9,sigmatxt,'fontsize',20);
        
        cellcount=length(vertchannel);
%         c=sprintf('Count = %2.0f',cellcount);
        c=sprintf('Count in gate = %2.0f',sum(tdataindexALL));
        text(3.45,1.6,c,'fontsize',20);
        p=sprintf('Fraction of parent = %0.2f',sum(tdataindexALL)/sum(singlets));
        text(3.1,1.3,p,'fontsize',20);
        
        xname=sprintf('log_1_0(%s)',analyzed.channel1{i});
        yname=sprintf('log_1_0(%s)',analyzed.channel2{i});
        xlabel(xname)
        ylabel(yname)
        axis([1 5.5 1 5.5])
        set(gca,'FontSize',24)
        set(gca,'LineWidth',2)
        title(samplenames(ind(i)),'interpreter','none')
        box on
        hold off
    end
    
    %% plot histogram of mu

    if args.muhisto
    setfig('GFP/mCherry');
    data(i).GFPmCherryratio=(tdvert./tdhorz);
    
    subplot(dimensionplotx,dimensionploty,i)
    
    [n,c]=hist(log10(data(i).GFPmCherryratio),200);
    area(c,n,'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
    set(gca,'linewidth',2)
    set(gca,'fontsize',14)
    xlabel('\mu')
%     set(gca,'YScale','log')
    title(samplenames(ind(i)),'interpreter','none')
    ylabel('count')
    xlim([-2 3])
    
    setfig('overlay');
    data(i).GFPmCherryratio=(tdvert./tdhorz);
    
%     subplot(dimensionplotx,dimensionploty,i)
    
    h=histogram(log10(data(i).GFPmCherryratio));
    h.Normalization = 'probability';
    h.BinWidth = 0.05;
    hold on
%     area(c,n,'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
    set(gca,'linewidth',2)
    set(gca,'fontsize',14)
    xlabel('log10(\mu)')
%     set(gca,'YScale','log')
    title('\mu = GFP/mCherry')
    ylabel('frequency')
    xlim([-2 3])
    
    analyzed.muratio(i).muratio=data(i).GFPmCherryratio;
    legend(samplenames,'interpreter','none','location','best');
    
    end
    
    
%% return data
if args.returndata
    analyzed.data=data;
end

    
end














end
