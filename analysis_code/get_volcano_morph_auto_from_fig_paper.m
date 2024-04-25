function get_volcano_morph_auto_from_fig_paper(xx,yy,zz,stage,stepsize,filecode)
% function get_volcano_morph(infile, filecode)
%   input
%       zz = zdata
%       filecode for saving results
%   output
%       H = average height of volcano
%       W = width of volcano base
%       T = width of crater
%       D = depth of crater
%

filetitle=replace(filecode,'_',' ');

% get center-crossing profiles to use
[m,n,p]=size(zz);
if m~=n
    fprintf('Warning: Volcano Landing Pad not square.%dx%dx%d \n',m,n,p)
end
icenter=round(m/2);
fprintf('icenter %d \n',icenter)
tphase=1:stepsize:stage;
Nt=length(tphase);
zpro=zeros(m,2,Nt);
for i=1:Nt
    % get the two perpendiculars for this phase
    zpro(:,1,i)=zz(:,icenter,tphase(i));
    zpro(:,2,i)=zz(icenter,:,tphase(i))';
end
[mp,np,pp]=size(zpro);
Npro=pp;
prolen=mp;
properp=np;
fprintf('Extracted %f profiles of length %d in %d directions\n',Npro,prolen,properp)
procol=parula(Npro);

% determine if there is a nonzero preexisting topo
% check tilt between all four corners of stage==1 
x11=xx(1); y11=yy(icenter); 
x12=xx(m); y12=yy(icenter);
x21=xx(icenter); y21=yy(1); 
x22=xx(icenter); y22=yy(n);
tilt1=atan2d((zz(m,icenter,1)-zz(1,icenter,1)),sqrt((x12-x11)^2+(y12-y11)^2));
tilt2=atan2d((zz(icenter,n,1)-zz(icenter,1,1)),sqrt((x21-x22)^2+(y21-y22)^2));
fprintf('tilt1 %f  and tilt2 %f  \n',tilt1,tilt2)
pretopo(:,1)=tand(tilt1).*xx';
pretopo(:,2)=tand(tilt2).*yy;
fprintf('size pretopo [%d,%d] min %f  max %f\n',size(pretopo),min(pretopo(:)),max(pretopo(:)))

%
% set up variables that will store morph info
baseLx=zeros(Nt,2); baseRx=zeros(Nt,2);
baseLe=zeros(Nt,2); baseRe=zeros(Nt,2);
bmsLx=zeros(Nt,2); bmsRx=zeros(Nt,2);
bmsLe=zeros(Nt,2); bmsRe=zeros(Nt,2);
bmcLx=zeros(Nt,2); bmcRx=zeros(Nt,2);
bmcLe=zeros(Nt,2); bmcRe=zeros(Nt,2);
rimLx=zeros(Nt,2); rimRx=zeros(Nt,2);
rimLe=zeros(Nt,2); rimRe=zeros(Nt,2);
crbote=zeros(Nt,2); crbotx=zeros(Nt,2);
minedg=zeros(Nt,2); maxedg=zeros(Nt,2);
minedga=zeros(Nt,2); maxedga=zeros(Nt,2);
minedgc=zeros(Nt,2); maxedgc=zeros(Nt,2);
minedgs=zeros(Nt,2); maxedgs=zeros(Nt,2);
minedgsl=zeros(Nt,2); maxedgsl=zeros(Nt,2);
%
% loop over profiles
for i=1:Npro
    for j=1:2
    % get current surface for phase i, profile j and adjust for pre-existing topo
    cursurf=zpro(:,j,i);
    %fprintf('size cursurf [%d,%d]\n',size(cursurf))
    switch j
        case 1
            x=xx';
        case 2
            x=yy;
    end
    %
    % check the ends to see if buildup at edges
    leftminc=min(cursurf(1:100));  rightminc=min(cursurf(900:1000));
    minedgc(i,j)=min([leftminc rightminc]);
    leftmaxc=max(cursurf(1:100));  rightmaxc=max(cursurf(900:1000));
    maxedgc(i,j)=max([leftmaxc rightmaxc]);
    % compute spatially averaged topo (19 m spans; quadratic fit)
    %   tried moving average (using smooth with 19 m window) but misses
    %   peaks and smears topo too much
    smsurf=smooth(zpro(:,j,i),19,'sgolay',2);
    %fprintf('size smsurf [%d,%d]\n',size(smsurf))
    % check the ends to see if buildup at edges
    leftmin=min(smsurf(1:100));   rightmin=min(smsurf(900:1000));
    minedg(i,j)=min([leftmin rightmin]);
    leftmax=max(smsurf(1:100));  rightmax=max(smsurf(900:1000));
    maxedg(i,j)=max([leftmax rightmax]);
    %
    smsurfa=smsurf-pretopo(:,j);  %  adjusted for pre-existing topo
    %fprintf('size smsurfa [%d,%d]\n',size(smsurfa))
    %size(smsurfa)
    % check the ends to see if buildup at edges
    leftmina=min(smsurfa(1:100));   rightmina=min(smsurfa(900:1000));
    minedga(i,j)=min([leftmina rightmina]);
    leftmaxa=max(smsurfa(1:100));   rightmaxa=max(smsurfa(900:1000));
    maxedga(i,j)=max([leftmaxa rightmaxa]);
    %
    % compute slope for smoothed topo
    gradsurf=gradient(smsurf);
    slopsurf=abs(gradsurf);   % geographic slope is the magnitude of the gradient
    aspsurf=sign(gradsurf);   % in 3D this would be the direction of the gradient
    % check the ends to see if buildup at edges
    leftmins=min(slopsurf(10:100));
    rightmins=min(slopsurf(900:990));
    minedgs(i,j)=min([leftmins rightmins]);
    leftmaxs=max(slopsurf(10:100));
    rightmaxs=max(slopsurf(900:990));
    maxedgs(i,j)=max([leftmaxs rightmaxs]);
    leftminsl=min(slopsurf(1:100));
    rightminsl=min(slopsurf(900:1000));
    minedgsl(i,j)=min([leftminsl rightminsl]);
    leftmaxsl=max(slopsurf(1:100));
    rightmaxsl=max(slopsurf(900:1000));
    maxedgsl(i,j)=max([leftmaxsl rightmaxsl]);
    
    % compute surface curvature
    curvsurf=gradient(slopsurf);
    %procurv=curvsurf/(gradsurf.^2+1).^(3/2);
% estimate base edge locations
%   define cone edge as one of the following possibilities
%       - first or last place with elevation > 1 grain + original (usually 0)
%       - first or last place with slope > 0
%       - first or last place with |curvature| > 0
%basei=find(smsurfa>0.1); % 1 grain too small as buildup at edges
%basei=find(smsurfa>0.4); % this works with if statement
%if isempty(basei)  
%    basei=find(smsurfa>0.1); 
%end
basei=find(smsurfa(10:990)>0.4); % exclude outer edges of model (edge effects)
    % changing level critieria from 0.2 to 0.4 because of sloping pre-topo
    % affect that creates a sediment apron on the downslope side
                                 % 1 grain still outside cone
if isempty(basei)
    basei=find(smsurfa(10:990)>0.3); %relax criteria if nothing
end
if isempty(basei)
    basei=find(smsurfa(10:990)>0.2); %relax criteria further if nothing
end
basei=basei+9;
%fprintf('size basei %d %d\n',size(basei))
%maxsma=0.015*max(smsurfa); if maxsma<0.1,maxsma=0.1; end
%basei=find(smsurfa(10:990)>maxsma);  basei=basei+9;
%fprintf('1st basei %d and last basei %d\n',basei(1),basei(end))
baseLx(i,j)=x(basei(1)); baseLe(i,j)=smsurfa(basei(1)); 
baseRx(i,j)=x(basei(end)); baseRe(i,j)=smsurfa(basei(end));

%bmsi=find(slopsurf>0.08); % this works but is very tight to cone
bmsi=find(slopsurf(10:990)>0.06); 
if isempty(bmsi)
    bmsi=find(slopsurf(10:990)>0.03); 
end
bmsi=bmsi+9;
bmsLx(i,j)=x(bmsi(1)); bmsRx(i,j)=x(bmsi(end));
bmsLe(i,j)=smsurfa(bmsi(1)); bmsRe(i,j)=smsurfa(bmsi(end));
bmsLs=slopsurf(bmsi(1)); bmsRs=slopsurf(bmsi(end));

bmci=find(abs(curvsurf(10:990))>0.01); bmci=bmci+9;
bmcLx(i,j)=x(bmci(1)); bmcRx(i,j)=x(bmci(end));
bmcLe(i,j)=smsurfa(bmci(1)); bmcRe(i,j)=smsurfa(bmci(end));

% estimate crater rim locations
%   define crater rim as one of the following possibilities
%       - place < (or >) 1000 (assumed but could vary) with highest
%       elevation
%       - first or last place with slope dropping back to 0
%       - first or last place with |curvature| << 0
[rimLe(i,j),rimLi]=max(smsurfa(10:fix(m/2))); 
    rimLi=rimLi+10;
[rimRe(i,j),rimRi]=max(smsurfa(fix(m/2):end-10)); 
    rimRi=rimRi;  % does not need to add 10 as starts at center
%size(rimLi)
%size(rimRi)
%fprintf('phase %d: rimLi %d rimRi %d\n',i, rimLi,rimRi+fix(m/2))
%size(x)
rimLx(i,j)=x(rimLi); rimRx(i,j)=x(fix(m/2)+rimRi);

%rmsi=find(slopsurf>0.01); 
%rmsLx(i)=x(bmsi(1)); rmsRx(i)=x(bmsi(end));
%rmsLe(i)=smsurf(bmsi(1)); rmsRe(i)=smsurf(bmsi(end));

%rmci=find(curvsurf>0);
%rmcLx(i)=x(bmci(1)); rmcRx(i)=x(bmci(end));
%rmcLe(i)=smsurf(bmci(1)); rmcRe(i)=smsurf(bmci(end));

% estimate crater bottom location
[crbote(i,j),crboti]=min(smsurfa(rimLi:(fix(m/2)+rimRi)));
crbotx(i,j)=x(rimLi+crboti);
%display([crbotx(i) crbote(i)])

figure(4)
subplot(311)
    plot(cursurf,'-','Color',procol(i,:))  % actual topo(:,i)
    hold on
        plot(smsurf,'--','LineWidth',2,'Color',procol(i,:))
        plot(baseLx(i,j),baseLe(i,j),'o','Color',procol(i,:),'MarkerSize',10)
        plot(baseRx(i,j),baseRe(i,j),'o','Color',procol(i,:),'MarkerSize',10)
        plot(rimLx(i,j),rimLe(i,j),'^','Color',procol(i,:))
        plot(rimRx(i,j),rimRe(i,j),'^','Color',procol(i,:))
        plot(crbotx(i,j),crbote(i,j),'s','Color',procol(i,:))
    xlabel('Distance (m)')
    ylabel('Elevation (m)')
    %xlim([600 1400])
subplot(312)
    plot(slopsurf,'-','LineWidth',2,'Color',procol(i,:))
    hold on
        %plot(aspsurf,'x','LineWidth',1.5,'Color',procol(i,:))
        plot(zeros(length(zpro),1),':k')
        %plot(bmsLx(i,j),bmsLe(i,j),'o','Color',procol(i,:))
        plot(bmsLx(i,j),bmsLs,'o','Color',procol(i,:),'MarkerSize',10)
        %plot(bmsRx(i,j),bmsRe(i,j),'o','Color',procol(i,:))
        plot(bmsRx(i,j),bmsRs,'o','Color',procol(i,:),'MarkerSize',10)
    xlabel('Distance (m)')
    ylabel('Slope (-)')
    %xlim([600 1400])
subplot(313)
    plot(curvsurf,'-.','LineWidth',2,'Color',procol(i,:))
    hold on
       % plot(procurv,'--b','LineWidth',1.5)
        plot(bmcLx(i,j),bmcLe(i,j),'o','Color',procol(i,:))
        plot(bmcRx(i,j),bmcRe(i,j),'o','Color',procol(i,:))
    xlabel('Distance (m)')
    ylabel('Curvature (-)')
    %xlim([600 1400])       
    
figure(10)
subplot(111)
if i>1
    hold on
end
    plot(smsurfa,'-','LineWidth',2,'Color',procol(i,:))
    hold on
        plot(baseLx(i,j),baseLe(i,j),'o','Color',procol(i,:),'MarkerSize',10)
        plot(baseRx(i,j),baseRe(i,j),'o','Color',procol(i,:),'MarkerSize',10)
        plot(rimLx(i,j),rimLe(i,j),'^','Color',procol(i,:))
        plot(rimRx(i,j),rimRe(i,j),'^','Color',procol(i,:))
        plot(crbotx(i,j),crbote(i,j),'s','Color',procol(i,:))
    xlabel('Distance (m)')
    ylabel('Elevation (m)')
    %xlim([600 1400])
    
    end % loop over profiles
end % loop over phases
figure(4)
subplot(311)
hold off
subplot(312)
hold off
subplot(313)
hold off
legend('curvature','edge picks')

figure(10)
hold off
xlim([300 700])
title(filetitle)
daspect([1 1 1])
ax=gca;
ax.FontSize=12;
% after manual reshaping use
%print('code#_morph_profiles.png','-dpng','-r600')

% also save to a fig file so can do separate overlay
figname=[filecode '_morph_profiles.fig'];
disp(figname)
savefig(figname)

%% compute morphological parameters
% basal width
baseW=baseRx-baseLx;
bmsW=bmsRx-bmsLx;
bmcW=bmcRx-bmcLx;
% crater width
rimW=rimRx-rimLx;
% volcano height
htR=rimRe-baseRe; htL=rimLe-baseLe;
%htH=mean([htR htL],2);
htH=(htR+htL)/2;
% crater depth
rimAE=0.5*(rimRe+rimLe);
Dcr=rimAE-crbote;
%fprintf('BaseWidth %f %f \n',baseW)
%fprintf('Height %f %f\n',htH)
%fprintf('CraterWidth  %f %f\n',rimW)
%fprintf('CraterDepth %f %f\n\n',Dcr)

whos
% non-dim parameters
flatness=rimW./baseW;
steepness=2*htH./(baseW-rimW);
crsteep=2*Dcr./rimW;
rcd=Dcr./htH;
%fprintf('flatness steepness crsteep rcd\n')
%fprintf('%f %f %f %f\n\n',flatness,steepness,crsteep,rcd)

%% summary plot of morphological parameters
figure(5)  % change in basal width with growth phase
hw=plot(tphase,bmsW,'xb','MarkerSize',10);
hold on
    ht=plot(tphase,rimW,'>m','MarkerSize',10);
hold off
xlabel('phase')
ylabel('width (m)')
ax=gca;
ax.FontSize=18;
legend([hw(1),ht(1)],'Wco: dir 1 & 2','Wcr: dir 1 & 2',...
    'location','northwest','fontsize',14)
title(filetitle)
%print([filecode '_morph_fig5_width.png'],'-dpng','-r600')

figure(6)
hh=plot(tphase',htH,'pc','MarkerSize',10);
hold on
hd=plot(tphase',Dcr,'>m','MarkerSize',10);
hold off
xlabel('phase')
ylabel('height (m)')
ax=gca;
ax.FontSize=18;
legend([hh(1),hd(1)],'Hco: dir 1 & 2','Dcr: dir 1 & 2',...
    'Location','NorthWest','FontSize',14)
xlim([0 stage])
title(filetitle)
%print([filecode '_morph_fig6_height.png'],'-dpng','-r600')

figure(7)
hf=plot(tphase, flatness,'pc','MarkerSize',10);
hold on
hs=plot(tphase,steepness,'+g','MarkerSize',10);
hr=plot(tphase,rcd,'>b','MarkerSize',10);
hc=plot(tphase,crsteep,'om','MarkerSize',10);
hold off
xlabel('profile number')
ylabel('flatness or steepness (m)')
ax=gca;
ax.FontSize=18;
legend([hf(1),hs(1),hr(1),hc(1)],'flatness dir 1 & 2','steepness dir 1 & 2',...
    'rcd dir 1 & 2','crsteep dir 1 & 2','Location','SouthEast','fontsize',14)
xlim([0 stage])
ylim([0 1])
title(filetitle)
%print([filecode '_morph_fig7_shape.png'],'-dpng','-r600')

disp(hh)
disp(hd)

% save data to file
outname=['testing_larger_baselevelcriteria_morph_' filecode '.mat'];
save(outname,'tphase','baseW','rimW','htH','Dcr','flatness','steepness','crsteep','rcd')

return
% end of new auto code
