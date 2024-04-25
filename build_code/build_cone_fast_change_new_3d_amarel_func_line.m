function build_cone_fast_change_new_3d_amarel_func(v0,nphases,hc,kmu,key,thetaW,Uo)
% function build_cone_fast_change_new_3d_amarel_func(v0,hc,kmu)
%      input
%           v0 = initial velocity (m/s); default 30 m/s
%      nphases = number of phases to run; default 100
%           hc = critical slope or smu; default 0.6
%          kmu = kinematic coef of fric; default 0.4
%
fprintf('settings: v0=%d; nphases=%d; hc=%3.1f; kmu=%3.1f \n',v0,nphases,hc,kmu)
fprintf('key: %s \n',key)

%% 3D Volcanic Eruptions - No Drag, Iterating with time
% modified for 3D modelling !!!!  7-17-17 kgb
%clearvars
%close all
dograinflow=1;

% Parameters that change frequently:
if isempty(v0)
v0=30; % initial velocity in meters per second
       % assumed uniform for each burst (phase)
end
w=1000; % width of landing space in meters
N=100000; % number of particles per phase
if isempty(nphases)
nphases=100; % amount of phases
end
delstage=1;  % controls how often topo and falltopo saved
crw=100; % crater width for setting vent level
%sphases=1:100:nphases; % phase numbers at which to save or analyze topo
% specification of ejection directions
%   bursttype = annulus
%   assumes particle ejection has a uniform distribution in azimuth
%    azimuth = az = 360*rand(1,N);
%    angle from vertical = av = Rav*randn(1,N) + Mav
%   where Mav is the mean angle from vertical and Rav the standard
%   deviation about that angle in a normal distribution
%
%   bursttype = directed
%   assumes that burst has narrow range in both azimuth and angle from
%   vertical
%       az = Raz*randn(1,N) + Maz
%   where Maz is the mean azimuth and Raz is the std dev about the mean
%       av = Rav*randn(1,N) + Mav
%   where Mav is the mean angle from the vertical and Rav is the std dev about the mean
%
bursttype='annulus';  % either 'annulus' or 'directed'  *annulus
Maz=0; % controls direction (set to 0 for annulus)   *0
Raz=360; % controls range of directions (set to 360 for annulus)   *360
Mav=55; % controls angle from vertical (value must be btwn 0 and 90)   *15
        % actual value IS angle from horizontal
        % max distance traveled for given velocity at Mav=45deg
Rav=10; % controls range of angle from vertical (value must be btwn 0 and 90)    *5
change='none';  % eventually will implement changing direction for directed burst  (set to 'none' for annulus)

%Constants
gftime=240; % grain flow time (std = 120)
g=9.8; % modified gravity
gs=0.1; % grain size in meters
r=w/2; % horizontal distance (radius from the cone)
n=0; % n=number of iterations in loop computing particle paths
%binsize=1; % 1 meter bin sizes
%zatx=zeros(w,w); % 1 meter bin values with length = width (initial topogrpahy surface)
    %zatx=0.000001*(1:w)'*(1:w); % 0.0001 give 100m elevation diff across
				% try 0.000001 so is slope but not important
    %  linind7 testing if angle of original topo matters -YES, with exactly
    %  zero slope cone is perfect and symmetrical (note internal label may
    %  be linind6 as forgot to update filename in code)
    % linind8 try a different slope direction
    %zatx=0.000001*(1:w)'*(1:w);
    %zatx=rot90(zatx);
    % linind9 try 180 from original
    %zatx=0.000001*(1:w)'*(1:w);
    %zatx=rot90(zatx,2);
    % linind10 try 270 from original with even shallower slope
    %zatx=0.00000001*(1:w)'*(1:w);
    %zatx=rot90(zatx,3);
    % linind11 back to start but saving slopes so can make direction random
    % when slope small
    zatx=0.000001*(1:w)'*(1:w);
x0c=r; % crater is in the middle
x0range=0; % separation of vents 
y0=r;
z0=zatx(round(r),round(r)); % initial z position is zero
%z=0; % starting z position
if isempty(hc)
hc=0.6; % minimum criteria for grain flow (~0.6 is classic)
end
if isempty(kmu)
kmu=0.4; % coefficient of kinetic friction  (~0.4 most runs)
end
%smu=0.62; % coefficient of static frition (this is really same as hc)
% Notes: must have kmu<=hc or not physically reasonable
%   coefftest 1 hc=1.0 kmu=1.0
%   coefftest 2 hc=1.0 kmu=0.8
%   coefftest 3 hc=1.0 kmu=0.6
%   coefftest 4 hc=1.0 kmu=0.4
%
dtb=0.05; % time step for ballistic computations (seconds)
dtg=2; % time step for grain flow computations (seconds?)

% lava flow parameters
vL=0; % lava flow velocity at vent
dL=crw; % lava flow width at vent set to crater width
dtl=1; % lava flow time scale (seconds)
thetaL=0; % lava flow direction (initally setting to positive x-axis direction)

% wind parameters
if exist('thetaW')
	fprintf('running with wind direction = %f degrees',thetaW)
else
	thetaW=90; % wind direction
	fprintf('using default wind direction = %f degrees',Uo)
end
if exist('Uo')
	fprintf('runnig with wind speed = %f m/s',Uo)
else
	Uo=0; % wind speed
	fprintf('using default wind speed = %d m/s',Uo)
end

% put all the above into a meta data structure for later storage
meta.conditions.Nphases=nphases;
%meta.conditions.Sphases=sphases;
meta.conditions.Npp=N;
meta.conditions.Vel=v0;
meta.conditions.BurstType=bursttype;
meta.conditions.BurstDir=Maz;
meta.conditions.BurstDirStdDev=Raz;
meta.conditions.BurstAngle=Mav;
meta.conditions.BurstAngleStdDev=Rav;
meta.conditions.BurstChange=change;
meta.conditions.GrainSize=gs;
meta.conditions.KMU=kmu;
meta.conditions.SMU=hc;
meta.conditions.LavaVel=vL;
meta.conditions.LavaWd=dL;
meta.conditions.LavaTime=dtl;
meta.conditions.LavaDir=thetaL;
meta.conditions.WindSpeed=Uo;
meta.conditions.WindDir=thetaW;
meta.conditions.GFtime=gftime;
meta.conditions.TimeStepBallistic=dtb;
meta.conditions.TimeStepGrainflow=dtg;
meta.conditions.VentWidth=crw;

% set file name to save data and figure
volfilename=[key '_volcano_vel' num2str(v0) '_NppP' num2str(N) ...
    '_NP' num2str(nphases) '_gs' num2str(gs) '_' bursttype ...
    num2str(Maz) 'dir' num2str(Raz) 'std' num2str(Mav) 'dir' ...
    num2str(Rav) '_vL' num2str(vL) '_rv3d_' datestr(now,'dd-mm-yy_HHMMSS')];
meta.runinfo.timestamp=now;
meta.runinfo.version='rv3d_test';

% this is a drastic rewrite by Karen bc Rajvi's version took an hour to
% complete 10-20 phases .... way too long when we need 500 phases
% -- rewrite for 3D further modifies the code (kgb, 7-18-17)

% setup topo array to save results in
%topo=zeros(w,w,1);  % right now only saving final results
falltopo=zeros(w,w,1); % right now only saving final results
%allpts=1:binsize:w;
%ptside=ones(size(allpts));
%ptside(allpts<r)=-1;  % normally set to -1 

% setup lava flow locations
iLx=round(w/2):w;  % binsize assumed to be 1
iLy=round((w-dL)/2):round((w+dL)/2);

% Loop over Phases for Successive Layers:
stagenum=0;
stages=zeros(1000,1000,nphases);
fallstages=zeros(1000,1000,nphases);
numtoosteep=zeros(nphases,1);
keepxland=zeros(N,nphases);
keepyland=zeros(N,nphases);
keeplandtime=zeros(N,nphases);
keepmeanxflow=zeros(nphases,gftime/dtg);
keepmeanyflow=zeros(nphases,gftime/dtg);
keepgpoints=cell(nphases,gftime/dtg);
keepaz=zeros(N,nphases);
keepav=zeros(N,nphases);
keepS=zeros(1000,1000,nphases);
keepA=zeros(1000,1000,nphases);
keephflow=zeros(1000,1000,nphases);
savx0=zeros(nphases);
for phases=1:nphases
    fprintf('starting phase %d\n',phases)    
    
    % Ballistics stage -- find landing places based on ballistic paths
    %   this part did not seem all that slow so will leave alone for now
    % SETUP Ejection Angles
    switch bursttype
        case 'annulus'
            az = Maz+Raz*rand(1,N);  % note that this case uses uniform distribution here
            av = Mav+Rav*randn(1,N);
        case 'directed'
            az = Maz+Raz*randn(1,N);
            av = Mav+Rav*randn(1,N);
    end
    keepaz(:,phases)=az;
    keepav(:,phases)=av;
    
    v0x=v0*cosd(av).*cosd(az)+Uo*cosd(thetaW); %initial x-velocity (adj by "wind")
    v0y=v0*cosd(av).*sind(az)+Uo*sind(thetaW); %initial y-velocity (adj by "wind")
    v0z=v0*sind(av); %initial z-velocity

    x0=x0c - 3*x0range + 2*x0range*randi(2); % provides one or two distinct locations per phase
    savx0(phases)=x0;

    z=@(t) z0+v0z.*t-0.5*g.*t.^2; %vertical particle path
    x=@(t) x0+v0x.*t; % horizontal particle path
    y=@(t) y0+v0y.*t; % horizontal particle path
    vz=@(t) v0z-g*t; % vertical velocity along path
    
    zpath=zeros(5/dtb,N); %makes empty matrix as placeholder
    xpath=zeros(5/dtb,N); %makes empty matrix as placeholder
    ypath=zeros(5/dtb,N); %makes empty matrix as placeholder
    vzpath=zeros(5/dtb,N); %makes empty matrix as placeholder
    for t=0:dtb:30 %looping with time with 0.5 second increments
        n=n+1; %counts number of loops
        zpath(n,:)=z(t);
        xpath(n,:)=x(t);
        ypath(n,:)=y(t);
        vzpath(n,:)=vz(t);
    end
    
    % Finds the time each particle hits the ground
    % modified version to find landing spot -- kgb 7-19-17
    %zindex = zeros(size(zpath(1,:))); % orig version
    zindex = zeros(1,N);  % modified so #rows = 1; #col = #particles directly
    for p=1:N
        % first - get z values of topo along xpath for given P
        %xpathindex=round(xpath(:,p)); 
        xpathindex=fix(xpath(:,p))+1; % 7/27/22 changing to fix()+1 to control bin better
            xpathindex(xpathindex>w)=w;
            xpathindex(xpathindex<1)=1;
        %ypathindex=round(ypath(:,p)); 
        ypathindex=fix(ypath(:,p))+1; % same here
            ypathindex(ypathindex>w)=w;
            ypathindex(ypathindex<1)=1;
        zatxpath=zeros(size(zpath(:,p)));
        for in=1:n
        %zatxpath(in)=zatx(xpathindex(in),ypathindex(in)); % original
        % flipping row & col (lintest5)
        zatxpath(in)=zatx(ypathindex(in),xpathindex(in));
        end
        % then - look for where zpath drops below zatxapath
        %zindex(p)=find((zpath(:,p)-zatxpath)<0, 1 );  % original
        %tempzindex=find((zpath(:,p)-zatxpath)<0, 2);  % should output 2 values
        %if sign(vzpath(tempzindex(1)))==-1
        %    zindex(p)=tempzindex(1);  %  picks 1st intersection if going down
        %else
        %    zindex(p)=tempzindex(2);  % picks 2rd intersection if going up
        %end
        tempzindex=find((zpath(:,p)-zatxpath)<0, 3 );  % should output 3 values
        if sign(vzpath(tempzindex(1)))==-1
            zindex(p)=tempzindex(1);  %  picks 1st intersection if going down
        elseif sign(vzpath(tempzindex(2)))==-1
            zindex(p)=tempzindex(2);  %  picks 2nd intersection if going down
        else
            zindex(p)=tempzindex(3);  % picks 3rd intersection if going up
        end
    end
    tg=zindex*dtb; %time particle hits the ground. Time before the y-distance becomes negative
    xland=x0+v0x.*tg; %finds the horizontal landing distance from crater
    yland=y0+v0y.*tg; %finds the horizontal landing distance from crater
    xlandbin=round(xland); %+1 %integer number for horizontal landing distance. Index in bin
    ylandbin=round(yland); %+1 %integer number for horizontal landing distance. Index in bin
    % check for out of bounds
    xlandbin(xlandbin<1)=1;
    xlandbin(xlandbin>w)=w;
    ylandbin(xlandbin<1)=1;
    ylandbin(xlandbin>w)=w;   

%     figure
%     topo=histogram(a,w); %histogram of angles (for last layer...)

    %landings=[xlandbin' ylandbin'];
    % flipping x & y in estimating zatxpath didn't change things so try
    % here (linind6)
    landings=[ylandbin' xlandbin'];
    [unique_val,~,ic]=unique(landings,'rows');
    %[occurrences]=hist(ic,unique(ic)); %finds any repitition of landing location
    occurrences = accumarray(ic,1); % should do what the above line did without deprecated code
    uni_ind=sub2ind(size(zatx),unique_val(:,1),unique_val(:,2));
    zatx(uni_ind)=zatx(uni_ind)+(occurrences.*gs); %counts particles that fall into same bin
    n=0;

    falltopo(uni_ind)=occurrences.*gs;
    fprintf('... finished ballistics, starting grainflow \n')  
    fprintf('       mean landing time %f s\n',mean(tg(:)))

    if dograinflow
    %Part 2 - Grain Flow
    
    % get right and left differences in height
    % ARGH -- actually need to work with gradient instead for 3D --- ARGH
    %rightdiff=zatx(1:end-1)-zatx(2:end);  % equiv to x(i)-x(i+1)
    %temp=zatx(2:end)-zatx(1:end-1);  % equiv to x(i)-x(i-1)     
    %leftdiff=[0 temp(1:end-1)]; % shift over so match up 
    
    % estimate slope and aspect (direction of slope)
    [S,A]=local_gradm_simple(zatx);
    %absS=abs(S); % nope - don't want to use abs(S) as negative S implies
    %ascent not descent
    keepS(:,:,phases)=S;
    keepA(:,:,phases)=A;

    % get indices of steep places
    gpoints=find(S>=hc);
        if ~isempty(gpoints)
            fprintf('... %d places too steep\n',length(gpoints))
        else
            fprintf('... no places too steep\n')
        end

    %pause  
    t2=0;
    while ~isempty(gpoints)
            % zero out slopes too small for grainflow
            S(S<hc)=0;  % note that the "slope" here is dL/dh rather than an angle

            % set direction of flow
            % this is in matrix A as the aspect

            % precompute most values
            %   change S=dL/dh to aflow=angle of slope
            aflow=atand(S);
            %   set default movement in index units to zero
            indflow=zeros(w,w);
            %   compute total horizontal movement for one time unit
            %       added factor of 0.5  (kgb 072518)
            xflow=0.5*(dtg.^2)*g*(sind(aflow)-kmu*cosd(aflow));  % in meters!
            % test whether larger xflow helps - just increases spread
            %xflow=1.0*(dtg.^2)*g*(sind(aflow)-kmu*cosd(aflow));  % in meters!
                keephflow(:,:,phases)=xflow;
            %   break into x-y components
            %   should have a conversion from meters to indicies but right
            %   now binsize is 1 m so they are numerically the same!
            x_xflow=round(xflow.*cosd(A));
            y_xflow=round(xflow.*sind(A));
                % save values for diagnosis
                keepmeanxflow(phases,round(t2/dtg)+1)=mean(x_xflow(:));
                keepmeanyflow(phases,round(t2/dtg)+1)=mean(y_xflow(:));
            % computer number of bins to move
            %indflow(S>=hc)=x_xflow(S>=hc)*w+y_xflow(S>=hc); % original
            % kgb 8/7/22 modified to align with correct Matlab linear
            % indexing - starting with keeping assumuption that x tracks
            % columns and y rows (lintest1 - weird result)
            %indflow(S>=hc)=w*(x_xflow(S>=hc)-1)+y_xflow(S>=hc);
            % try opposite - nope, original was most likely correct
            % (lintest2)
            %indflow(S>=hc)=x_xflow(S>=hc)+w*(y_xflow(S>=hc)-1);
            % this flips the original - result looks okay (lintest3)
            %indflow(S>=hc)=x_xflow(S>=hc)+y_xflow(S>=hc)*w;
            % testing original again (lintest4)
            indflow(S>=hc)=x_xflow(S>=hc)*w+y_xflow(S>=hc); % original
            % for lintest 5 will flip row,col in landing path estimate on
            %     line ~188
            % lintest 12 - changed no grain flow to use absolute value of
            %     slope rather than actual slope (S can be large negative)
            %   ARGH! stopped at mid stage 9 as memory issue potentially - the 
            %     temporarty result is nasty - gpoints was geting worse not better
            %     horizontal flow both pos and neg
            % lintest 13 - additionally adjusted aflow to use abs(S) in
            %     hopes that makes horizontal flow always positive
            %       better but still goofy
            %     checked algorithm -- nope, abs(S) is wrong because neg S
            %     implies ascent not descent
            % lintest 14 - back to original!!!!! will also modify
            %     local_gradm_simple to randomize directions for very small
            %     slopes which looks like S<0.0001 or something (this
            %     really only randomized when multiple directions were
            %     steepest)
            % lintest 15 - sticking with initial eqn but adjusting slope,
            %     aspect to ignore steepest line of desecent when S<0.001
            %     and instead pick random direction - will do just 30
            %     phases to see if matters - doesn't look like it does
            % lintest 16 - going to double xflow and see if that helps -
            %     simply changed 0.5 to 1.0 at beginning of eqn - need to
            %     dig deeper into this and think about timing as values of
            %     xflow often below 1 which means no movement - running
            %     overnight so doing 100 phases - this is artifical test
            %     but want to see how much matters - need time to go
            %     through equations --> everything farther out but basic
            %     shape is the same
            % lintest 17 - will revert to 0.5 - will also revert to not
            %     fully random even at low slope - will look at result of
            %     changing kmu to be less than smu (= hc) ==> eek! it's
            %     square and xflow gooofy
            % lintest 18 - changing kmu to >smu - this looks better - hmm,
            %     effects are really sensitive to values and ratio
            % lintest 19 - setting kmu=0.4 smu=0.6 which is where most
            %     models seem to have been run - see if that still looks
            %     okay

            % check number of moving points
            %Nmove=sum(xflow>0);
            %if Nmove > length(gpoints)
            %    fprintf('too many moving\n')
            %end

            % shift points from old to new locations
            %       this probably still assumes that never hit edges (keep w big)
            % remove particle from old spot
            zatx(gpoints)=zatx(gpoints)-gs; 
            % check for out of index points
                % oops this puts out of index points at the edge
                % might be better to let them diappear! Fix this. 12-30-20
            gtemp=gpoints+indflow(gpoints);
            gtemp(gtemp>w*w)=w*w;
            gtemp(gtemp<1)=1;
            % add particle to new spot
            zatx(gtemp)=zatx(gtemp)+gs.*ones(length(gpoints),1); %adds new particle location to topography array
        
             %pause
             t2=t2+dtg;
             %fprintf('time %4.1f \n',t2)
             if(t2>gftime)
                 fprintf('t2 hit gftime \n')
                 break
             end
        % start next loop by seeing if still steep places
        % get right and left differences in height
        %rightdiff=x(1:end-1)-x(2:end);  % equiv to x(i)-x(i+1)
        %temp=x(2:end)-x(1:end-1);  % equiv to x(i)-x(i-1) 
        %leftdiff=[0 temp(1:end-1)];
        % estimate slope and aspect (direction of slope)
        [S,A]=local_gradm_simple(zatx);
        %absS=abs(S); - nope, this is wrong - neg S is ascent not descent
        
        % get indices of steep places
        gpoints=find(S>=hc);
        if ~isempty(gpoints)
            fprintf('... %d places too steep\n',length(gpoints))
        else
            fprintf('... no places too steep\n')
        end
        keepgpoints{phases,round(t2/dtg)+1}=gpoints;
        %fprintf('now %d steep points \n',length(gpoints))
        %fprintf('now %d topo < 0 \n',sum(x<0))
        %fprintf('before %4.1f %4.1f \n',x(gpoints))

        %display(x(gpoints))
        %pause
    end
    fprintf('t2 is %4.1f\n',t2)
    %fprintf('at end of phase %d ... %d places too steep\n',phases,length(gpoints))
    numtoosteep(phases)=length(gpoints);
    end
    % lava flow
    % for now 
    %   all lava movement left to end of phase; probably should
    %       implement to have (a) lava movement after ballistics and  
    %       (b)lava movement during grainflow loops
    %   lava moves in postive x direction only
    %   lava movement intiates at crater
    if vL~=0  % fiat no effect if zero velocity
    dxL=round(vL*(mean(tg(:))+t2)); % distance moved in meters (and bins)
    newtopo=zatx(iLx,iLy);  % initialize to current state
    for ixl=1:length(iLx)
        for iyl=1:length(iLy)
            if iLx(ixl)==round(w/2) 
                % at center, keep same elevation
                newtopo(ixl,iyl)=zatx(iLx(ixl),iLy(iyl));
            elseif iLx(ixl)<= (round(w/2)+crw/2)
                % in crater, keep same elevation
                newtopo(ixl,iyl)=zatx(iLx(ixl),iLy(iyl));
            else
                % new height should be same as that dxL towards crater
                newtopo(ixl,iyl)=zatx(iLx(ixl)-dxL,iLy(iyl));                
            end
        end
    end
    zatx(iLx,iLy)=newtopo;
    end  % end of lava section
    
    % save topo for later plotting
    if rem(phases,delstage)==0
        stagenum=stagenum+1;
        fprintf('... saving phase %d as stage %d \n',stagenum,phases)
        stages(:,:,stagenum)=zatx;
        fallstages(:,:,stagenum)=falltopo;
        keepxland(:,stagenum)=xland';
        keepyland(:,stagenum)=yland';
        keeplandtime(:,stagenum)=tg;
    end
    % adjust vent level - set to mean center height
    ixcen=round(w./2); iycen=round(w./2);
    zcrater=zatx(ixcen-crw:ixcen+crw,iycen-crw:iycen+crw);
    z0=mean(zcrater(:));
    % reset falltopo
    falltopo=zeros(w,w,1);
end
topo(:,:,1)=zatx;

%% Graph of volcano
fig1=figure(1);
    %plot(topo)
    surf(topo)
    shading flat
    %hold on
    %xlabel('Topography Width (meters)')
    %ylabel('Height (meters)')
    %axis([0 w 0 w/2])
    %axis equal
    %axis([0 w 0 w])
    %axis auto
    %axis image
    daspect([1 1 1])
    if v0==30
        axis([300 700 300 700 0 400])
    end
    
%% save
%save([volfilename '.mat'], 'topo', 'falltopo', 'meta','stages','fallstages','numtoosteep',...
%    'keepxland','keepyland','keepmeanxflow','keepmeanyflow','keepgpoints',...
%    'keepaz','keepav','-v7')
save([volfilename '.mat'], 'topo', 'falltopo', 'meta','stages','fallstages','numtoosteep',...
    'keepxland','keepyland','keepmeanxflow','keepmeanyflow','keepgpoints',...
    'keepaz','keepav','keeplandtime','-v7.3')
fprintf('main data file saved \n')
%save([volfilename '_slopetracking.mat'], 'keepS', 'keepA','keephflow','-v7.3')
%printf('slope tracking data file saved \n')
savefig(fig1,[volfilename '.fig'])
print([volfilename '.png'],'-dpng')
fprintf('volcano figure saved \n')

% shape diagnostics 
fprintf('skipping diagnostics \n')
%plot_and_slice_diagnostics_cons_wide([volfilename '.mat'],nphases)
%   8/29/22 adjusted next command to save this and other figures
%plot_and_slice_diagnostics_cons_wide_sptest([volfilename '.mat'],nphases,v0)
%fprintf('diagnostics complete \n')

end
