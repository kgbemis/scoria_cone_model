function plot_and_slice_diagnostics_cons_wide_sptest_func(stages,meta,fallstages,stage,v0)
% function plot_and_slice_volcano(matfilename)
%   loads output of volcano model then plots in 3D and 2D
%

%fprintf('starting diagnostics\n')

% pass variables instead of loading
%load(matfilename,'stages','meta','fallstages') 
%load('volcano_vel40_NppP5000_NP100_gs0.1_annulus0dir360std45dir5_rv3d_26-07-18_011742.mat')
% variables:    falltopo
%               meta
%               topo
%               stages

stagetopo=stages(:,:,stage);
[S,A]=local_gradm_simple(stagetopo);

figure(3)
sp1=subplot(221);
    contour(stagetopo,0:5:100)
    if v0==30
    axis([300 700 300 700 0 400])
    else
    axis([0 1000 0 1000])
    end
    caxis([0 75])
    colorbar
    daspect(sp1,[1 1 1])
    xlabel('Distance (meters)')
    ylabel('Distance (meters)')
sp2=subplot(222);
    pcolor(atand(S))
    %caxis([-45 45])
    caxis([0 45])
    colorbar
    shading flat
    title('slope')
    if v0==30
    axis([300 700 300 700 0 400])
    else
    axis([0 1000 0 1000])
    end
    daspect(sp2,[1 1 1])
sp3=subplot(223);
    pcolor(A)
    caxis([-180 180])
    colorbar
    shading flat
    title('aspect')
    if v0==30
    axis([300 700 300 700 0 400])
    else
    axis([0 1000 0 1000])
    end
    daspect(sp3,[1 1 1])
sp4=subplot(224);
    plot(stages(:,500,stage))
    xlabel('Distance (m)')
    ylabel('Height (m)')
    if v0==30
    axis([300 700 0 100])
    else
    axis([0 1000 0 100])
    end
    daspect(sp4,[2 1 1])
    vpar=meta.conditions;
    text(201,155,['Vel=' num2str(vpar.Vel) '; Nparticles=' num2str(vpar.Npp) ...
    '; Nphases=' num2str(stage)]);
    text(201,125,['Mav=' num2str(vpar.BurstAngle) ...
    '; lava vel=' num2str(vpar.LavaVel) '; Uo=' num2str(vpar.WindSpeed) ...
    '; KMU-' num2str(vpar.KMU)]);

outname=[matfilename(1:end-4) '_diagnostics_cons_wide.png'];
print(outname,'-dpng')

figure(4)
    contour(stagetopo,0:5:100)
    if v0==30
    axis([300 700 300 700 0 400])
    else
    axis([0 1000 0 1000])
    end
    caxis([0 75])
    colorbar
    daspect([1 1 1])
    xlabel('Distance (meters)')
    ylabel('Distance (meters)')
    grid on
outname4=[matfilename(1:end-4) '_con.png'];
print(outname4,'-dpng')
    
figure(5) % look at falltopo
    sumfall=sum(fallstages,3);
    maxfall=100*ceil(max(sumfall(:))/100); 
    contour(sumfall,0:5:maxfall)
    if v0==30
    axis([300 700 300 700 0 400])
    else
    axis([0 1000 0 1000])
    end
    caxis([0 maxfall])
    colorbar
    daspect([1 1 1])
    xlabel('Distance (meters)')
    ylabel('Distance (meters)')
    title('summed fall')
    %xlim([300 700])
    %ylim([300 700])
    grid on
outname5=[matfilename(1:end-4) '_summedfall.png'];
print(outname5,'-dpng')

if stage>100
figure(6) % look at falltopo
    sumfall=sum(fallstages(:,:,150:200),3);
    maxfall=100*ceil(max(sumfall(:))/100); 
    contour(sumfall,0:5:maxfall)
    if v0==30
    axis([300 700 300 700 0 400])
    else
    axis([0 1000 0 1000])
    end
    caxis([0 maxfall])
    colorbar
    daspect(sp1,[1 1 1])
    xlabel('Distance (meters)')
    ylabel('Distance (meters)')
    title('summed fall over fourth 50 phases')
    %xlim([300 700])
    %ylim([300 700])   
    
    %xlim([200 800])
    %ylim([200 800])
outname6=[matfilename(1:end-4) '_summedfall_4thquarter.png'];
print(outname6,'-dpng')
end

    
figure(7)
    pcolor(atand(S))
    caxis([0 50])
    colorbar
    shading flat
    title('slope')
    if v0==30
    axis([300 700 300 700 0 400])
    else
    axis([0 1000 0 1000])
    end
    daspect([1 1 1])
outname7=[matfilename(1:end-4) '_slope.png'];
print(outname7,'-dpng')

end