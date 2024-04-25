Readme for scoria cone model
(NOTE: most up-to-date information will be in the readme_scoria_cone_model.m file)

this file describes the code for running the scoria cone model
      
note that code assumes Matlab R2020a version but has been run in R2023b with no obvious issues

The two folders:
1. build_cone contains the main simulation code:
   build_cone_fast_change_new_3d_amarel_func_line.m   - this is the main function
   local_gradm_simple.m                               - this is called by the main function
   plot_and_slice_diagnostics_cons_wide_sptest.m      - this could be called by the main function (but is currently commented out)
   runfunc.matlab.build				   - this is an example of a batch script for a HPC scheduler
2. analysis_tools contains some key functions for looking at the results:
   plot_and_slice_diagnostics_cons_wide_sptest_func.m  
   local_gradm_simple.m                                
   get_volcano_morph_auto_from_fig_paper.m             
   run_get_morph_paper.m

    
command notes for build code
basic command:
   build_cone_fast_change_new_3d_amarel_func(v0,nphases,hc,kmu,key,thetaW,Uo)
   
input:
      v0 = ejection velocity for run (tested values 30, 40, 50)
      nphases = number of phases to run
      hc = coefficient of static friction (smu)
      kmu = coefficient of kinetic friction
      key = name to identify run in output filenames
      thetaW = wind direction (note that drag force is not implemented; this
          value has no effect when Uo=0)
      Uo = wind speed (note that drag force is not implemented so this
          should be set to 0; results are invalid otherwise)
          
output:
  most output saved to files whose name includes a timestamp to prevent
	  ever overwriting past runs
	final topography, stages of topography, fall deposits, and various 
	  other values are saved in a .mat file at the end of a run
  a stream of diagnostic comments are output to either stdout or the 
    Matlab interactive command screen
  one 3D figure of final topography is created and saved to both .fig 
		and .png formats

hard coded parameters that can be changed
      w = width of the simulation space
      N = number of particles per phase 
      crw = vent width (mostly controls lava flow width; also effects 
          region reset to average vent level at end of each phase)
      bursttype = shape of particle ejection (usually set to "annulus"
          and only this option really tested)
      Mav,Rav = mean and standard deviation of angle from the horizontal
          at which particles are ejected
          (typically use 55 and 5 respectively but this can be varied)
     Maz,Raz = mean and standard deviation of azimuth at which particles
          are ejected (usually set to 0 and 360 respectively which
          produces an annulus of ejection)
     gftime = maximum in model time of interations within a grainflow
          phase (generally 240 s)
     gs = grain size (only controls the height of the landed blocks and
          lava flow for now)
     zatx = intial topography
           currently default is zatx=0.000001*(1:w)'*(1:w);
     x0range = separation of vents for line experiments (if set to 0
          only one vent is used; otherwise 2 vents separated by
          2*x0range are used)
     dtb,dtg,dtl= time steps for ballistics, grain flow and lava flow
          respectively
     vL = lava flow speed
      thetaL = lava flow direction


the command below will run from a script or the command line or from a batch script 
sent to HPC scheuduler (an example such script is included) 

build_cone_fast_change_new_3d_amarel_func_line(30,5,0.6,0.4,'onevent',90,0)

some notes on the analysis_tools 
 a. plot_and_slice_diagnostics_cons_wide_sptest_func.m creates several plots from the .mat
    file output including a 3D view, a diagnostics page with contours, slope, aspect, 
    and profile, and a  seperate contour plot.  It calls local_gradm_simple.m to do this.                           
 b. run_get_morph_paper.m	calls get_volcano_morph_auto_from_fig_paper.m to extract
    morphometric parameters from the output of the build code.     	
