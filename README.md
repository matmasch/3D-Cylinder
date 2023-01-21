# 3D-Cylinder
This matlab code will synthesize a CNT forest around the perimeter of a cylinde, representing a microfiber.

The main CNT forest synthesis function is ThreeD_Cylinder.m. You will want to modify the fname variable to match an output directory on your system.

If you do not already have fsparse installed in your working directory, run make.m FIRST.


clear
E=1e12; rout=5e-9; ri=rout*0.85; A=pi*(rout^2-ri^2); Ic=pi/2*(rout^4-ri^4); EA=E*A; EI=E*Ic; G=5e11; 
title=input('Name of Run ','s');

fname = 'C:\Users\maschmannm\Documents\MATLAB\Image Folder\Junk\';


steps=200;% Number of time steps
diam=8;        h_span=diam*1e-6;    %      ate lengthwise(micron)
length=10;      v_span=length*1e-6   % Span of the substrate breadthwise(micron)

numberBeamsx=120;% Number of CNTs along perimeter
numberBeamsy=50;% Number of CNTs to be modeled breadthwise

plotint=10;

ang_stdev=5;% Choose 0,5,10, percent for parametric study
rate_stdev=10;% Choose 0, 10, 20 percent for parametric study
avgRate=60e-9;% Average growth per time step (meters)
numberBeams=numberBeamsx*numberBeamsy;% Total number of Beams on the Substrate
ContinuousPlot=1;% 0=plotting off.  1=plotting on
PeriodicBoundary=0;% 0=off, 1=on
