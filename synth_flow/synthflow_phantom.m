%SYNTHFLOW_PHANTOM   synthetic flow phantom wrapper script
%
%   SYNTHFLOW_PHANTOM( )
%
%   instructions:
%       This is an editable wrapper script for synthflow_acquire.m, which 
%       allows the user to adjust various parameters of the flow phantom 
%       and the simulated scanner, such as the gradient first moments. All 
%       output files are saved in the working directory.
%
%   output:
%       simulated flow phantom .nii
%       stack of phase slices through phantom
%       gradient first moments .txt files associated with stack
%       .fig/.png files
%
%   see also: synthflow_waterpipe, synthflow_acquire

% Tom Roberts (t.roberts@kcl.ac.uk)


%% Save Data Directory
saveDataDir = pwd;


%% Flow Phantom Volume
%- defines a cube containing the flow phantom

% volume in px
PIX.x = 128; PIX.y = 128; PIX.z = 80;


%% Define Peak Pipe Flow Velocities
PIPE1.vel = 100e-2; %metres/sec
PIPE2.vel = 80e-2;
PIPE3.vel = 25e-2;

PIPE4.vel = -100e-2;
PIPE5.vel = -80e-2;
PIPE6.vel = -25e-2;


%% Slice Dimension Properties
%- defines the stack of slices its orientation in the scanner bore

% slice dimensions
xl = PIX.x;
yl = PIX.y;
zl = 3;     % no. of slices

% rotation
% NB: non-rotated slice is defined as transverse within scanner, i.e:
% vector normal points in z-direction (FH)
xTheta = 20;
yTheta = 80;
zTheta = 10;

% translation (from origin)
tx_offset = 0;
ty_offset = 0;
tz_offset = 0;


%% Slice First Moment Properties
%- For conventional SPGR, only Gsl is nonzero
% Gro = 0;
% Gpe = 0;
% Gsl = 10;

%- For bSSFP, Gro and Gsl are nonzero (see: Nielsen 2009).
Gro = 9.95;
Gpe = 0;
Gsl = -10.98;


%% Scanner Volume
%- defines the dimensions of the scanner around the phantom
%- must be bigger than the flow phantom

%+x = AP
SCN.a = 0;
SCN.p = 400; % px

%+y = RL
SCN.r = 0;
SCN.l = 400;

%+z = FH
SCN.f = 0;
SCN.h = 600;


%% Acquire Slices in 'Scanner'
synthflow_acquire;


% synthflow_phantom(...)