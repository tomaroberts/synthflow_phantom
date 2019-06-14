%SYNTHFLOW_ACQUIRE   acquires a stack of slices through a water flow phantom
%
%   SYNTHFLOW_ACQUIRE()
%
%   requires: parameter definitions in wrapper script: synthflow_phantom
%
%   output:
%       simulated flow phantom .nii
%       stack of phase slices through phantom
%       gradient first moments .txt files associated with stack
%       .fig/.png files
%
%   see also: synthflow_phantom, synthflow_waterpipe

% Tom Roberts (t.roberts@kcl.ac.uk)

%% NOTES

%   MATLAB visualisation of slices in flow phantom slightly buggy. Alphamap
%   means regions of phase = 0 are rendered transparent. This does not
%   impact the final .nii files - purely a visualisation issue.

% TODO:
%   - currently 6 orthogoal pipes as default - make more general
%   - implement voxel scaling - currently only 1 mm isotropic
%   - slight misalignment in final .nii due to rotation (see below)


%% Create Pipes
PIPE1.V = synthflow_waterpipe( PIPE1.vel, [PIX.x, PIX.y, PIX.z] );

PIPE2.V = synthflow_waterpipe( PIPE2.vel, [PIX.x, PIX.y, PIX.z] , ...
    'pipeSize', [round(0.5*PIX.x/12), round(0.5*PIX.y/12)], ...
    'pipeDir', 'x', ...
    'pipePos',  [round(PIX.x * -0.3), round(PIX.y * -0.35)] );

PIPE3.V = synthflow_waterpipe( PIPE3.vel, [PIX.x, PIX.y, PIX.z] , ...
    'pipeSize', [round(0.5*PIX.x/10), round(0.5*PIX.y/10)], ...
    'pipeDir', 'y', ...
    'pipePos',  [round(PIX.x * 0), round(PIX.y * 0.4)] );

PIPE4.V = synthflow_waterpipe( PIPE4.vel, [PIX.x, PIX.y, PIX.z] , ...
    'pipeSize', [round(0.5*PIX.x/10), round(0.5*PIX.y/10)], ...
    'pipeDir', 'z', ...
    'pipePos',  [round(PIX.x * 0), round(PIX.y * -0.2)] );

PIPE5.V = synthflow_waterpipe( PIPE5.vel, [PIX.x, PIX.y, PIX.z] , ...
    'pipeSize', [round(0.5*PIX.x/12), round(0.5*PIX.y/12)], ...
    'pipeDir', 'x', ...
    'pipePos',  [round(PIX.x * -0.3), round(PIX.y * -0.18)] );

PIPE6.V = synthflow_waterpipe( PIPE6.vel, [PIX.x, PIX.y, PIX.z] , ...
    'pipeSize', [round(0.5*PIX.x/10), round(0.5*PIX.y/10)], ...
    'pipeDir', 'y', ...
    'pipePos',  [round(PIX.x * 0), round(PIX.y * 0.18)] );

% permute x/y pipes
PIPE2.V = permute(PIPE2.V,[1,3,2]);
PIPE3.V = permute(PIPE3.V,[3,2,1]);
PIPE5.V = permute(PIPE5.V,[1,3,2]);
PIPE6.V = permute(PIPE6.V,[3,2,1]);


%% Create Flow Phantom Velocity Volume and Component Volumes
Vz = PIPE1.V + PIPE4.V;
Vx = PIPE2.V + PIPE5.V;
Vy = PIPE3.V + PIPE6.V;
V = Vx + Vy + Vz;


%% View Flow Phantom Velocity Volume
h1 = vol3d('cdata',V,'texture','3D');
view(45+180,30);
vol3d(h1);
am = abs(linspace(min(V(:)),max(V(:)),64)); am(31:33) = 0; % zero = transparent
alphamap(am);
set(gca,'CameraViewAngleMode','manual');
axis image;
colormap(jet(512)); c = colorbar; c.Label.String = 'Velocity [cm/s]';
title('Flow Phantom Velocity Volume');
xlabel('x-axis (AP)');
ylabel('y-axis (RL)');
zlabel('z-axis (FH)');


%% Define Centre Points of Scanner Axes
SCN.xCentre = SCN.p/2;
SCN.yCentre = SCN.l/2;
SCN.zCentre = SCN.h/2;


%% Slice Affine Matrix Processing
% NB: must create slice in 3D space first because phase of slice is dependent
% on its orientation relative to the scanner gradients

% centre point of slice
xcp = round(xl/2);
ycp = round(yl/2);
zcp = round(zl/2);

% slice translation correction
tx = tx_offset + SCN.xCentre - xcp + 1; % +1 seems correct to line-up with PHIpad
ty = ty_offset + SCN.yCentre - ycp + 1;
tz = tz_offset + SCN.zCentre - zcp + 1;

% scaling (convert from px to mm)
sx = 1;
sy = 1;
sz = 1;

% rotation
rx = rotxtar(deg2rad(xTheta));
ry = rotytar(deg2rad(yTheta));
rz = rotztar(deg2rad(zTheta));

% translation matrix
T = [1 0 0 tx; ...
     0 1 0 ty; ...
     0 0 1 tz; ...
     0 0 0 1];

% arbitrary point matrix (ie: centre of slice)
P = [1 0 0 xcp; ...
     0 1 0 ycp; ...
     0 0 1 zcp; ...
     0 0 0 1];

% minus arbitrary point (necessary for rotating around a point)
% --- see: http://www.euclideanspace.com/maths/geometry/affine/aroundPoint/
Pminus = [1 0 0 -xcp; ...
          0 1 0 -ycp; ...
          0 0 1 -zcp; ...
          0 0 0 1]; 

% scaling matrix
S = [sx 0 0 0;
     0 sy 0 0;
     0 0 sz 0;
     0 0 0 1];

% rotation matrices
rx(:,4) = 0; ry(:,4) = 0; rz(:,4) = 0;
rx(4,:) = [0 0 0 1]; ry(4,:) = [0 0 0 1]; rz(4,:) = [0 0 0 1];
R = rx*ry*rz;

% matlab coordinates adjustment matrix
% --- see: http://gru.stanford.edu/doku.php/mrTools/coordinateTransforms
M = [1 0 0 -1; ...
     0 1 0 -1; ...
     0 0 1 -1; ...
     0 0 0 1];
 
% final affine transformation matrix
% --- followed here: http://gru.stanford.edu/doku.php/mrTools/coordinateTransforms
% --- and here: http://www.euclideanspace.com/maths/geometry/affine/aroundPoint/

% AFF = M*R*S*T;            % rotation only around (0,0,0)
AFF = M*T*P*R*S*Pminus;     % rotation around arbitrary point P (ie: centre of slice)


%% Gradient First Moment Processing

% scale gradient first moments into correct units
% currently based on Philips Ingenia 1.5T scanner
Gscn = [Gro, Gpe, Gsl] .* (1e-3)^2; %msec^2.mT/m 
% Gscn = [Gro, Gpe, Gsl]; % no First Moment scaling (values as in scanner GVE)

% transform from xyz to scanner coordinates
Cprime = [0 -1 0; 1 0 0; 0 0 1]; % Philips
G = Gscn * Cprime;

% rotate gradient moments with slice rotation
G = R(1:3,1:3) * G';

GAMMA = 42577; %Hz/mT
% GAMMA = 2 * pi * 42577; %Hz/mT


%% Measure Phase of Flow Phantom in Scanner Bore
%- uses gradient first moments defined above
PHI = GAMMA .* ( G(1).*Vx + G(2).*Vy + G(3).*Vz );

% Set Flow Phantom in Scanner Bore
PHIpad = padarray( PHI, [(SCN.p-size(PHI,1))/2 , ...
    (SCN.l-size(PHI,2))/2 , ...
    (SCN.h-size(PHI,3))/2], ...
    0, 'both' );
    

%% Render Phase Representation of Flow Phantom in Bore
% NB: this must be rendered for surf to work when extracting slices

disp('Scanning Flow Phantom ...');

figure;
ax1 = axes;
h2 = vol3d('cdata',PHIpad,'texture','2D');
view(45+180,30);
vol3d(h2);
am = abs(linspace(min(PHIpad(:)),max(PHIpad(:)),64)); % zero = transparent
alphamap(am);
axis image;
ax1.CameraViewAngleMode = 'manual';
grid on;
hold on;
colormap(gray(512)); c = colorbar; c.Label.String = 'Phase [rads]';

xlabel('x-axis (AP)');
ylabel('y-axis (RL)');
zlabel('z-axis (FH)');

title(['x = ' num2str(xTheta) ...
       ' y = ' num2str(yTheta) ...
       ' z = ' num2str(zTheta) ' (degrees)']);


%% Make Stack of Slices and Put in World Coordinates

% make grid for slices
[gx,gy,gz]=meshgrid(1:xl,1:yl,1:zl);

% apply affine transformation to convert i,j,k into world coordinates
for xx = 1:yl
    for yy = 1:xl
        for zz = 1:zl
            
            % get coordinate
            currCoord = [gx(xx,yy,zz);
                         gy(xx,yy,zz);
                         gz(xx,yy,zz);
                         1];
            
            % apply affine to coordinate
            currAFFcoord = AFF*currCoord;
                     
            % save world coordinates
            cx(xx,yy,zz) = currAFFcoord(1);
            cy(xx,yy,zz) = currAFFcoord(2);
            cz(xx,yy,zz) = currAFFcoord(3);
            
        end
    end
end
  
   
%% Acquire Stack of Slices through Flow Phantom Phase Volume

disp('Acquiring Slices ...');

% make grid for scanner volume
[xvol,yvol,zvol] = meshgrid( 1:size(PHIpad,1) , 1:size(PHIpad,2), 1:size(PHIpad,3) );

% create stack of phase images:
% - transformed slices put into scanner volume and phase extracted for each
STACK = zeros(xl,yl,zl);

for ii = 1:zl
    hslice = slice(xvol,yvol,zvol,PHIpad,cx(:,:,ii),cy(:,:,ii),cz(:,:,ii));
    hslice.FaceColor = 'k';
    hslice.EdgeColor = 'none';
    
    % make finer alphamap and adjust to render zero phase regions as transparent
    if ii == 1
        % colormap
        newCM = linspace(h2.parent.ALim(1),h2.parent.ALim(2),512);
        
        % alphamap
        fract = 10^-6; idx = newCM<=fract; % find zero point in colormap        
        newAM = abs(linspace(h2.parent.ALim(1),-1*h2.parent.ALim(1),512));
        newAM(find(idx,1,'last'):find(idx,1,'last')+1) = 0; alphamap(newAM);
    end
    
    drawnow;
    STACK(:,:,ii) = hslice.CData;
    
    disp(['Slice ' num2str(ii) ' ...']);
end

disp('Stack Acquired ...');


%% Add Random Noise to Images

mu = mean(STACK(:));

% measured SNR in water flow phantom acquired with bSSFP:
% - ROI in water = 137
% - ROI in noise around phantom = 2
% - 137/2 = 68.5

sigma = 0.05 ./ 13; % this gives SNR = 68 in synthetic phantom image

N = normrnd(0,sigma,size(STACK));
STACKnoisy = STACK + N;
STACK = STACKnoisy;


%% Save Stack as .nii

nii = make_nii(permute(STACK,[2,1,3])); %NB: permute as MATLAB flips x/y

% nb: 3 below is super important!! Otherwise .nii doesn't render in-plane in MITK
nii.hdr.dime.dim = [3 size(nii.img,1) size(nii.img,2) size(nii.img,3) 1 1 1 1];
nii.hdr.dime.pixdim(4) = 1; %give some thickness...

nii.hdr.hist.srow_x = AFF(2,:); %NB: MATLAB flips x/y - need AFF to match phi permutation
nii.hdr.hist.srow_y = AFF(1,:);
nii.hdr.hist.srow_z = AFF(3,:);
nii.hdr.hist.sform_code = 1;


%% Rotation Offset Correction
% TODO: improve this
%- Rotations cause small misalignments when viewed in MITK
%- I manually realigned the stacks using small adjustments below
%- I think the origin of this issue is to do with how different software
%  interprets whether grid points correspond to the corners or centre of
%  voxels.

if (xTheta == 45 && yTheta == 0 && zTheta == 0) || (xTheta == 50 && yTheta == 0 && zTheta == 0) 
    warning('Applying -1 offset to srow_x in .nii ... ');
    nii.hdr.hist.srow_x(4) = nii.hdr.hist.srow_x(4) - 1; %offset for 45 degree x-rotated stacks
end

if (xTheta == 90 && yTheta == 0 && zTheta == 0) 
    warning('Applying -2 offset to srow_x in .nii ... ');
    nii.hdr.hist.srow_x(4) = nii.hdr.hist.srow_x(4) - 2; %offset for 90 degree x-rotated stacks
end

if (xTheta == 0 && yTheta == 0 && zTheta == 90)
    warning('Applying -2 offset to srow_y in .nii ... ');
    nii.hdr.hist.srow_y(4) = nii.hdr.hist.srow_y(4) - 2; %offset for 90 degree z-rotated stacks
end

if (xTheta == 0 && yTheta == 90 && zTheta == 0)
    warning('Applying -2 offset to srow_z in .nii ... ');
    nii.hdr.hist.srow_z(4) = nii.hdr.hist.srow_z(4) - 1; %offset for 90 degree y-rotated stacks
end

if (xTheta == 30 && yTheta == 60 && zTheta == 0)
    warning('Applying -1 offset to srow_z in .nii ... ');
    nii.hdr.hist.srow_z(4) = nii.hdr.hist.srow_z(4) - 1;
end

if (xTheta == 0 && yTheta == 0 && zTheta == 150)
    warning('Applying -1 offset to srow_x in .nii ... ');
    nii.hdr.hist.srow_x(4) = nii.hdr.hist.srow_x(4) - 1;
    
    warning('Applying -2 offset to srow_y in .nii ... ');
    nii.hdr.hist.srow_y(4) = nii.hdr.hist.srow_y(4) - 2;
end

if (xTheta == 0 && yTheta == 0 && zTheta == 130)
    warning('Applying -1 offset to srow_x in .nii ... ');
    nii.hdr.hist.srow_x(4) = nii.hdr.hist.srow_x(4) - 1;
    
    warning('Applying -1.5 offset to srow_y in .nii ... ');
    nii.hdr.hist.srow_y(4) = nii.hdr.hist.srow_y(4) - 1.5;
end

if (xTheta == 0 && yTheta == 0 && zTheta == 110)
    warning('Applying -1 offset to srow_x in .nii ... ');
    nii.hdr.hist.srow_x(4) = nii.hdr.hist.srow_x(4) - 0.5;
    
    warning('Applying -1.5 offset to srow_y in .nii ... ');
    nii.hdr.hist.srow_y(4) = nii.hdr.hist.srow_y(4) - 2;
end

fileName = ['stack_rx' num2str(xTheta) ...
                 '_ry' num2str(yTheta) ...
                 '_rz' num2str(zTheta) ...
                 '_tx' num2str(tx_offset) ...
                 '_ty' num2str(ty_offset) ...
                 '_tz' num2str(tz_offset) ...
                 '_slices' num2str(zl) ...
                 '.nii.gz'];
             
save_nii(nii,[saveDataDir '\' fileName]);

disp('Saving stack as ... ');
disp(fileName);
disp('DONE ...');


%% Save Flow Phantom Volume within Scanner

Vpad = padarray(V,[(SCN.p-size(V,1))/2 , (SCN.l-size(V,2))/2 , (SCN.h-size(V,3))/2],0,'both');
% 
% % Vpad = padarray(Vx,[(SCN.p-size(V,1))/2 , (SCN.l-size(V,2))/2 , (SCN.h-size(V,3))/2],0,'both');
% % Vpad = padarray(Vy,[(SCN.p-size(V,1))/2 , (SCN.l-size(V,2))/2 , (SCN.h-size(V,3))/2],0,'both');
% % Vpad = padarray(Vz,[(SCN.p-size(V,1))/2 , (SCN.l-size(V,2))/2 , (SCN.h-size(V,3))/2],0,'both');
% 
niiV = make_nii(Vpad);
niiV.hdr.dime.dim    = [3 size(niiV.img,1) size(niiV.img,2) size(niiV.img,3) 1 1 1 1];
niiV.hdr.dime.pixdim = [0 1 1 1 1 1 1 1]; %make voxels 1mm isotropic

% set srow(4) = 0
niiV.hdr.hist.srow_x = [1 0 0 0];
niiV.hdr.hist.srow_y = [0 1 0 0];
niiV.hdr.hist.srow_z = [0 0 1 0];
niiV.hdr.hist.sform_code = 1;

save_nii(niiV,[saveDataDir '\flow_phantom_vel_volume.nii.gz']);

disp('Saving Flow Phantom Volume as ... ');
disp('flow_phantom_vel_volume.nii.gz');
disp('DONE ...');


%% Save Gradient First Moments to .txt Files

disp('Saving gradient moments as ... ');
disp([fileName(1:end-7) '_grad_moments.txt']);

% revert G back to scale on scanner
G = G ./ (1e-3)^2;

Gunit = G./norm(G);
Gnorm = norm(G);

fileID = fopen([saveDataDir '\' fileName(1:end-7) '_scanner_grad_moments.txt'],'w');
fprintf(fileID,'%.4f\t',G(:));
fclose(fileID);

fileID = fopen([saveDataDir '\' fileName(1:end-7) '_grad_moment_dirs.txt'],'w');
fprintf(fileID,'%.4f\t',Gunit(:));
fclose(fileID);

fileID = fopen([saveDataDir '\' fileName(1:end-7) '_grad_moment_vals.txt'],'w');
fprintf(fileID,'%.4f\t',Gnorm(:));
fclose(fileID);

disp('DONE ...');


%% Save Visualisations

disp('Saving stacks/velocity visualisation as ... ');
disp([fileName(1:end-7) '.fig']);
saveas(gcf,[saveDataDir '\' fileName(1:end-7) '.fig']);

disp('Saving stacks/velocity visualisation render as ... ');
disp([fileName(1:end-7) '.png']);
saveas(gcf,[saveDataDir '\' fileName(1:end-7) '.png']);

disp('DONE ...');


%% End Message
disp('Flow Phantom Acquisition finished ...');


% synthflow_acquire(...)