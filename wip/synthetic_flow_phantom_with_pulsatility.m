
tic;

%% studyPath

sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Synthetic_Flow_Phantom';
cd(sp)

% dataDir = pwd;
saveDataDir = 'bSSFP_6pipes_w_pulsatility/data/';
mkdir(saveDataDir);


%% PHANTOM VOLUME PROPERTIES

% volume in px
PIX.x = 128; PIX.y = 128; PIX.z = 80;

% size of voxel (mm)
VOX.x = 1.5; VOX.y = 1.5; VOX.z = 3.0;

% resolution (mm/px)
RES.x = VOX.x/PIX.x; RES.y = VOX.y/PIX.y; RES.z = VOX.z/PIX.z;

% FOV (mm)
FOV.x = VOX.x*PIX.x; FOV.y = VOX.y*PIX.y; FOV.z = VOX.z*PIX.z;

% define flow velocities
PIPE1.vel = 100e-2; %metres/sec
PIPE2.vel = 80e-2;
PIPE3.vel = 25e-2;

PIPE4.vel = -100e-2;
PIPE5.vel = -80e-2;
PIPE6.vel = -25e-2;


%% CREATE PULSATILITY FUNCTIONS
nPhases = 30;  % Basically TR?
x = linspace(1,20,nPhases);
y = normpdf(x,5,1);
yfac = 1/max(y);
pulseProf1 = yfac .* y; %sharp, arterial-ish profile
pulseProf1(pulseProf1 < 0.3) = 0.3; %slowest flow value
% figure; plot(x,pulseProf1,'o-'); axis([1 20, 0 1]);

y = normpdf(x,10,9);
yfac = 1/max(y);
pulseProf2 = circshift(yfac .* y,0); % more constant, venous profile
% figure; plot(x,pulseProf2,'o-'); axis([1 20, 0 1]);


%% CREATE PIPES
%%% PIPE 1 - z-direction
PIPE1.xsize = round(0.5*PIX.x/10);
PIPE1.ysize = round(0.5*PIX.y/10);
PIPE1.xpos = round(PIX.x * 0);
PIPE1.ypos = round(PIX.y * 0);

gauss = customgauss([PIX.x,PIX.y], PIPE1.xsize, PIPE1.ysize, 0, 0, PIPE1.vel, [PIPE1.xpos,PIPE1.ypos]);
outsidePipeIdx = find(gauss < 0.2*PIPE1.vel); % clip minimum flow at 20% of peak
gauss(outsidePipeIdx) = 0;
% imtar(gauss);

Vz = repmat(gauss,1,1,PIX.z); %extend gauss to make a pipe!

% Make pipe pulsatile
for nn = 1:nPhases; Vz1(:,:,:,nn) = Vz .* pulseProf1(nn); end

%%% PIPE 2 - x-direction
PIPE2.xsize = round(0.5*PIX.x/12);
PIPE2.ysize = round(0.5*PIX.y/12);
PIPE2.xpos = round(PIX.x * -0.3);
PIPE2.ypos = round(PIX.y * -0.35);

gauss = customgauss([PIX.x,PIX.y], PIPE2.xsize, PIPE2.ysize, 0, 0, PIPE2.vel, [PIPE2.xpos,PIPE2.ypos]);
outsidePipeIdx = find(gauss < 0.2*PIPE2.vel); % clip minimum flow at 20% of peak
gauss(outsidePipeIdx) = 0;
gauss(:,81:end)=[];
% imtar(gauss);

Vx = repmat(gauss,1,1,PIX.x); %extend gauss to make a pipe!

% Make pipe pulsatile
for nn = 1:nPhases; Vx1(:,:,:,nn) = Vx .* pulseProf1(nn); end

%%% PIPE 3 - y-direction
PIPE3.xsize = round(0.5*PIX.x/10);
PIPE3.ysize = round(0.5*PIX.y/10);
PIPE3.xpos = round(PIX.x * 0);
PIPE3.ypos = round(PIX.y * 0.4);

gauss = customgauss([PIX.x,PIX.y], PIPE3.xsize, PIPE3.ysize, 0, 0, PIPE3.vel, [PIPE3.xpos,PIPE3.ypos]);
outsidePipeIdx = find(gauss < 0.2*PIPE3.vel); % clip minimum flow at 20% of peak
gauss(outsidePipeIdx) = 0;
gauss(81:end,:)=[];
% imtar(gauss);

Vy = repmat(gauss,1,1,PIX.y); %extend gauss to make a pipe!

% Make pipe pulsatile
for nn = 1:nPhases; Vy1(:,:,:,nn) = Vy .* pulseProf1(nn); end

%%% PIPE 4 - negative z-direction
PIPE4.xsize = round(0.5*PIX.x/10);
PIPE4.ysize = round(0.5*PIX.y/10);
PIPE4.xpos = round(PIX.x * 0);
PIPE4.ypos = round(PIX.y * -0.2);

gauss = customgauss([PIX.x,PIX.y], PIPE4.xsize, PIPE4.ysize, 0, 0, PIPE4.vel, [PIPE4.xpos,PIPE4.ypos]);
outsidePipeIdx = find(abs(gauss) < 0.2*abs(PIPE4.vel)); % clip minimum flow at 20% of peak
gauss(outsidePipeIdx) = 0;
% imtar(gauss);

Vz = repmat(gauss,1,1,PIX.z); %extend gauss to make a pipe!

% Make pipe pulsatile
for nn = 1:nPhases; Vz2(:,:,:,nn) = Vz .* pulseProf2(nn); end

%%% PIPE 5 - negative x-direction
PIPE5.xsize = round(0.5*PIX.x/12);
PIPE5.ysize = round(0.5*PIX.y/12);
PIPE5.xpos = round(PIX.x * -0.3);
PIPE5.ypos = round(PIX.y * -0.18);

gauss = customgauss([PIX.x,PIX.y], PIPE5.xsize, PIPE5.ysize, 0, 0, PIPE5.vel, [PIPE5.xpos,PIPE5.ypos]);
outsidePipeIdx = find(abs(gauss) < 0.2*abs(PIPE5.vel)); % clip minimum flow at 20% of peak
gauss(outsidePipeIdx) = 0;
gauss(:,81:end)=[];
% imtar(gauss);

Vx = repmat(gauss,1,1,PIX.x); %extend gauss to make a pipe!

% Make pipe pulsatile
for nn = 1:nPhases; Vx2(:,:,:,nn) = Vx .* pulseProf2(nn); end

%%% PIPE 6 - negative y-direction
PIPE6.xsize = round(0.5*PIX.x/10);
PIPE6.ysize = round(0.5*PIX.y/10);
PIPE6.xpos = round(PIX.x * 0);
PIPE6.ypos = round(PIX.y * 0.18);

gauss = customgauss([PIX.x,PIX.y], PIPE6.xsize, PIPE6.ysize, 0, 0, PIPE6.vel, [PIPE6.xpos,PIPE6.ypos]);
outsidePipeIdx = find(abs(gauss) < 0.2*abs(PIPE6.vel)); % clip minimum flow at 20% of peak
gauss(outsidePipeIdx) = 0;
gauss(81:end,:)=[];
% imtar(gauss);

Vy = repmat(gauss,1,1,PIX.y); %extend gauss to make a pipe!

% Make pipe pulsatile
for nn = 1:nPhases; Vy2(:,:,:,nn) = Vy .* pulseProf2(nn); end


%% Permute all pipes into same volume
Vz1 = Vz1;
Vz2 = Vz2;
Vx1 = permute(Vx1,[1,3,2,4]);
Vx2 = permute(Vx2,[1,3,2,4]);
Vy1 = permute(Vy1,[3,2,1,4]);
Vy2 = permute(Vy2,[3,2,1,4]);

Vx = Vx1 + Vx2;
Vy = Vy1 + Vy2;
Vz = Vz1 + Vz2;

clear Vx1 Vx2 Vy1 Vy2 Vz1 Vz2


%% CREATE FLOW PHANTOM VELOCITY VOLUME
V = Vz + Vx + Vy;
% implay_RR(V(:,:,60,:));


%% VIEW PIPES (V)
% h = vol3d('cdata',abs(V(:,:,:,1)),'texture','3D'); %abs purely for visualisation because of alphamap
% view(3); 
% % Update view since 'texture' = '2D'
% vol3d(h);
% set(gca,'CameraViewAngleMode','manual');
% axis square;
% xlabel('x-axis (AP)');
% ylabel('y-axis (RL)');
% zlabel('z-axis (FH)');


%% SCANNER VOLUME

%+x = AP
SCN.a = 0;
SCN.p = 400;

%+y = RL
SCN.r = 0;
SCN.l = 400;

%+z = FH
SCN.f = 0;
SCN.h = 600;

% centre points of scanner axes
SCN.xCentre = SCN.p/2;
SCN.yCentre = SCN.l/2;
SCN.zCentre = SCN.h/2;


%% STACK
% Must create slice properties before working out the phase values of the
% flow phantom because PHI depends on the bSSFP sequence. i.e.: different
% gradient orientations will result in different PHI.

% slice dimensions
xl = PIX.x; % i
yl = PIX.y; % j
zl = 100;    % k (no. of slices)
dyn = 5;   % no. of dynamics 

STACK = zeros(xl,yl,zl,dyn);

% rotation
% NB: non-rotated slice is defined as transverse within scanner, i.e:
% vector normal points in z-direction (FH)
xTheta = 45;
yTheta = 0;
zTheta = 0;

% centre point of slice
xcp = round(xl/2);
ycp = round(yl/2);
zcp = round(zl/2);

% translation (from origin)
tx_offset = 10;
ty_offset = 5;
tz_offset = 0;

tx = tx_offset + SCN.xCentre - xcp + 1; % +1 seems correct to line-up with PHIpad
ty = ty_offset + SCN.yCentre - ycp + 1;
tz = tz_offset + SCN.zCentre - zcp + 1;

% scaling (voxel size/slice thickness)
sx = 1;
sy = 1;
sz = 1;


%% DETERMINE AFFINE MATRIX

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


%% GRADIENT MOMENTS

% bSSFP slice gradient moment definitions
Gro = 9.95;
Gpe = 0;
Gsl = -10.98;
% Gscn = [Gro, Gpe, Gsl]; % no First Moment scaling (values as in scanner GVE)
Gscn = [Gro, Gpe, Gsl] .* (1e-3)^2; %msec^2.mT/m --- First Moment scaling into correct units

% % QFlow slice gradient moment definitions
% warning('QFlow velocity encoding ...');
% Gro = 0;
% Gpe = 0;
% Gsl = 5; % guessed at this - can be more precise using GVE, if needs be
% % Gscn = [Gro, Gpe, Gsl]; % no First Moment scaling (values as in scanner GVE)
% Gscn = [Gro, Gpe, Gsl] .* (1e-3)^2; %msec^2.mT/m --- First Moment scaling into correct units


% transform from xyz to scanner coordinates
Cprime = [0 -1 0; 1 0 0; 0 0 1];
G = Gscn * Cprime;

% rotate gradient moments with slice rotation
G = R(1:3,1:3) * G';

GAMMA = 42577; %Hz/mT
% GAMMA = 2 * pi * 42577; %Hz/mT


%% PULSATILITY LOOP - acquire one slice per loop

% Pulsatility counter to keep track of where we are in RR-interval
phaseCtr = repmat(1:nPhases,1,10000); phaseCtr = phaseCtr(1:zl*dyn);
Ctr = 1;

disp('Applying slices to velocity volume ...');

for ZZ = 1:zl
    
    for DD = 1:dyn
    
    currPhase = phaseCtr(Ctr);
    
    %% FLOW PHANTOM PHASE VOLUME
    PHI = GAMMA .* ( G(1).*Vx(:,:,:,currPhase) + G(2).*Vy(:,:,:,currPhase) + G(3).*Vz(:,:,:,currPhase) );

    % pad around PHI with zeros:
    PHIpad = padarray(PHI,[(SCN.p-size(PHI,1))/2 , (SCN.l-size(PHI,2))/2 , (SCN.h-size(PHI,3))/2],0,'both');

    % check to see if PHIpad ranges from -ve to 0
    % if so, invert, because vol3d likes value of empty pixels < render region
    % flagInvertPHIpad = 'no';
    % if max(PHIpad(:)) == 0
    %     PHIpad = -1 * PHIpad;
    %     flagInvertPHIpad = 'yes';
    % end


    %% View VOLUME
    % NB: this must be rendered for surf to work when extracting slices

    figure;
    ax1 = axes;
    h = vol3d('cdata',PHIpad,'texture','2D');
    view(3); 
    % Update view since 'texture' = '2D'
    vol3d(h);
    am = abs(linspace(min(PHIpad(:)),max(PHIpad(:)),64)); % zero = transparent
    alphamap(am);
    axis image;
    ax1.CameraViewAngleMode = 'manual';
    grid on;
    hold on;

%     xlabel('x-axis (AP)');
%     ylabel('y-axis (RL)');
%     zlabel('z-axis (FH)');


    %% MAKE SLICES + APPLY TO VELOCITY VOLUME

    % make coordinate grid to transform with AFF
    placeholderSlice = ones([xl,yl,1]); % initalise slices for surf
    % imtar(placeholderSlice(:,:,1));

    % make i,j,k slice coordinate grid
    [gx,gy,gz]=meshgrid(1:xl,1:yl,ZZ);
    g1 = ones(size(gx));


    % apply affine transformation to convert i,j,k into world coordinates
    for xx = 1:yl
        for yy = 1:xl
            for zz = 1

                % get coordinate
                currCoord = [gx(xx,yy,zz);
                             gy(xx,yy,zz);
                             gz(xx,yy,zz);
                             1];

                % apply affine to coordinate
                currAFFcoord = AFF*currCoord;

                % save transformed coordinate
                cx(xx,yy,zz) = currAFFcoord(1);
                cy(xx,yy,zz) = currAFFcoord(2);
                cz(xx,yy,zz) = currAFFcoord(3);

            end
        end
    end

    % permute: because of MATLAB x/y convention
    placeholderSlicePermute = permute(placeholderSlice,[1,2,3]);

    % get transformed slice coordinates using surf
    for ii = ZZ

        hslice = surf(cx,cy,cz,placeholderSlicePermute);

        %get the coordinates
        xslice(:,:,ii) = get(hslice,'XData');
        yslice(:,:,ii) = get(hslice,'YData');
        zslice(:,:,ii) = get(hslice,'ZData');

        delete(hslice); % visualisation comes later

    end
%     view(3);

%     title(['x = ' num2str(xTheta) ...
%            ' y = ' num2str(yTheta) ...
%            ' z = ' num2str(zTheta) ' (degrees)']);

    % make a grid the size of the Scanner VOLUME
    [xvol,yvol,zvol] = meshgrid( 1:size(PHIpad,1) , 1:size(PHIpad,2), 1:size(PHIpad,3) );

    % put the transformed slice into Scanner VOLUME grid and extract VOLUME values
    for ii = ZZ

        hslice = slice(xvol,yvol,zvol,PHIpad,xslice(:,:,ii),yslice(:,:,ii),zslice(:,:,ii));
        hslice.FaceColor = 'texturemap';
        hslice.EdgeColor = 'none';
    %     drawnow;
        STACK(:,:,ii,DD) = hslice.CData;

        disp(['Slice ' num2str(ii) ' , Dynamic ' num2str(DD) ' ...']);

        delete(hslice);
    end

    close;
    
    % update pulsatility counter (ie: next phase in cardiac cycle)
    Ctr = Ctr + 1;

    end % dynamics loop
    
end % end slices loop

disp('DONE ...');


%% BIT OF A BODGE --- VISUALISE STACKS THROUGH VOLUME
% NB: this is purely for visualisation.
% The code above extracts the slice data.
% Separate code is needed for visualisation because plotting surf on and
% existing vol3d, when surf contains negative AND positive numbers causes
% havoc with the alphamaps.
% TO DO: fix this so that I can view the stacks through the velocity
% volume, without having to render abs(PHIpad).

% disp('Rendering visualisation of stack through velocity volume ...');
% 
% figure;
% ax1 = axes;
% 
% %%%%% crucial here: abs() --- needed to for visualisation
% % without abs, alphamap doesn't render zeros as transparent when using surf
% h = vol3d('cdata',abs(PHIpad),'texture','2D'); 
% view(3); 
% % Update view since 'texture' = '2D'
% vol3d(h);
% axis image;
% ax1.CameraViewAngleMode = 'manual';
% grid on;
% hold on;
% 
% xlabel('x-axis (AP)');
% ylabel('y-axis (RL)');
% zlabel('z-axis (FH)');
% 
% title(['x = ' num2str(xTheta) ...
%        ' y = ' num2str(yTheta) ...
%        ' z = ' num2str(zTheta) ' (degrees)']);
% 
% for ii = 1:zl
%     hslice = slice(xvol,yvol,zvol,abs(PHIpad),xslice(:,:,ii),yslice(:,:,ii),zslice(:,:,ii));
%     hslice.FaceColor = 'k';
%     hslice.EdgeColor = 'none';
%     drawnow;
%     
%     disp(['Slice ' num2str(ii) ' ...']);
% end
% 
% % CMAP = colormap;
% % CBAR = colorbar;
% % cSiz = size(CMAP,1);
% % cMin = CBAR.Limits(1);
% % cMax = CBAR.Limits(2);
% % am2 = linspace(cMin,cMax,cSiz);
% % alphamap(abs(am2));
% 
% disp('DONE ...');


%% ADD NOISE TO STACK

mu = mean(STACK(:));
% sigma = std(V(:));
% measured bSSFP flow phantom:
% - ROI in water = 137
% - ROI in noise around phantom = 2
% - 137/2 = 68.5
sigma = 0.05 ./ 13; % this gives SNR = 68 in synthetic phantom image

N = normrnd(0,sigma,size(STACK));
STACKnoisy = STACK + N;
STACK = STACKnoisy;

%% Save stack as nii

nii = make_nii(permute(STACK,[2,1,3,4])); %NB: permute as MATLAB flips x/y

% nb: 4 below is super important!! Otherwise .nii doesn't render in-plane in MITK
nii.hdr.dime.dim = [4 size(nii.img,1) size(nii.img,2) size(nii.img,3) dyn 1 1 1];
nii.hdr.dime.pixdim(4) = 1; %give some thickness...

nii.hdr.hist.srow_x = AFF(2,:); %NB: MATLAB flips x/y - need AFF to match phi permutation
nii.hdr.hist.srow_y = AFF(1,:);
nii.hdr.hist.srow_z = AFF(3,:);
nii.hdr.hist.sform_code = 1;

%%% Offset corrections
% TO DO: fix this!
% -2 offset required for some rotations... not sure why?!
% -1 offset required for 45 degree rotation...

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


% TO DO: automate pixdim once I implement voxel scaling

fileName = ['stack_rx' num2str(xTheta) ...
                 '_ry' num2str(yTheta) ...
                 '_rz' num2str(zTheta) ...
                 '_tx' num2str(tx_offset) ...
                 '_ty' num2str(ty_offset) ...
                 '_tz' num2str(tz_offset) ...
                 '_slices' num2str(zl) ...
                 '_dyn' num2str(dyn) ...
                 '_pulsatile.nii.gz'];
             
save_nii(nii,[saveDataDir fileName]);

disp(['Saving stack as ... ' fileName]);
disp('DONE ...');


%% Save GRADIENT MOMENTS

disp(['Saving gradient moments as ... ' fileName(1:end-7) '_grad_moments.txt']);

% revert G back to scale on scanner
G = G ./ (1e-3)^2;

Gunit = G./norm(G);
Gnorm = norm(G);

fileID = fopen([saveDataDir fileName(1:end-7) '_scanner_grad_moments.txt'],'w');
fprintf(fileID,'%.4f\t',G(:));
fclose(fileID);

fileID = fopen([saveDataDir fileName(1:end-7) '_grad_moment_dirs.txt'],'w');
fprintf(fileID,'%.4f\t',Gunit(:));
fclose(fileID);

fileID = fopen([saveDataDir fileName(1:end-7) '_grad_moment_vals.txt'],'w');
fprintf(fileID,'%.4f\t',Gnorm(:));
fclose(fileID);

disp('DONE ...');


%% Save FIGURE

% disp(['Saving stacks/velocity visualisation as ... ' fileName(1:end-7) '.fig']);
% saveas(gcf,[saveDataDir fileName(1:end-7) '.fig']);
% 
% disp(['Saving stacks/velocity visualisation render as ... ' fileName(1:end-7) '.png']);
% saveas(gcf,[saveDataDir fileName(1:end-7) '.png']);
% 
% disp('DONE ...');
% 
% 
% toc;


%% SAVE VELOCITY VOLUME

% clear niiV
% 
% Vpad = padarray(V,[(SCN.p-size(V,1))/2 , (SCN.l-size(V,2))/2 , (SCN.h-size(V,3))/2],0,'both');
% 
% % Vpad = padarray(Vx,[(SCN.p-size(V,1))/2 , (SCN.l-size(V,2))/2 , (SCN.h-size(V,3))/2],0,'both');
% % Vpad = padarray(Vy,[(SCN.p-size(V,1))/2 , (SCN.l-size(V,2))/2 , (SCN.h-size(V,3))/2],0,'both');
% % Vpad = padarray(Vz,[(SCN.p-size(V,1))/2 , (SCN.l-size(V,2))/2 , (SCN.h-size(V,3))/2],0,'both');
% 
% niiV = make_nii(Vpad);
% niiV.hdr.dime.dim    = [3 size(niiV.img,1) size(niiV.img,2) size(niiV.img,3) 1 1 1 1];
% niiV.hdr.dime.pixdim = [0 1 1 1 1 1 1 1]; %make voxels 1mm isotropic
% 
% % srow(4) = 0
% niiV.hdr.hist.srow_x = [1 0 0 0];
% niiV.hdr.hist.srow_y = [0 1 0 0];
% niiV.hdr.hist.srow_z = [0 0 1 0];
% niiV.hdr.hist.sform_code = 1;
% 
% save_nii(niiV,'VEL_VOLUME_srow4_equal0.nii.gz');
% 
% disp('Made Velocity Volume .nii ...');


%% Save permutations

%%%%%%%%%% NB: need to correct below. Moved from STACK.phi to STACK


%- testing permutations of STACK

% STACK.phi123 = permute(STACK,[1,2,3]);
% STACK.phi132 = permute(STACK,[1,3,2]);
% 
% STACK.phi213 = permute(STACK,[2,1,3]);
% STACK.phi231 = permute(STACK,[2,3,1]);
% 
% STACK.phi312 = permute(STACK,[3,1,2]);
% STACK.phi321 = permute(STACK,[3,2,1]);
% 
% %%% 1) perm123
% 
% nii = make_nii(STACK.phi123);
% 
% % nb: 3 below is super important!! Otherwise .nii doesn't render in-plane in MITK
% nii.hdr.dime.dim    = [3 size(nii.img,1) size(nii.img,2) size(nii.img,3) 1 1 1 1];
% % nii.hdr.dime.dim    = [3 size(STACK.phi,1) size(STACK.phi,2) size(STACK.phi,3) 1 1 1 1];
% nii.hdr.dime.pixdim(4) = 1; %give some thickness... 
% 
% nii.hdr.hist.srow_x = AFF(2,:);
% nii.hdr.hist.srow_y = AFF(1,:);
% nii.hdr.hist.srow_z = AFF(3,:);
% nii.hdr.hist.sform_code = 1;
% 
% %TO DO: fill in pixdim and glmax/glmin etc.
% 
% fileName = ['stack_x' num2str(xTheta) ...
%                  '_y' num2str(yTheta) ...
%                  '_z' num2str(zTheta) ...
%                  '_tx' num2str(manualtx) ...
%                  '_ty' num2str(manualty) ...
%                  '_tz' num2str(manualtz) ...
%                  '_test_perm123.nii.gz'];
%              
% save_nii(nii,fileName);
% 
% disp(['Saved stack as ... ' fileName]);
% 
% 
% %%% 2) perm132
% 
% nii = make_nii(STACK.phi132);
% 
% % nb: 3 below is super important!! Otherwise .nii doesn't render in-plane in MITK
% nii.hdr.dime.dim    = [3 size(nii.img,1) size(nii.img,2) size(nii.img,3) 1 1 1 1];
% % nii.hdr.dime.dim    = [3 size(STACK.phi,1) size(STACK.phi,2) size(STACK.phi,3) 1 1 1 1];
% nii.hdr.dime.pixdim(4) = 1; %give some thickness... 
% 
% nii.hdr.hist.srow_x = AFF(2,:);
% nii.hdr.hist.srow_y = AFF(1,:);
% nii.hdr.hist.srow_z = AFF(3,:);
% nii.hdr.hist.sform_code = 1;
% 
% %TO DO: fill in pixdim and glmax/glmin etc.
% 
% fileName = ['stack_x' num2str(xTheta) ...
%                  '_y' num2str(yTheta) ...
%                  '_z' num2str(zTheta) ...
%                  '_tx' num2str(manualtx) ...
%                  '_ty' num2str(manualty) ...
%                  '_tz' num2str(manualtz) ...
%                  '_test_perm132.nii.gz'];
%              
% save_nii(nii,fileName);
% 
% disp(['Saved stack as ... ' fileName]);
% 
% 
% %%% 3) perm213
% 
% nii = make_nii(STACK.phi213);
% 
% % nb: 3 below is super important!! Otherwise .nii doesn't render in-plane in MITK
% nii.hdr.dime.dim    = [3 size(nii.img,1) size(nii.img,2) size(nii.img,3) 1 1 1 1];
% % nii.hdr.dime.dim    = [3 size(STACK.phi,1) size(STACK.phi,2) size(STACK.phi,3) 1 1 1 1];
% nii.hdr.dime.pixdim(4) = 1; %give some thickness... 
% 
% nii.hdr.hist.srow_x = AFF(2,:);
% nii.hdr.hist.srow_y = AFF(1,:);
% nii.hdr.hist.srow_z = AFF(3,:);
% nii.hdr.hist.sform_code = 1;
% 
% %TO DO: fill in pixdim and glmax/glmin etc.
% 
% fileName = ['stack_x' num2str(xTheta) ...
%                  '_y' num2str(yTheta) ...
%                  '_z' num2str(zTheta) ...
%                  '_tx' num2str(manualtx) ...
%                  '_ty' num2str(manualty) ...
%                  '_tz' num2str(manualtz) ...
%                  '_test_perm213.nii.gz'];
%              
% save_nii(nii,fileName);
% 
% disp(['Saved stack as ... ' fileName]);
% 
% 
% %%% 4) perm231
% 
% nii = make_nii(STACK.phi231);
% 
% % nb: 3 below is super important!! Otherwise .nii doesn't render in-plane in MITK
% nii.hdr.dime.dim    = [3 size(nii.img,1) size(nii.img,2) size(nii.img,3) 1 1 1 1];
% % nii.hdr.dime.dim    = [3 size(STACK.phi,1) size(STACK.phi,2) size(STACK.phi,3) 1 1 1 1];
% nii.hdr.dime.pixdim(4) = 1; %give some thickness... 
% 
% nii.hdr.hist.srow_x = AFF(2,:);
% nii.hdr.hist.srow_y = AFF(1,:);
% nii.hdr.hist.srow_z = AFF(3,:);
% nii.hdr.hist.sform_code = 1;
% 
% %TO DO: fill in pixdim and glmax/glmin etc.
% 
% fileName = ['stack_x' num2str(xTheta) ...
%                  '_y' num2str(yTheta) ...
%                  '_z' num2str(zTheta) ...
%                  '_tx' num2str(manualtx) ...
%                  '_ty' num2str(manualty) ...
%                  '_tz' num2str(manualtz) ...
%                  '_test_perm231.nii.gz'];
%              
% save_nii(nii,fileName);
% 
% disp(['Saved stack as ... ' fileName]);
% 
% 
% %%% 5) perm312
% 
% nii = make_nii(STACK.phi312);
% 
% % nb: 3 below is super important!! Otherwise .nii doesn't render in-plane in MITK
% nii.hdr.dime.dim    = [3 size(nii.img,1) size(nii.img,2) size(nii.img,3) 1 1 1 1];
% % nii.hdr.dime.dim    = [3 size(STACK.phi,1) size(STACK.phi,2) size(STACK.phi,3) 1 1 1 1];
% nii.hdr.dime.pixdim(4) = 1; %give some thickness... 
% 
% nii.hdr.hist.srow_x = AFF(2,:);
% nii.hdr.hist.srow_y = AFF(1,:);
% nii.hdr.hist.srow_z = AFF(3,:);
% nii.hdr.hist.sform_code = 1;
% 
% %TO DO: fill in pixdim and glmax/glmin etc.
% 
% fileName = ['stack_x' num2str(xTheta) ...
%                  '_y' num2str(yTheta) ...
%                  '_z' num2str(zTheta) ...
%                  '_tx' num2str(manualtx) ...
%                  '_ty' num2str(manualty) ...
%                  '_tz' num2str(manualtz) ...
%                  '_test_perm312.nii.gz'];
%              
% save_nii(nii,fileName);
% 
% disp(['Saved stack as ... ' fileName]);
% 
% 
% %%% 6) perm321
% 
% nii = make_nii(STACK.phi321);
% 
% % nb: 3 below is super important!! Otherwise .nii doesn't render in-plane in MITK
% nii.hdr.dime.dim    = [3 size(nii.img,1) size(nii.img,2) size(nii.img,3) 1 1 1 1];
% % nii.hdr.dime.dim    = [3 size(STACK.phi,1) size(STACK.phi,2) size(STACK.phi,3) 1 1 1 1];
% nii.hdr.dime.pixdim(4) = 1; %give some thickness... 
% 
% nii.hdr.hist.srow_x = AFF(2,:);
% nii.hdr.hist.srow_y = AFF(1,:);
% nii.hdr.hist.srow_z = AFF(3,:);
% nii.hdr.hist.sform_code = 1;
% 
% %TO DO: fill in pixdim and glmax/glmin etc.
% 
% fileName = ['stack_x' num2str(xTheta) ...
%                  '_y' num2str(yTheta) ...
%                  '_z' num2str(zTheta) ...
%                  '_tx' num2str(manualtx) ...
%                  '_ty' num2str(manualty) ...
%                  '_tz' num2str(manualtz) ...
%                  '_test_perm321.nii.gz'];
%              
% save_nii(nii,fileName);
% 
% disp(['Saved stack as ... ' fileName]);


