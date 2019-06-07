%% synthetic flow phantom - rotate stack, re-analyse velocity
% - 13/03/2019 --- Alena asked me to investigate rotating one of the phase
% stacks along with the gradient moments, then redoing the analytical
% solution for velocity. Think she wants to check it matches MIRTK
%
% 14/03/2019 --- This ended up being quite useful. Rotating the gradient
% moments together with the phase stack worked, as expected. I extended
% this script to examine what happens to the final velocity values if the
% gradient moments are NOT rotated with stack rotation. Used 5/15/25
% degrees rotation of grad_moments. The error increased as expected.
% POWERPOINT file with comparison: 
% C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Synthetic_Flow_Phantom\Gradient_moment_rotation_error_comparison.pptx


%% STEP 1) stacks affine transformed into same grid using:
%- stacks_to_common_grid_script.bash
%- this is just so that the voxels are all aligned
%- NB: this IRTK process causes the nifti headers to become qform rather 
%  than sform!!! Bit annoying...


%% studypath
sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Synthetic_Flow_Phantom';
cd(sp)


%% bSSFP dataset

gradmomDir = 'bSSFP_6pipes/data/';
dataDirNoisy = 'bSSFP_6pipes/data_noisy_snr89/';

stackNames = {...
    'stack_rx0_ry0_rz90_tx10_ty5_tz0_slices100',...
    'stack_rx0_ry90_rz0_tx10_ty5_tz0_slices100',...
    'stack_rx90_ry0_rz0_tx10_ty5_tz0_slices100',...
    'stack_rx45_ry0_rz0_tx10_ty5_tz0_slices100',...
    'stack_rx30_ry60_rz0_tx10_ty5_tz0_slices100',...
    'stack_rx50_ry0_rz0_tx10_ty5_tz0_slices100' ...
    };


%% load stacks and gradient moments
for ii = 1:numel(stackNames)
    tempNii = load_untouch_nii([dataDirNoisy stackNames{ii} '_affine.nii.gz']);
    STACKS(:,:,:,ii) = tempNii.img;
    
    G(ii,:) = load([gradmomDir stackNames{ii} '_scanner_grad_moments.txt']);
    clear tempNii
end


%% Rotate stack by 5 degrees and rotate grad_moments
%- start with rx45 stack
rx45 = load_untouch_nii([dataDirNoisy stackNames{4} '.nii.gz']);
rx45a = load_untouch_nii([dataDirNoisy stackNames{4} '_affine.nii.gz']);

rx50 = load_untouch_nii([dataDirNoisy stackNames{6} '.nii.gz']);
rx50a = load_untouch_nii([dataDirNoisy stackNames{6} '_affine.nii.gz']);


% % get affine:
% % -NB: x/y flipped in original synthetic flow phantom script
% % -might be important for Alena?
% AFForig45(2,:) = rx45.hdr.hist.srow_x; 
% AFForig45(1,:) = rx45.hdr.hist.srow_y;
% AFForig45(3,:) = rx45.hdr.hist.srow_z;
% AFForig45(4,:) = [0 0 0 1];
% 
% AFForig50(2,:) = rx50.hdr.hist.srow_x; 
% AFForig50(1,:) = rx50.hdr.hist.srow_y;
% AFForig50(3,:) = rx50.hdr.hist.srow_z;
% AFForig50(4,:) = [0 0 0 1];

% make rotation matrix:
xTheta = 5;
yTheta = 0;
zTheta = 0;

rx = rotxtar(deg2rad(xTheta));
ry = rotytar(deg2rad(yTheta));
rz = rotztar(deg2rad(zTheta));

rx(:,4) = 0; ry(:,4) = 0; rz(:,4) = 0;
rx(4,:) = [0 0 0 1]; ry(4,:) = [0 0 0 1]; rz(4,:) = [0 0 0 1];
R = rx*ry*rz;

% % new affine matrix
% AFFrot = R * AFForig45;
% 
% % apply rotation to 3D volume
% % AFFrot(:,4) = [0;0;0;1]; % necessary for imwarp
% tform = affine3d(R);
% rx45_rot.img = imwarp(rx45.img,tform);
% rx45a_rot.img = imwarp(rx45a.img,tform);


% view rx45
hFig1 = figure;
hAx1  = axes;
slice(double(rx45.img),size(rx45.img,1)/2,size(rx45.img,2)/2,size(rx45.img,3)/2);
grid on, shading interp, colormap gray; xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis'); title('rx45');

% view rx50
hFig2 = figure;
hAx2  = axes;
slice(double(rx50.img),size(rx50.img,1)/2,size(rx50.img,2)/2,size(rx50.img,3)/2);
grid on, shading interp, colormap gray; xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis'); title('rx50');

% rotate rx45 by 5 degrees in x
% ---NB: imrotate3 goes in opposite direction to my code, ie: +5 degree
% rotation in synthetic_phantom.m = -5 degree rotation in imrotate3
rx45_rot.img = imrotate3(rx45.img,-5,[0 1 0],'crop','FillValues',0);

hFig3 = figure;
hAx3  = axes;
slice(double(rx45_rot.img),size(rx45_rot.img,1)/2,size(rx45_rot.img,2)/2,size(rx45_rot.img,3)/2);
grid on, shading interp, colormap gray; xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis'); title('rx45 rot/trans');

% view subtraction
% implay_RR(rx50.img - rx45_rot.img);



%% Rotate grad_moments
Grx45 = G(4,:);
Grx50 = G(6,:);

Grot = R(1:3,1:3) * Grx45'; % should be equal Grx50

% Grx50
% Grot'






%% Vectorise phase stacks
P.stack1 = STACKS(:,:,:,1); P.stack1 = P.stack1(:);
P.stack2 = STACKS(:,:,:,2); P.stack2 = P.stack2(:);
P.stack3 = STACKS(:,:,:,3); P.stack3 = P.stack3(:);
P.stack4 = STACKS(:,:,:,4); P.stack4 = P.stack4(:);
P.stack5 = STACKS(:,:,:,5); P.stack5 = P.stack5(:);
P.stack6 = STACKS(:,:,:,6); P.stack6 = P.stack6(:);


%% Scaling
gamma = 42577; %Hz/mT
% gamma = 2 .* pi .* 42577; %Hz/mT
Mscaling = (1e-3).^2;  %ms^2.mT/m --- First Moment scaling into correct units


%% 5 stacks
%- x90 / y90 / z90 / rx45 / rx30-ry60
Pmat_orig = [P.stack1, P.stack2, P.stack3, P.stack4, P.stack5]';
G_world_orig = G([1,2,3,4,5],:);


% Solve simultaenous equations for each voxel
for ii = 1:length(Pmat_orig)
    Vworld_vec(:,ii) = gamma .* Mscaling .* G_world_orig\Pmat_orig(:,ii); %X = A\B => Vxyz*Vxyz_vec = Pmat
end

disp('Solved.');


% Un-vectorise and convert to velocity stacks
subV=ind2subV(size(STACKS(:,:,:,1)),1:numel(P.stack1));

for ii = 1:length(subV)
    i = subV(ii,1);
    j = subV(ii,2);
    k = subV(ii,3);
    
    %world coords
    VEL_orig5.rl(i,j,k) = Vworld_vec(1,ii); % metres/sec
    VEL_orig5.ap(i,j,k) = Vworld_vec(2,ii);
    VEL_orig5.fh(i,j,k) = Vworld_vec(3,ii);
end

VEL3D_orig5 = VEL_orig5.rl + VEL_orig5.ap + VEL_orig5.fh;

disp('VEL done');

% %%% View
% implay_RR([VEL_orig5.rl, VEL_orig5.ap, VEL_orig5.fh],'jet',[-1,1]);
% implay_RR(VEL3D_orig5, 'jet', [-1,1] );


% 3 stacks
%- same as above but 3 stacks
%- x90 / y90 /  rx45
%%
% Pmat_orig = [P.stack1, P.stack2, P.stack4]';
% G_world_orig = G([1,2,4],:);
% 
% 
% % Solve simultaenous equations for each voxel
% for ii = 1:length(Pmat_orig)
%     Vworld_vec(:,ii) = gamma .* Mscaling .* G_world_orig\Pmat_orig(:,ii); %X = A\B => Vxyz*Vxyz_vec = Pmat
% end
% 
% disp('Solved.');
% 
% 
% % Un-vectorise and convert to velocity stacks
% subV=ind2subV(size(STACKS(:,:,:,1)),1:numel(P.stack1));
% 
% for ii = 1:length(subV)
%     i = subV(ii,1);
%     j = subV(ii,2);
%     k = subV(ii,3);
%     
%     %world coords
%     VEL_orig3.rl(i,j,k) = Vworld_vec(1,ii); % metres/sec
%     VEL_orig3.ap(i,j,k) = Vworld_vec(2,ii);
%     VEL_orig3.fh(i,j,k) = Vworld_vec(3,ii);
% end
% 
% VEL3D_orig3 = VEL_orig3.rl + VEL_orig3.ap + VEL_orig3.fh;
% 
% disp('VEL done');
% 
% % %%% View
% % implay_RR([VEL_orig3.rl, VEL_orig3.ap, VEL_orig3.fh],'jet',[-1,1]);
% % implay_RR(VEL3D_orig3, 'jet', [-1,1] );



%% 5 stacks
%- replace rx45 with rx50, plus Grx45 rotated by 5 degrees
%- x90 / y90 / z90 / rx50, but using Grot instead / rx30-ry60
Pmat_rot = [P.stack1, P.stack2, P.stack3, P.stack6, P.stack5]';
G_world_rot = G([1,2,3,4,5],:);
G_world_rot(4,:) = Grot'; %replace Grx45 with Grot

% Solve simultaenous equations for each voxel
for ii = 1:length(Pmat_rot)
    Vworld_vec(:,ii) = gamma .* Mscaling .* G_world_rot\Pmat_rot(:,ii); %X = A\B => Vxyz*Vxyz_vec = Pmat
end

disp('Solved.');


% Un-vectorise and convert to velocity stacks
subV=ind2subV(size(STACKS(:,:,:,1)),1:numel(P.stack1));

for ii = 1:length(subV)
    i = subV(ii,1);
    j = subV(ii,2);
    k = subV(ii,3);
    
    %world coords
    VEL_rot.rl(i,j,k) = Vworld_vec(1,ii); % metres/sec
    VEL_rot.ap(i,j,k) = Vworld_vec(2,ii);
    VEL_rot.fh(i,j,k) = Vworld_vec(3,ii);
end

VEL3D_rot = VEL_rot.rl + VEL_rot.ap + VEL_rot.fh;

disp('VEL done');

% %%% View
% implay_RR([VEL_rot.rl, VEL_rot.ap, VEL_rot.fh],'jet',[-1,1]);
% implay_RR(VEL3D_rot, 'jet', [-1,1] );


%% 5 stacks
%- replace rx45 with rx50, plus Grx45 rotated by 5 degrees
%- x90 / y90 / z90 / rx50, but using Grot instead / rx30-ry60
Pmat_rot = [P.stack1, P.stack2, P.stack3, P.stack6, P.stack5]';
G_world_rot = G([1,2,3,4,5],:);
G_world_rot(4,:) = Grot'; %replace Grx45 with Grot

% Solve simultaenous equations for each voxel
for ii = 1:length(Pmat_rot)
    Vworld_vec(:,ii) = gamma .* Mscaling .* G_world_rot\Pmat_rot(:,ii); %X = A\B => Vxyz*Vxyz_vec = Pmat
end

disp('Solved.');


% Un-vectorise and convert to velocity stacks
subV=ind2subV(size(STACKS(:,:,:,1)),1:numel(P.stack1));

for ii = 1:length(subV)
    i = subV(ii,1);
    j = subV(ii,2);
    k = subV(ii,3);
    
    %world coords
    VEL_rot.rl(i,j,k) = Vworld_vec(1,ii); % metres/sec
    VEL_rot.ap(i,j,k) = Vworld_vec(2,ii);
    VEL_rot.fh(i,j,k) = Vworld_vec(3,ii);
end

VEL3D_rot = VEL_rot.rl + VEL_rot.ap + VEL_rot.fh;

disp('VEL done');

% %%% View
% implay_RR([VEL_rot.rl, VEL_rot.ap, VEL_rot.fh],'jet',[-1,1]);
% implay_RR(VEL3D_rot, 'jet', [-1,1] );









%% 5 stacks - current SVRTK method
%--- rotates stacks without updating phase
Pmat_svrtk = [P.stack1, P.stack2, P.stack3, P.stack4, P.stack5]';
G_world_svrtk = G([1,2,3,6,5],:);

%%% 5 degree rotation of Grx45
G_world_svrtk(4,:) = [0    2.0154  -14.6799];

% Solve simultaenous equations for each voxel
for ii = 1:length(Pmat_svrtk)
    Vworld_vec(:,ii) = gamma .* Mscaling .* G_world_svrtk\Pmat_svrtk(:,ii); %X = A\B => Vxyz*Vxyz_vec = Pmat
end

disp('Solved.');


% Un-vectorise and convert to velocity stacks
subV=ind2subV(size(STACKS(:,:,:,1)),1:numel(P.stack1));

for ii = 1:length(subV)
    i = subV(ii,1);
    j = subV(ii,2);
    k = subV(ii,3);
    
    %world coords
    VEL_svrtk5.rl(i,j,k) = Vworld_vec(1,ii); % metres/sec
    VEL_svrtk5.ap(i,j,k) = Vworld_vec(2,ii);
    VEL_svrtk5.fh(i,j,k) = Vworld_vec(3,ii);
end

VEL3D_svrtk5 = VEL_svrtk5.rl + VEL_svrtk5.ap + VEL_svrtk5.fh;

disp('VEL done');



%% 5 stacks - current SVRTK method
%--- rotates stacks without updating phase
Pmat_svrtk = [P.stack1, P.stack2, P.stack3, P.stack4, P.stack5]';
G_world_svrtk = G([1,2,3,6,5],:);

%%% 15 degree rotation of Grx45
G_world_svrtk(4,:) = [0    4.5339  -14.1069];

% Solve simultaenous equations for each voxel
for ii = 1:length(Pmat_svrtk)
    Vworld_vec(:,ii) = gamma .* Mscaling .* G_world_svrtk\Pmat_svrtk(:,ii); %X = A\B => Vxyz*Vxyz_vec = Pmat
end

disp('Solved.');


% Un-vectorise and convert to velocity stacks
subV=ind2subV(size(STACKS(:,:,:,1)),1:numel(P.stack1));

for ii = 1:length(subV)
    i = subV(ii,1);
    j = subV(ii,2);
    k = subV(ii,3);
    
    %world coords
    VEL_svrtk15.rl(i,j,k) = Vworld_vec(1,ii); % metres/sec
    VEL_svrtk15.ap(i,j,k) = Vworld_vec(2,ii);
    VEL_svrtk15.fh(i,j,k) = Vworld_vec(3,ii);
end

VEL3D_svrtk15 = VEL_svrtk15.rl + VEL_svrtk15.ap + VEL_svrtk15.fh;

disp('VEL done');



%% 5 stacks - current SVRTK method
%--- rotates stacks without updating phase
Pmat_svrtk = [P.stack1, P.stack2, P.stack3, P.stack4, P.stack5]';
G_world_svrtk = G([1,2,3,6,5],:);

%%% 25 degree rotation of Grx45
G_world_svrtk(4,:) = [0    6.9147  -13.1053];

% Solve simultaenous equations for each voxel
for ii = 1:length(Pmat_svrtk)
    Vworld_vec(:,ii) = gamma .* Mscaling .* G_world_svrtk\Pmat_svrtk(:,ii); %X = A\B => Vxyz*Vxyz_vec = Pmat
end

disp('Solved.');


% Un-vectorise and convert to velocity stacks
subV=ind2subV(size(STACKS(:,:,:,1)),1:numel(P.stack1));

for ii = 1:length(subV)
    i = subV(ii,1);
    j = subV(ii,2);
    k = subV(ii,3);
    
    %world coords
    VEL_svrtk25.rl(i,j,k) = Vworld_vec(1,ii); % metres/sec
    VEL_svrtk25.ap(i,j,k) = Vworld_vec(2,ii);
    VEL_svrtk25.fh(i,j,k) = Vworld_vec(3,ii);
end

VEL3D_svrtk25 = VEL_svrtk25.rl + VEL_svrtk25.ap + VEL_svrtk25.fh;

disp('VEL done');




%% View
% compare correct vs. different levels of gradient moment error
implay_RR([VEL_orig5.rl, VEL_orig5.ap, VEL_orig5.fh;
           VEL_svrtk5.rl, VEL_svrtk5.ap, VEL_svrtk5.fh
           VEL_svrtk15.rl, VEL_svrtk15.ap, VEL_svrtk15.fh
           VEL_svrtk25.rl, VEL_svrtk25.ap, VEL_svrtk25.fh],'jet',[-1,1]);

% view difference
implay_RR([VEL_orig5.rl-VEL_svrtk5.rl,...
           VEL_orig5.ap-VEL_svrtk5.ap,...
           VEL_orig5.fh-VEL_svrtk5.fh,...
           VEL3D_orig5 - VEL3D_svrtk5;
           VEL_orig5.rl-VEL_rot.rl,...
           VEL_orig5.ap-VEL_rot.ap,...
           VEL_orig5.fh-VEL_rot.fh,...
           VEL3D_orig5 - VEL3D_rot;
           VEL_orig5.rl-VEL_svrtk15.rl,...
           VEL_orig5.ap-VEL_svrtk15.ap,...
           VEL_orig5.fh-VEL_svrtk15.fh,...
           VEL3D_orig5 - VEL3D_svrtk15;
           VEL_orig5.rl-VEL_svrtk25.rl,...
           VEL_orig5.ap-VEL_svrtk25.ap,...
           VEL_orig5.fh-VEL_svrtk25.fh,...
           VEL3D_orig5 - VEL3D_svrtk25;],'jet',[-0.1,0.1]);

% view summation
implay_RR([VEL3D_orig5, VEL3D_svrtk3], 'jet', [-1,1] );
implay_RR(VEL3D_orig5 - VEL3D_svrtk3, 'jet', [-0.05,0.05] );


