%% synthetic_flow_phantom_recon.m
%
% - reconstruct velocity volumes using synthetic flow phantom data


%% transform stacks into same grid using irtk
% Run IRTK commands in: stacks_to_common_grid_script.txt


%% studypath
sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Synthetic_Flow_Phantom';
cd(sp)


%% QFlow dataset

dataDir = 'QFlow_6pipes/data/';

stackNames = {...
    'stack_rx0_ry0_rz90_tx10_ty5_tz0_slices100',...
    'stack_rx0_ry90_rz0_tx10_ty5_tz0_slices100',...
    'stack_rx90_ry0_rz0_tx10_ty5_tz0_slices100',...
    };


%% bSSFP dataset

dataDir = 'bSSFP_6pipes/data/';

stackNames = {...
    'stack_rx0_ry0_rz90_tx10_ty5_tz0_slices100',...
    'stack_rx0_ry90_rz0_tx10_ty5_tz0_slices100',...
    'stack_rx90_ry0_rz0_tx10_ty5_tz0_slices100',...
    'stack_rx45_ry0_rz0_tx10_ty5_tz0_slices100',...
    'stack_rx30_ry60_rz0_tx10_ty5_tz0_slices100' ...
    };


%% load stacks and gradient moments
for ii = 1:numel(stackNames)
    tempNii = load_untouch_nii([dataDir stackNames{ii} '_affine.nii.gz']);
    STACKS(:,:,:,ii) = tempNii.img;
    
    G(ii,:) = load([dataDir stackNames{ii} '_grad_moments.txt']);
    clear tempNii
end



%% Vectorise phase stacks
P.stack1 = STACKS(:,:,:,1); P.stack1 = P.stack1(:);
P.stack2 = STACKS(:,:,:,2); P.stack2 = P.stack2(:);
P.stack3 = STACKS(:,:,:,3); P.stack3 = P.stack3(:);
P.stack4 = STACKS(:,:,:,4); P.stack4 = P.stack4(:);
P.stack5 = STACKS(:,:,:,5); P.stack5 = P.stack5(:);


%% Scaling
gamma = 42577; %Hz/mT
% gamma = 2 .* pi .* 42577; %Hz/mT
Mscaling = (1e-3).^2;  %ms^2.mT/m --- First Moment scaling into correct units


%% 3 orthogonal stacks
Pmat = [P.stack1, P.stack2, P.stack3]';
G_world = G([1,2,3],:);
% G_world(1,:) = G_world(1,:) * -1;
% G_world(2,:) = G_world(2,:) * -1;
% G_world(3,:) = G_world(3,:) * -1;

% G_world(:,1) = G_world(:,1) * -1;
% G_world(:,2) = G_world(:,2) * -1;
% G_world(:,3) = G_world(:,3) * -1;


%% 3 stacks: 1 orth, 2 oblique
% Pmat = [P.stack1, P.stack4, P.stack5]';   % 1 orthogonal + 2 oblique
% G_world = G([1,4,5],:);

%% 5 stacks
% P.stack4 = STACKS(:,:,:,4); P.stack4 = P.stack4(:);
% P.stack5 = STACKS(:,:,:,5); P.stack5 = P.stack5(:);
% Pmat = [P.stack1, P.stack2, P.stack3, P.stack4, P.stack5]';   % orthogonal stacks only + 2 oblique
% G_world = G; % NB: no need for Mscaling as grad_moments defined: .* (1e-3)^2; %msec^2.mT/m in original synthetic data  


%% Solve simultaenous equations for each voxel
for ii = 1:length(Pmat)
    Vworld_vec(:,ii) = gamma .* Mscaling .* G_world\Pmat(:,ii); %X = A\B => Vxyz*Vxyz_vec = Pmat
end

disp('Solved.');


%% Un-vectorise and convert to velocity stacks
subV=ind2subV(size(STACKS(:,:,:,1)),1:numel(P.stack1));

for ii = 1:length(subV)
    i = subV(ii,1);
    j = subV(ii,2);
    k = subV(ii,3);
    
    %world coords
    VEL.rl(i,j,k) = Vworld_vec(1,ii); % metres/sec
    VEL.ap(i,j,k) = Vworld_vec(2,ii);
    VEL.fh(i,j,k) = Vworld_vec(3,ii);
end

VEL3D = VEL.rl + VEL.ap + VEL.fh;

disp('VEL done');


%% View
implay_RR([VEL.rl, VEL.ap, VEL.fh],'jet',[-1,1]);
implay_RR(VEL3D, 'jet', [-1,1] );


%% load original synthetic phantom
Vpad = load_untouch_nii('VEL_VOLUME_srow4_equal0.nii.gz');
implay_RR(Vpad.img(:,:,250:350),'jet',[-1,1]);

% imtar(VEL3D(:,:,31)); colormap('jet');
% imtar(Vpad.img(:,:,250+31)); colormap('jet');




