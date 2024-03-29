function [ synthPipe ] = synthflow_waterpipe( vel, volPx, varargin )
%SYNTHFLOW_WATERPIPE   create a synthetic pipe of constant flowing water
%
%   synthPipe = SYNTHFLOW_WATERPIPE( vel, volPx )
%
%   SYNTHFLOW_WATERPIPE( ..., 'param', val)
%
%   output:
%       synthPipe           synthetic water pipe            (m/s)
%
%   input:
%       vel                 peak velocity of water in pipe  (m/s)
%       volPx               3D volume dimensions            (xpixels,  ypixels,  zpixels)
%
%   optional parameter-value pairs:
%       pipeSize            size of pipe in pixels          (xsize,  ysize)
%       pipeDir             direction of pipe               ('x' / 'y' / 'z')
%       pipePos             position of pipe                (xpos,  ypos)
%
%   see also: synthflow_phantom, synthflow_acquire

% Tom Roberts (t.roberts@kcl.ac.uk)

%% Notes

% Function creates a truncated 2D Gaussian flow profile, then extrudes this
% into a 3D cylindrical pipe. By default, pipe is centred in middle of volume.

% TO DO: improve extrusion. x/y/z direction is a bit clunky.


%% Parse Input

default.pipeSize          = [round(0.5*volPx(1)/10) , round(0.5*volPx(1)/10)];
default.pipeDir           = 'z';
default.pipePos           = [0, 0];

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(  p, 'vel' );

addRequired(  p, 'volPx' );

add_param_fn( p, 'pipeSize', default.pipeSize, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {}, mfilename ) );

add_param_fn( p, 'pipeDir', default.pipeDir, ...
        @(x) validateattributes( x, {'char'}, ...
        {}, mfilename ) );

add_param_fn( p, 'pipePos', default.pipePos, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {}, mfilename ) );

parse( p, vel, volPx, varargin{:} );

pipeSize   = p.Results.pipeSize;
pipeDir    = p.Results.pipeDir;
pipePos    = p.Results.pipePos; 


%% Gaussian Flow Profile
gaussianProfile = customgauss( [volPx(1), volPx(2)], pipeSize(1), pipeSize(2), 0, 0, vel, pipePos );


%% Truncate Gaussian Profile at Edges of Pipe
outsidePipeIdx = find( abs(gaussianProfile) < 0.2 * abs(vel) ); % abs accounts for negative flows
gaussianProfile( outsidePipeIdx ) = 0;


%% Extrude Profile into 3D Pipe
if strcmp(pipeDir,'x')
    gaussianProfile(:,volPx(3)+1:end)=[];
    synthPipe = repmat( gaussianProfile, 1, 1, volPx(1) );
elseif strcmp(pipeDir,'y')
    gaussianProfile(volPx(3)+1:end,:)=[];
    synthPipe = repmat( gaussianProfile, 1, 1, volPx(2) );
elseif strcmp(pipeDir,'z')   
    synthPipe = repmat( gaussianProfile, 1, 1, volPx(3) );
end


end  % synthflow_waterpipe(...)