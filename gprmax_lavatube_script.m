%% Modelling Lava Tubes in gprMax
% Author: Stephanie Bringeland
% The objective of this script is to take lava tube parameters and write 
% .in files for gprMax 
% The user will input individual or ranges of parameters and the script
% will write the corresponding number of .in files

% Issues:
% Not sure what relative permeability and magnetic loss refer to. Based on
% examples from gprMax it seems that non-magnetic materials have values of
% 1 and 0 respectively
% 
% Because of the voxelization function, the depth/horizontal and vertical scaling 
% factors no longer cycle through the for loop. Working on making this possible.

run_files = 1; % If you want to run the files through the OS, set this = 1.
% Beware this will take time if there are many files to run.
clear all
close all
%% Step 1 - Input Parameters
% Essential model elements 
domain = [10 10 10]; % model size in metres
dx_dy_dz = [0.002 0.002 0.002]; % model discretization in metres
time_window = 0.000000003; % total required simulated time in seconds

% General elements
messages ='y'; % y will allow the screen to display information while the 
% file is running, 'n' will suppress this
output_dir = "C:\Users\steph\Documents\Work\2020\ESA\matlab"; % output file 
% directory
title = "Lava Tube GPR Model";

% Material parameters - repeat for each material in the model
% Material 1
material1 = "lunarbasalt"; % material identifier
rel_permittivity1 = 10; % relative dielectric permittivity of material
%rel_permittivity1 = [7 8 9 10];
cond1 = 0.00000001; % conductivity in Siemens/m
rel_permeability1 = 1; % relative permeability 
mag_loss1 = 0; % magnetic loss in Ohms/metre (?)
% Material 2
material2 = "emptyspace";
rel_permittivity2 = 1;
cond2 = 0;
rel_permeability2 = 1;
mag_loss2 = 0; % need to check on these last 2 parameters for a vacuum, 
% not sure if these values make sense

% Source/Receiver
wave_type = "ricker"; % waveform type (source)
max_amp_scale = 1; % scaling of the maximum amplitude of the waveform 
centre_fq = [400000000]; % centre frequency of the 
% waveform in Hertz
wave_name = "source_wave"; % identifier
rx_pos = [0 8 0]; % receiver starting position
src_steps = [0.002 0 0]; % how much to move the source by each time (metres)
rx_steps = [0.002 0 0]; % how much to move the receiver by each time (metres)

%% Step 2 - Object Construction
% I'm initially going to write this for a simple buried cylinder in the 
% centre of the model, but eventually hope to expand it so more complex 
% shapes can be created. The cylinder is default horizontal
simple_shape = 0; % define the shape of a simple buried cylinder, no voxels
complex_shape = 1; % import coordinates to create a complex voxel model
cyl_mat = "emptyspace"; % tFhe material the void is made of

depth = 5; % depth of the centrepoint in metres
%depth = [5 10 15 20 25 50];
radius = 64; % radius of the cylinder in metres
%radius = [1 2 4 8 16 32 64 128];

File_directory = 'C:\Users\steph\Documents\Work\2020\ESA\matlab\galapagos-lava-tube.stl';
Model_name = 'v1.0_steph_test';

% The file either is an csv-file (use importdata) or a .stl file (use stlread)

% 'Cave_info' consists of faces and vertices
% 'stlread' is a function from the MatLab file exchange, see https://nl.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader
Cave_info = stlread(File_directory);
Cavity_coordinates = Cave_info.Points;

if complex_shape == true
    Voxel_filename = 'voxel_text.txt';
end

% If the number of points is much more than you need, you can reduce the
% number with this parameter:
Reduce_factor = 1;

% It takes some time and right now we don't use it, so it is commented out
% Cavity_coordinates = Reduce_coordinates(Cavity_coordinates, Reduce_factor);

% What this does; make the cavity Sizing_factor times larger, or stretch it
% by a factor Stretching_factor_h horizontally.
% Vertical stretching is a bit trickier, or shape-preserving stretching
% Alternative; stretch it until the average radius/length is .../... meter
Sizing_factor = 1;
Stretching_factor_h = 1;

Cavity_coordinates = Sizing_factor*Cavity_coordinates;
Cavity_coordinates = Stretching_factor_h*Cavity_coordinates;

% A grid will be made, and here you can choose the spacing between points
% (same in all directions)
% It is multiplied with the sizing factor for the cavity automatically
Grid_spacing = 5;
Grid_spacing = Grid_spacing*Sizing_factor;
Rock_density = 2.5; % density in g/cm³
Cavity_density = 0; % density in g/cm³
Model_density = abs(Rock_density - Cavity_density);
Model_name = [Model_name]; % The name you want your model to have
Edging_factor = 0.50; % How much additional spacing you want around the cavity, in order to prevent edge effects to play a dominant role. 
% This can be set 0 without too many problems, but from about 0.25 on you
% probably won't see the edge effects anymore. This does not affect the
% performance anymore, and just places extra stations
Cavity_depth = 10; % How much space filled with rock there is between the top of the cavity and the stations (on the ground)

if Grid_spacing <=0
    error('Error: grid spacing must be positive')
end
if ( (max(Cavity_coordinates(:,1))-min(Cavity_coordinates(:,1)))/Grid_spacing * (max(Cavity_coordinates(:,2))-min(Cavity_coordinates(:,2)))/Grid_spacing * (max(Cavity_coordinates(:,3))-min(Cavity_coordinates(:,3)))/Grid_spacing > 10000000) 
    warning('Warning: many voxels present (> 10 000 000), modelling process can take a while')
elseif  ( (max(Cavity_coordinates(:,1))-min(Cavity_coordinates(:,1)))/Grid_spacing * (max(Cavity_coordinates(:,2))-min(Cavity_coordinates(:,2)))/Grid_spacing * (max(Cavity_coordinates(:,3))-min(Cavity_coordinates(:,3)))/Grid_spacing > 30000000) 
    error('Error: too many voxels present (> 30 000 000). Reduce by adjusting the grid spacing')
end

if( (max(Cavity_coordinates(:,1))-min(Cavity_coordinates(:,1)))/Grid_spacing < 1 || (max(Cavity_coordinates(:,2))-min(Cavity_coordinates(:,2)))/Grid_spacing < 1 || (max(Cavity_coordinates(:,3))-min(Cavity_coordinates(:,3)))/Grid_spacing < 1)
    error('Error: grid spacing too high; zero voxels detected in one of the directions')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gridded_coords = {};
Sizing_addition = [];
Model_limits = [];
Cavity_limits = [];
for a=1:3 % In this loop (3 or 3x2 coordinates) the grid is created and the limits for the cavity and the model (with edging factor) are calculated
    Gridded_coords(a) = {round(min(Cavity_coordinates(:,a))):Grid_spacing:(round(max(Cavity_coordinates(:,a)))-Grid_spacing)};
    Sizing_addition(a) = round((Gridded_coords{1,a}(end)-Gridded_coords{1,a}(1))*Edging_factor/Grid_spacing)*Grid_spacing;
    Cavity_limits(2*a-1) = round(min(Cavity_coordinates(:,a)));
    Cavity_limits(2*a) = round(max(Cavity_coordinates(:,a)));
    if a<3
        Model_limits(2*a-1) = round(min(Cavity_coordinates(:,a)))-Sizing_addition(a);
        Model_limits(2*a) = round(max(Cavity_coordinates(:,a)))+Sizing_addition(a);
    else
        Model_limits(2*a-1) = round(min(Cavity_coordinates(:,a)))-Cavity_depth;
        Model_limits(2*a) = round(max(Cavity_coordinates(:,a)))+Cavity_depth;
    end
end

% Shift everything lower/higher with the maximum z-coordinate of the model,
% such that the model ends at z=0, the cavity begins at z=0-Cavity_depth
% and all limits are changed accordingly
Gridded_coords{1,3} = Gridded_coords{1,3} - Model_limits(6);
Cavity_coordinates(:,3) = Cavity_coordinates(:,3) - Model_limits(6);
Cavity_limits(5) = Cavity_limits(5) - Model_limits(6);
Cavity_limits(6) = Cavity_limits(6) - Model_limits(6);
Model_limits(5) = Model_limits(5) - Model_limits(6);
Model_limits(6) = Model_limits(6) - Model_limits(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is still to create a grid, but in the format of an n x n x n
% meshgrid. This makes it easier to create 'logicals', which speed uo the
% process later on

% Not the addition of half a grid spacing! This is since voxel coordinates
% actually signifiy the locations of the centers of voxels.
[X,Y,Z] = meshgrid(Gridded_coords{1,1}+Grid_spacing/2,Gridded_coords{1,2}+Grid_spacing/2,Gridded_coords{1,3}+Grid_spacing/2);
Logical = true(length(Gridded_coords{1,2}),length(Gridded_coords{1,1}),length(Gridded_coords{1,3}));
% The Grid_matrix_3n is an nx3 matrix consisting of every coordinate
% considered in the model
Grid_matrix_3n = [X(Logical),Y(Logical),Z(Logical)];

% Now that the grid has been created and is in the right format and the limits
% are known, the voxel models can be made

File_directory = [Model_name,'.txt']; % The model name defined earlier is added
Cavity_bool = true; % This is to use the following function twice, instead of making separate functions for it
[inside_box,outside_box] = Voxel_Creation(material1,material2,Cavity_coordinates,Grid_matrix_3n,Grid_spacing,Voxel_filename);


%% Step 3 - Write model files
% Not the most efficient way to do this but makes it so that a user can
% input any number of values for any of the variables that can be changed
% and the loop will create an .in file for each possible combination of
% parameters. The .in files created can be opened as text files to view the
% contents/syntax, and the documentation for the commands is at 
% http://docs.gprmax.com/en/latest/input.html

%inputfilenames = strings(1,10000); 
i = 1;
for a = 1:length(rel_permittivity1)
    for b = 1:length(rel_permittivity2)
        for c = 1:length(cond1)
            for d = 1:length(cond2)
                for e = 1:length(depth)
                    for f = 1:length(radius)
                        for g = 1:length(max_amp_scale)
                            for h = 1:length(centre_fq)
                                file_name = sprintf('gprmax_in\\gprMax_modeltest_fq%i_pv1%i_d%i_r%i.in',centre_fq(h),rel_permittivity1(a),depth(e),radius(f));
                                fid = fopen(file_name,'w');
                              
                                % The following line is SUPER long but when
                                % I tried to break it up it wouldn't work
                                % for some reason
                                fprintf(fid, '#title: %s \r\n#domain: %.4f %.4f %.4f \r\n#dx_dy_dz: %.4f %.4f %.4f \r\n#time_window: %.15f \r\n#material: %.3f %.3f %.3f %.3f %s \r\n#material: %.3f %.3f %.3f %.3f %s \r\n#waveform: %s %f %f %s \r\n#rx: %.4f %.4f %.4f\r\n#rx_steps: %.4f %.4f %.4f\r\n#src_steps: %.4f %.4f %.4f\r\n#box: %.4f %.4f %.4f %.4f %.4f %.4f %s\r\n#geometry_view: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f geometryview%s n',...
                                    title,domain(1),domain(2),domain(3),dx_dy_dz(1),dx_dy_dz(2),dx_dy_dz(3),...
                                    time_window,rel_permittivity1(a),cond1(c),rel_permeability1,mag_loss1,material1,...
                                    rel_permittivity2(b),cond2(d),rel_permeability2,mag_loss2,material2,...
                                    wave_type,max_amp_scale(g),centre_fq(h),wave_name,rx_pos(1),rx_pos(2),rx_pos(3),...
                                    rx_steps(1),rx_steps(2),rx_steps(3),src_steps(1),src_steps(2),src_steps(3),...
                                    0,0,0,domain(1),rx_pos(2),domain(3),material1,0,0,0,domain(1),...
                                    domain(2),domain(3),dx_dy_dz(1),dx_dy_dz(2),dx_dy_dz(3),file_name);
                                % if statement that prints either the voxels created by the voxel function or
                                % the cylinder defined by the depth and radius specified by the user
                                % (the order of the input file only matters in the sense that the lattermost object
                                % will overwrite all previous objects, so it doesn't matter that this is at the very end
                                if complex_shape == 1
                                    for k = (1:length(inside_box(:,1)))
                                        fprintf(fid,'#box %f %f %f %f %f %f %s\r\n',inside_box(k,1),inside_box(k,2),inside_box(k,3),inside_box(k,4),inside_box(k,5),inside_box(k,6),material2);
                                    end

                                    for k = (1:length(outside_box(:,1)))
                                        fprintf(fid,'#box %f %f %f %f %f %f %s\r\n',outside_box(k,1),outside_box(k,2),outside_box(k,3),outside_box(k,4),outside_box(k,5),outside_box(k,6),material1);
                                    end
                                elseif simple_shape == 1
                                    fprintf(fid,'#cylinder: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %s \r\n',...
                                    0.5*domain(1),rx_pos(2)-depth(e),0,...
                                    0.5*domain(1),rx_pos(2)-depth(e),domain(3),radius(f),material2);
                                end
                                fclose(fid);
                                inputfilenames(i) = string(file_name);
                                i = i+1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Step 4 - Run files through OS
% If the program creates a lot of files, beware this will likely take a
% while. One file took 5.6 seconds.
if run_files = true;
cmd1 = 'conda activate gprMax';

for k = 1:length(inputfilenames)
    [path(k),inputname(k),ext(k)] = fileparts(inputfilenames(k));
end

for j = 1:length(inputfilenames)
    cmd2 = sprintf('python -m gprMax gprmax_in\\%s -n 60',inputfilenames(j));
    cmd3 = sprintf('python -m tools.outputfiles_merge gprmax_in\\%s',inputname(j));
    command = sprintf('%s & %s & %s',cmd1,cmd2,cmd3);
    status = system(command);
end
end  
%% Function adapted from Frank's IGMAS voxelization function

 function [inside_box,outside_box] = Voxel_Creation(Rock_material,Void_material,Cavity_coordinates,Grid_matrix_3n,Grid_spacing,File_name)
    % The standard MatLab functionality of making 'alphaShapes' is used
    % here. If needed, they can be easily visualised with
    % 'plot(Cavity_Alphashape)'
    Cavity_Alphashape = alphaShape(Cavity_coordinates);
    Cavity_Alphashape.Alpha = Cavity_Alphashape.Alpha*2; % Some security which might has to be improved later on
    % The reason for the use of alphashapes (earlier we already had the
    % full model) is the built-in efficient way of checking whether
    % points are in the interior. This is done for all points together
    % and is much faster than alternatives
    In_bool = inShape(Cavity_Alphashape,Grid_matrix_3n(:,1),Grid_matrix_3n(:,2),Grid_matrix_3n(:,3));
    
    % A way to sort the inside and outside points and to add the density as
    % a fourth column
    Inside_points = [Grid_matrix_3n(In_bool,1),Grid_matrix_3n(In_bool,2),Grid_matrix_3n(In_bool,3)];
    %Inside_density = Model_density*ones(length(Inside_points(:,1)),1);
    
    Void_material = repelem(Void_material,length(Inside_points(:,1))).';
    Inside_points = struct('Inside_points',Inside_points,'Void_material',Void_material);

    Outside_points = [Grid_matrix_3n(~In_bool,1),Grid_matrix_3n(~In_bool,2),Grid_matrix_3n(~In_bool,3)];
    %Outside_density = zeros(length(Outside_points(:,1)),1);
    
    Rock_material = repelem(Rock_material,length(Outside_points(:,1))).';
    Outside_points = struct('Outside_points',Outside_points,'Rock_material',Rock_material);
    
    % 'Export' is the voxelised model and consists of interior and exterior
    % points with their respective densities.
    %Export = [Inside_points;Outside_points];
    
    if(length(Inside_points(:,1)) < 1)
        error('Error: no inside points detected. Decrease grid spacing, or increase alpha of alphaShape')
    end
    % This sorting is to have both files sorted the same way. It is useful
    % for checking, but should not affect performance
    %[~,idx] = sort(Export(:,3)); % sort just the first column
    %Export = Export(idx,:); 
    %[~,idx] = sort(Export(:,2)); % sort just the first column
    %Export = Export(idx,:); 
    %[~,idx] = sort(Export(:,1)); % sort just the first column
    %Export = Export(idx,:); 

    % The actual process of writing the voxel file
    %fileID = fopen(File_name,'w');
    %sprintf('x\t y\t z\t cellValue\n');
    for i = 1:length(Inside_points.Inside_points(:,1))
      in_lowleft_x(i) = (-0.5*Grid_spacing)+Inside_points.Inside_points(i,1);
      in_lowleft_y(i) = (-0.5*Grid_spacing)+Inside_points.Inside_points(i,2);
      in_lowleft_z(i) = (-0.5*Grid_spacing)+Inside_points.Inside_points(i,3);
      in_upright_x(i) = (0.5*Grid_spacing)+Inside_points.Inside_points(i,1);
      in_upright_y(i) = (0.5*Grid_spacing)+Inside_points.Inside_points(i,2);
      in_upright_z(i) = (0.5*Grid_spacing)+Inside_points.Inside_points(i,3);
    end
    for j = 1:length(Outside_points.Outside_points(:,1))
      
      out_lowleft_x(j) = (-0.5*Grid_spacing)+Outside_points.Outside_points(j,1);
      out_lowleft_y(j) = (-0.5*Grid_spacing)+Outside_points.Outside_points(j,2);
      out_lowleft_z(j) = (-0.5*Grid_spacing)+Outside_points.Outside_points(j,3);
      out_upright_x(j) = (0.5*Grid_spacing)+Outside_points.Outside_points(j,1);
      out_upright_y(j) = (0.5*Grid_spacing)+Outside_points.Outside_points(j,2);
      out_upright_z(j) = (0.5*Grid_spacing)+Outside_points.Outside_points(j,3);
      %inside_vox = sprintf('%d %d %d %s\n', Export(i,1),Export(i,2),Export(i,3),Export(i,4));
    end
    
inside_box = [(in_lowleft_x).',(in_lowleft_y).',(in_lowleft_z).',(in_upright_x).',(in_upright_y).',(in_upright_z).'];
outside_box = [(out_lowleft_x).',(out_lowleft_y).',(out_lowleft_z).',(out_upright_x).',(out_upright_y).',(out_upright_z).'];

end   


                            
                                


