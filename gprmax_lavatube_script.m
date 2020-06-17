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

% Object Construction
% I'm initially going to write this for a simple buried cylinder in the 
% centre of the model, but eventually hope to expand it so more complex 
% shapes can be created. The cylinder is default horizontal
depth = 5; % depth of the centrepoint in metres
radius = 50; % radius of the cylinder in metres
cyl_mat = "emptyspace"; % the material the cylinder is made of

% Source/Receiver
wave_type = "ricker"; % waveform type (source)
max_amp_scale = 1; % scaling of the maximum amplitude of the waveform 
centre_fq = [400000000, 500000000, 600000000]; % centre frequency of the 
% waveform in Hertz
wave_name = "source_wave"; % identifier
rx_pos = [0 8 0]; % receiver starting position
src_steps = [0.002 0 0]; % how much to move the source by each time (metres)
rx_steps = [0.002 0 0]; % how much to move the receiver by each time (metres)

%% Step 2 - Write model files
% Not the most efficient way to do this but makes it so that a user can
% input any number of values for any of the variables that can be changed
% and the loop will create an .in file for each possible combination of
% parameters. The .in files created can be opened as text files to view the
% contents/syntax, and the documentation for the commands is at 
% http://docs.gprmax.com/en/latest/input.html

inputfilenames = strings(1,10000); 
i = 1;
for a = 1:length(rel_permittivity1)
    for b = 1:length(rel_permittivity2)
        for c = 1:length(cond1)
            for d = 1:length(cond2)
                for e = 1:length(depth)
                    for f = 1:length(radius)
                        for g = 1:length(max_amp_scale)
                            for h = 1:length(centre_fq)
                                file_name = sprintf('gprMax_modeltest_fq%i.in',centre_fq(h));
                                fid = fopen(file_name,'w');
                                text
                                fprintf(fid, '#title: %s \r\n#domain: %.4f %.4f %.4f \r\n#dx_dy_dz: %.4f %.4f %.4f \r\n#time_window: %.15f \r\n#material: %.3f %.3f %.3f %.3f %s \r\n#material: %.3f %.3f %.3f %.3f %s \r\n#waveform: %s %f %f %s \r\n#rx: %.4f %.4f %.4f\r\n#rx_steps: %.4f %.4f %.4f\r\n#src_steps: %.4f %.4f %.4f\r\n#box: %.4f %.4f %.4f %.4f %.4f %.4f %s \r\n#cylinder: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %s \r\n#geometry_view: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f geometryview%s n',...
                                    title,domain(1),domain(2),domain(3),dx_dy_dz(1),dx_dy_dz(2),dx_dy_dz(3),...
                                    time_window,rel_permittivity1(a),cond1(c),rel_permeability1,mag_loss1,material1,...
                                    rel_permittivity2(b),cond2(d),rel_permeability2,mag_loss2,material2,...
                                    wave_type,max_amp_scale(g),centre_fq(h),wave_name,rx_pos(1),rx_pos(2),rx_pos(3),...
                                    rx_steps(1),rx_steps(2),rx_steps(3),src_steps(1),src_steps(2),src_steps(3),...
                                    0,0,0,domain(1),rx_pos(2),domain(3),material1,0.5*domain(1),rx_pos(2)-depth(e),0,...
                                    0.5*domain(1),rx_pos(2)-depth(e),domain(3),radius(f),material2,0,0,0,domain(1),...
                                    domain(2),domain(3),dx_dy_dz(1),dx_dy_dz(2),dx_dy_dz(3),file_name);
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

%% Step 3 - Run files through OS

                            
                                


