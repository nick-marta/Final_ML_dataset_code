clc; clear; close all;
format short
%load DeepMIMO_FURI_Demo.mat % For DeepMIMO Comparison
for y = [1:2]
    for x = [1:5]
BW = 100e6; % spacing for DeepMIMO unit in Hz

%% Signal Parameters
freq = 2.4e9; % transmit frequency
num_subcarriers = 1024;

%% Math Constants
light_constant = 3e8; % speed of light constant

%% Grid Parameters
% top_left = [-150; -9; 1];
% top_right = [-150; 9; 1];
% bottom_left = [160; -9; 1];
% bottom_right = [160; 9; 1];
% num_row = 50;
% num_col = 10;
% num_users = num_col * num_row;
x_coord_list = [1:5];
x_coord_list_2 = [1:5];
y_coord_list = [1:2];

%0.5m
%x_coord_list = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10];
%x_coord_list_2 = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10];
%y_coord_list = [0.5,1,1.5,2,2.5];

%0.1m
%x_coord_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
%x_coord_list_2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
%y_coord_list = [0.1,0.2,0.3,0.4,0.5];

x_coord_list = [x_coord_list,x_coord_list_2];

y_coord_list = repelem(y_coord_list,5);

z_coord_list = 1;
z_coord_list = repelem(z_coord_list,10);

disp(x_coord_list)
disp(y_coord_list)
disp(z_coord_list)
%% Quick Grid Parameters
%second number is x axis, first is y axis
top_left = [-5; 5/2; 1]; %bottom left
top_right = [-5; 5/2; 1]; %top left
bottom_left = [10; 5/2; 1]; % top right
bottom_right = [10; 5/2; 1]; %bottom right
num_row = 3;
num_col = 1;
num_users = 11; %think it's the number of rx users + 1

%% Script Execution Parameters
load_environment = true; % do you want to load the stl environment
generate_user_coordinates = false; % do you want to generate user coordinates or load user coordinates
generate_tx_rx_objects = true;
is_perform_raytracing = true; % do you want to prefor the raytracing or load rays

%% Dataset Parameters
tx_users = [1];
rx_users = [1,3];

%% LOAD THE ENVIRONMNET
% choose to load the environment or to generate it
if load_environment == true
    % generate environmnet
    fprintf('... Generating the environment ...\n')
    TR = stlread("City_Buildings.stl"); % Read the STL file
    save city_buildings.mat TR
    viewer = siteviewer("SceneModel", TR); % View the STL environmnet
else
    % load environmnet
    fprintf('... Loading the environmnet ... \n')
    load city_buildings.mat
%     viewer = siteviewer("SceneModel", TR); % View the STL environmnet
end

%% GENERATE USER COORDINATES
% choose to load the user coordinates or to generate them
fprintf('... Generating user coordinates ...\n')
if generate_user_coordinates == true
    % generate user coordinates
    [x_axis_coordinates, y_axis_coordinates, z_axis_coordinates] = create_user_grid_pairs(top_left, top_right, bottom_left, num_row, num_col, num_users);
else 
    % load user coordinates
    fprintf('... Loading user coordinates ... \n')
    x_axis_coordinates = [10,-5];
    y_axis_coordinates = [2.5,2.5];
    z_axis_coordinates = [1,1];
    % TODO- Load user coordinates
end

%% GENERATE TRANSMITTER AND RECIEVER COMBINATIONS
fprintf('... Generating transmitter and reciever site combinations ...\n')
% choose to load the tx and rx objects or gernate them
if generate_tx_rx_objects == true
    % generate tx and rx site objects
    [tx,rx] = desinate_sites(x,y,x_coord_list, y_coord_list, z_coord_list, freq, light_constant);
    fprintf('... Saving transmitter and reciever coordinates ... \n')
    % TODO- Save transmitter and reciever coordinates
else
    % load the tx_rx_sites

    fprintf('... Loading transmitter and reciever coordinates ... \n')
    % TODO- Load transmitter and reciever coordinates
end

%% DEFINE PROPAGATION MODEL
% define propagation model
pm = propagationModel("raytracing", ... % Defining model as a ray tracing model
    "CoordinateSystem","cartesian", .... % Defining coordinate system as cartesian
    "Method","sbr",... % Defining method as sbr
    "MaxNumReflections",1, ... % Only record paths with certain number of reflections (0 indexed)
    "MaxNumDiffractions",1, ... % Only record paths with certain number of diffractions (0 indexed)
    "SurfaceMaterial","concrete"); % Define material of STL

%% PERFORM THE RAYTRACING
% choose to load the raytracing or perform the raytracing
fprintf('... Performing the raytracing ...\n')

if is_perform_raytracing == true
    % perform the raytracing
    rays = raytrace(tx,rx,pm); % Generates the rays from progation model 'pm', transmitter 'tx', and reciever 'rx'

else
    % load the raytracing
    load rays.mat
end

fprintf('... Raytracing done! ...\n')

%% Create Channel From Rays
% Pull desired rays
desired_rays = rays(tx_users,:);

% Create channels from desired rays
raytracing_channel_array = cell(1, length(desired_rays)); % Pre-allocate cell array
for i = 1:length(desired_rays)
    rays_temp = rays{1, i};
    if isempty(rays_temp) 
        continue;
    else
        tx_temp = tx(1); % Transmitter
        rx_temp = rx(i); % Receiver
        % Create a channel using the create_channel function
        raytracing_channel_temp = create_channel(rays_temp, tx_temp, rx_temp, BW);
        % Save the created channel to the cell array
        raytracing_channel_array{i} = raytracing_channel_temp;
    end
end
TxSig = zeros(num_subcarriers,1);
TxSig(1) = 1;

chan_dataset=zeros(num_users-1,1024); %make this genertic to match the nuner of users
for i = 2:num_users %make this generric
    if (~isempty(raytracing_channel_array{i-1})) %it would run into issues when tx and rx were in the same location, that is resolved with this
    chan_dataset(i-1,:)=fft(raytracing_channel_array{i-1}(TxSig)).';
    end
end

if (x == 1 && y == 1)
   writematrix(chan_dataset,'C:\Users\NickMarta\Documents\419_vs_code_projs\590_proj\week of 4_25\4_26_10_users.txt','Delimiter','tab')
else
    writematrix(chan_dataset,'C:\Users\NickMarta\Documents\419_vs_code_projs\590_proj\week of 4_25\4_26_10_users.txt','Delimiter','tab','Writemode','append')
end

    end
end
%% Perform Raytracing

function [x_axis_coordinates_temp, y_axis_coordinates_temp, z_axis_coordinates_temp] = create_user_grid_pairs(top_left, top_right, bottom_left, num_row, num_col, num_users)
    
    starting_row = top_left(1);
    ending_row = bottom_left(1);
    
    starting_col = top_left(2);
    ending_col = top_right(2);
   
    x_axis_temp = linspace(starting_row,ending_row,num_row);
    y_axis_temp = linspace(starting_col,ending_col, num_col);

    z_axis_coordinates_temp = ones(1,num_users);
    
    temp = ones(1,length(y_axis_temp));
    x_axis_coordinates_temp = kron(temp,x_axis_temp);
    
    temp = ones(1,length(x_axis_temp));
    y_axis_coordinates_temp = kron(temp,y_axis_temp);
    y_axis_coordinates_temp = sort(y_axis_coordinates_temp);

%     save user_coordinates.mat x_axis_coordinates_temp y_axis_coordinates_temp z_axis_coordinates_temp

end

function [tx, rx] = desinate_sites(x,y,x_coord_list, y_coord_list, z_coord_list, freq,c)
       
    tx = txsite("cartesian", ... % Define coordinate system
        "AntennaPosition",[x; y; 1], ... % Define transmitter antenna location
        "TransmitterFrequency",freq, ... % Define transmitter frequency
        "TransmitterPower",1/64); % In terms of Watts per sample (If I want 1 W with 10 sample each has to be 1/10 W)
    show(tx,"ShowAntennaHeight",false)
    
    cfgArray = arrayConfig("Size",[1 1],"ElementSpacing",(c/freq)/2);
    
    rx = rxsite("cartesian", ... % Define coordinate system
        "AntennaPosition",[x_coord_list; y_coord_list; z_coord_list], ... % Define transmitter antenna location
        "Antenna", cfgArray); %Define location of Recieving antenna
    show(rx,"ShowAntennaHeight",false)

end

function raytracing_channel = create_channel(rays, tx, rx, BW)

    raytracing_channel = comm.RayTracingChannel(rays,tx,rx);
    raytracing_channel.SampleRate = BW;
    raytracing_channel.ChannelFiltering = true; % filter an input signal
    raytracing_channel.ReceiverVirtualVelocity = [0; 0; 0]; % Velocity of reciever
    raytracing_channel.NormalizeChannelOutputs = false;
    raytracing_channel.NormalizeImpulseResponses = false;

end

% function output = perform_raytracing(channel)
% 
% end


%% OLD RAYTRACING LOOP
% if perform_raytracing == true
%     % perform the raytracing
%     rays = raytrace(tx_array,tx_array,pm); % Generates the rays from progation model 'pm', transmitter 'tx', and reciever 'rx'
% 
%     for i = 2:num_users
%         rays_temp = comm.Ray;
%         % rays_temp = raytrace(tx_array(i),rx_array(i,:),pm); % Generates the rays from progation model 'pm', transmitter 'tx', and reciever 'rx'
%         rays = [rays; rays_temp];
%     end
%     save rays.mat rays
% else
%     % load the raytracing
%     load rays.mat
% end
% fprintf('... Raytracing done! ...\n')

%% OLD SITE DESINATION (EXLUDE ONE USER AT A TIME E.G. 25X24)
%     x_axis_coordinates = zeros(num_users, length(x_axis_coordinates_temp)-1);
%     y_axis_coordinates = zeros(num_users, length(y_axis_coordinates_temp)-1);
%     z_axis_coordinates = zeros(num_users, length(z_axis_coordinates_temp)-1);
% 
%     for i = 1:num_users
%     % Exclude user i (pointer for excluded user)
%         
%         for k = 1:num_users
% 
%             if k ~= i
%             % Store x, y, z coordinates of recievers excluding user i
%                 if k < i
%                     x_axis_coordinates(i,k) = x_axis_coordinates_temp(k);
%                     y_axis_coordinates(i,k) = y_axis_coordinates_temp(k);
%                     z_axis_coordinates(i,k) = z_axis_coordinates_temp(k);
% 
%                 else
%                     x_axis_coordinates(i,k-1) = x_axis_coordinates_temp(k);
%                     y_axis_coordinates(i,k-1) = y_axis_coordinates_temp(k);
%                     z_axis_coordinates(i,k-1) = z_axis_coordinates_temp(k);
%                 end
%             end
%         end
%     end

%% OLD CREATE RAYTRACING CHANNEL
% rays = rays.'; % Stores rays
% 
% TODO- 64 sub carriers, 1 TX, 10 RX, channel data (10x64 complex),
%       NOT geographically close, different channels (classification),
%       generate dataset on the remote desktop, create shell script to
%       automate data sent to dropbox
%
% raytracing_channel = comm.RayTracingChannel(rays,tx_cell_array,rx_cell_array);
% raytracing_channel.SampleRate = BW;
% raytracing_channel.ChannelFiltering = true; % filter an input signal
% raytracing_channel.ReceiverVirtualVelocity = [0; 0; 0]; % Velocity of reciever
% raytracing_channel.NormalizeChannelOutputs = false;
% raytracing_channel.NormalizeImpulseResponses = false;
% 
% showProfile(raytracing_channel) % Plot delay spread, AoA, and AoD
% 
% %% Create transmitted signal
% tx_signal = zeros(1,num_ofdm);
% tx_signal(1,1) = ones(1,1);
% tx_signal = tx_signal(:);
% 
% y = raytracing_channel(tx_signal); % pass signal through channel
% h = abs(y);
% 
% %% Plot Raytracing Output
% figure("Name",'RtChan Output')
% stem(h)
% title('Channel Taps')
% xlabel('Samples')
% ylabel('Normalized Magnitude')

%% Pull Channel Parameters
% for user = 1:num_users
%     rays_temp = rays{user};
% 
%     if isempty(rays_tremp)
%         continue
%     else 
%         num_of_rays = length(rays_temp);
%         new_rays(1:len) = RayStructure; 
%     end
% end
