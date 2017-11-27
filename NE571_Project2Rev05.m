%% NE 571: Project 1.
% William Gerding, CJ Oldham, Austin Stanford.
% Due: 2nd Oct., 2017.

% This takes our previous project, and transfers it into a two-group
% problem. We should be able to mostly recycle the code from Project 1, but
% due to the complicated loop used, this code is being written mostly from
% scratch.

clearvars, clc, format compact

%% Section 1: Defining Constants & User Inputs.
% This is where we define the necessary constants, as well as gather the
% user's desired properties such as height and radius. This will mostly be
% recycled from William's Project 1.

global del_r del_z

% ---------------
% Group 1 (Fast)
% ---------------

% In the core.
nusigf_1 = 0.008476;    % /cm
sig_r1 = 0.02619;        % /cm
D1 = 1.2627;            % cm
sigs_12 = 0.0494;
% In the Reflector.
nusigf_1w = 0;
sig_r1w = sigs_12;
D1w = 1.13;

% ---------------
% Group 2 (Thermal)
% ---------------

% In the Core.
nusigf_2 = 0.18514;
siga_2 = 0.1210;
D2 = 0.3543;

% In the Reflector.
nusigf_2w = 0;
D2w = 0.16;
sig_2w = 0.0197;

% ---------------
% User Input

prompt = {'Number of Radial Nodes (N):', 'Number of Axial Nodes (M):', 'Core Radius (cm):', 'Core Height (cm):','Number of Reflector Nodes:'};
title = 'Input';
defans = {'5','6','10','20','2'};
outp = inputdlg(prompt,title,1,defans);
if isempty(outp)
    return
end

% Extracting the User Input.

N = str2double(outp{1});
M = str2double(outp{2});
R = str2double(outp{3});
H = str2double(outp{4});
Ref = str2double(outp{5});

%  Derived Variables.
del_r = R./(N-1);     % cm
del_z = H./(M-1);     % cm

% Assigning the values to ri.
ri = zeros(1,N-1);

for i=1:(N-2)
    ri(i+1) = del_r.*(i);
end

%% Diffusion equation matrices Fast
% List of node distances

for i = 1:N-1
    r = ri(i);
    for j = 1:M-2   
      if i<N-Ref
        if i == 1 && j == 1 %Bottom Left Node
            Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) =  ((D1*del_r*del_z^2)/2 + (D1*del_r^3)/4 + (sig_r1*del_r^3*del_z^2)/8);                    % (i,j)
            Af((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-(D1*del_r*del_z^2)/2);                                                          % (i+1,j)
            Af((j-1)*(N-1)+i,j*(N-1)+i) = (-(D1*del_r^3)/8);                                                                       % (i,j+1)
        elseif i == 1             
            if j == M-2 %Top Left Node
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((D1*del_r*del_z^2)/2 + (D1*del_r^3)/4 + (sig_r1*del_r^3*del_z^2)/8);                % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-(D1*del_r*del_z^2)/2);                                                       % (i+1,j)
                Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-(D1*del_r^3)/8);                                                               % (i,j-1)
            else % Left Most Column
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((D1*del_r*del_z^2)/2 + (D1*del_r^3)/4 + (sig_r1*del_r^3*del_z^2)/8);                % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-(D1*del_r*del_z^2)/2);                                                       % (i+1,j)
                Af((j-1)*(N-1)+i,j*(N-1)+i) = (-(D1*del_r^3)/8);                                                                   % (i,j+1)
                Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-(D1*del_r^3)/8);                                                               % (i,j-1)
            end            
        elseif i ~= 1 && j == 1            
            if i == N-1 %Bottom Right Node
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1*r*del_z^2) + (2*D1*r*del_r^2) + (sig_r1*r*del_r^2*del_z^2));  % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,j*(N-1)+i) = (-D1*r*(del_r^2));                                                              % (i,j+1)
            else %Bottom Row
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1*r*del_z^2) + (2*D1*r*del_r^2) + (sig_r1*r*del_r^2*del_z^2));    % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,j*(N-1)+i) = (-D1*r*(del_r^2));                                                               % (i,j+1)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D1*(r+(del_r/2))*(del_z^2));                                              % (i+1,j)
            end
        elseif i ~= 1 && j == M-2           
            if i == N-1 %Top Right Node
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1*r*del_z^2) + (2*D1*r*del_r^2) + (sig_r1*r*del_r^2*del_z^2));  % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D1*r*(del_r^2));                                                          % (i,j-1)
            else %Top Row
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1*r*del_z^2) + (2*D1*r*del_r^2) + (sig_r1*r*del_r^2*del_z^2));  % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D1*(r+(del_r/2))*(del_z^2));                                              % (i+1,j)
                Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D1*r*(del_r^2));                                                          % (i,j-1)
            end
        elseif j ~= 1 && i == N-1        
            if j == M-2 %Top Right Node. This is a repeat
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1*r*del_z^2) + (2*D1*r*del_r^2) + (sig_r1*r*del_r^2*del_z^2));  % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) =(-D1*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D1*r*(del_r^2));                                                          % (i,j-1)
            else %Far Right Side nodes
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1*r*del_z^2) + (2*D1*r*del_r^2) + (sig_r1*r*del_r^2*del_z^2));  % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D1*r*(del_r^2));                                                          % (i,j-1)
                Af((j-1)*(N-1)+i,j*(N-1)+i) = (-D1*r*(del_r^2));                                                               % (i,j+1)
            end
        elseif j <= M-2 % all of the interior Nodes
            Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1*r*del_z^2) + (2*D1*r*del_r^2) + (sig_r1*r*del_r^2*del_z^2));      % (i,j)
            Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1*(r-(del_r/2))*(del_z^2));                                                  % (i-1,j)
            Af((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D1*(r+(del_r/2))*(del_z^2));                                                  % (i+1,j)
            Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D1*r*(del_r^2));                                                              % (i,j-1)
            Af((j-1)*(N-1)+i,j*(N-1)+i) = (-D1*r*(del_r^2));                                                                   % (i,j+1
        else
        end
      else   %% Reflector       
        if i ~= 1 && j == 1            
            if i == N-1 %Bottom Right Node
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1w*r*del_z^2) + (2*D1w*r*del_r^2) + (sig_r1w*r*del_r^2*del_z^2));  % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,j*(N-1)+i) = (-D1w*r*(del_r^2));                                                              % (i,j+1)
            else %Bottom Row
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1w*r*del_z^2) + (2*D1w*r*del_r^2) + (sig_r1w*r*del_r^2*del_z^2));    % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,j*(N-1)+i) = (-D1w*r*(del_r^2));                                                               % (i,j+1)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D1w*(r+(del_r/2))*(del_z^2));                                              % (i+1,j)
            end
        elseif i ~= 1 && j == M-2           
            if i == N-1 %Top Right Node
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1w*r*del_z^2) + (2*D1w*r*del_r^2) + (sig_r1w*r*del_r^2*del_z^2));  % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D1w*r*(del_r^2));                                                          % (i,j-1)
            else %Top Row
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1w*r*del_z^2) + (2*D1w*r*del_r^2) + (sig_r1w*r*del_r^2*del_z^2));  % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D1w*(r+(del_r/2))*(del_z^2));                                              % (i+1,j)
                Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D1w*r*(del_r^2));                                                          % (i,j-1)
            end
        elseif j ~= 1 && i == N-1        
            if j == M-2 %Top Right Node. This is a repeat
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1w*r*del_z^2) + (2*D1w*r*del_r^2) + (sig_r1w*r*del_r^2*del_z^2));  % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) =(-D1w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D1w*r*(del_r^2));                                                          % (i,j-1)
            else %Far Right Side nodes
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1w*r*del_z^2) + (2*D1w*r*del_r^2) + (sig_r1w*r*del_r^2*del_z^2));  % (i,j)
                Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D1w*r*(del_r^2));                                                          % (i,j-1)
                Af((j-1)*(N-1)+i,j*(N-1)+i) = (-D1w*r*(del_r^2));                                                               % (i,j+1)
            end
        elseif j <= M-2 % all of the interior Nodes
            Af((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D1w*r*del_z^2) + (2*D1w*r*del_r^2) + (sig_r1w*r*del_r^2*del_z^2));      % (i,j)
            Af((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D1w*(r-(del_r/2))*(del_z^2));                                                  % (i-1,j)
            Af((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D1w*(r+(del_r/2))*(del_z^2));                                                  % (i+1,j)
            Af((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D1w*r*(del_r^2));                                                              % (i,j-1)
            Af((j-1)*(N-1)+i,j*(N-1)+i) = (-D1w*r*(del_r^2)); 
          end
      end
    end    
end


%% Diffusion equation matrices Thermal
% List of node distances

for i = 1:N-1
    r = ri(i);
    for j = 1:M-2
      if i<N-Ref
        if i == 1 && j == 1 %Bottom Left Node
            At((j-1)*(N-1)+i,(j-1)*(N-1)+i) =  ((D2*del_r*del_z^2)/2 + (D2*del_r^3)/4 + (siga_2*del_r^3*del_z^2)/8);                    % (i,j)
            At((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-(D2*del_r*del_z^2)/2);                                                          % (i+1,j)
            At((j-1)*(N-1)+i,j*(N-1)+i) = (-(D2*del_r^3)/8);                                                                       % (i,j+1)
        elseif i == 1             
            if j == M-2 %Top Left Node
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((D2*del_r*del_z^2)/2 + (D2*del_r^3)/4 + (siga_2*del_r^3*del_z^2)/8);                % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-(D2*del_r*del_z^2)/2);                                                       % (i+1,j)
                At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-(D2*del_r^3)/8);                                                               % (i,j-1)
            else % Left Most Column
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((D2*del_r*del_z^2)/2 + (D2*del_r^3)/4 + (siga_2*del_r^3*del_z^2)/8);                % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-(D2*del_r*del_z^2)/2);                                                       % (i+1,j)
                At((j-1)*(N-1)+i,j*(N-1)+i) = (-(D2*del_r^3)/8);                                                                   % (i,j+1)
                At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-(D2*del_r^3)/8);                                                               % (i,j-1)
            end            
        elseif i ~= 1 && j == 1            
            if i == N-1 %Bottom Right Node
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2*r*del_z^2) + (2*D2*r*del_r^2) + (siga_2*r*del_r^2*del_z^2));  % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,j*(N-1)+i) = (-D2*r*(del_r^2));                                                              % (i,j+1)
            else %Bottom Row
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2*r*del_z^2) + (2*D2*r*del_r^2) + (siga_2*r*del_r^2*del_z^2));    % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,j*(N-1)+i) = (-D2*r*(del_r^2));                                                               % (i,j+1)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D2*(r+(del_r/2))*(del_z^2));                                              % (i+1,j)
            end
        elseif i ~= 1 && j == M-2            
            if i == N-1 %Top Right Node
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2*r*del_z^2) + (2*D2*r*del_r^2) + (siga_2*r*del_r^2*del_z^2));  % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D2*r*(del_r^2));                                                          % (i,j-1)
            else %Top Row
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2*r*del_z^2) + (2*D2*r*del_r^2) + (siga_2*r*del_r^2*del_z^2));  % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D2*(r+(del_r/2))*(del_z^2));                                              % (i+1,j)
                At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D2*r*(del_r^2));                                                          % (i,j-1)
            end
        elseif j ~= 1 && i == N-1        
            if j == M-2 %Top Right Node. This is a repeat
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2*r*del_z^2) + (2*D2*r*del_r^2) + (siga_2*r*del_r^2*del_z^2));  % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) =(-D2*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D2*r*(del_r^2));                                                          % (i,j-1)
            else %Far Right Side nodes
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2*r*del_z^2) + (2*D2*r*del_r^2) + (siga_2*r*del_r^2*del_z^2));  % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D2*r*(del_r^2));                                                          % (i,j-1)
                At((j-1)*(N-1)+i,j*(N-1)+i) = (-D2*r*(del_r^2));                                                               % (i,j+1)
            end
        elseif j <= M-2 % all of the interior Nodes
            At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2*r*del_z^2) + (2*D2*r*del_r^2) + (siga_2*r*del_r^2*del_z^2));      % (i,j)
            At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2*(r-(del_r/2))*(del_z^2));                                                  % (i-1,j)
            At((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D2*(r+(del_r/2))*(del_z^2));                                                  % (i+1,j)
            At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D2*r*(del_r^2));                                                              % (i,j-1)
            At((j-1)*(N-1)+i,j*(N-1)+i) = (-D2*r*(del_r^2));                                                                   % (i,j+1)
        else
        end
    else   %% Reflector       
        if i ~= 1 && j == 1            
            if i == N-1 %Bottom Right Node
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2w*r*del_z^2) + (2*D2w*r*del_r^2) + (sig_2w*r*del_r^2*del_z^2));  % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,j*(N-1)+i) = (-D2w*r*(del_r^2));                                                              % (i,j+1)
            else %Bottom Row
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2w*r*del_z^2) + (2*D2w*r*del_r^2) + (sig_2w*r*del_r^2*del_z^2));    % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,j*(N-1)+i) = (-D2w*r*(del_r^2));                                                               % (i,j+1)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D2w*(r+(del_r/2))*(del_z^2));                                              % (i+1,j)
            end
        elseif i ~= 1 && j == M-2           
            if i == N-1 %Top Right Node
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2w*r*del_z^2) + (2*D2w*r*del_r^2) + (sig_2w*r*del_r^2*del_z^2));  % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D2w*r*(del_r^2));                                                          % (i,j-1)
            else %Top Row
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2w*r*del_z^2) + (2*D2w*r*del_r^2) + (sig_2w*r*del_r^2*del_z^2));  % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D2w*(r+(del_r/2))*(del_z^2));                                              % (i+1,j)
                At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D2w*r*(del_r^2));                                                          % (i,j-1)
            end
        elseif j ~= 1 && i == N-1        
            if j == M-2 %Top Right Node. This is a repeat
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2w*r*del_z^2) + (2*D2w*r*del_r^2) + (sig_2w*r*del_r^2*del_z^2));  % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) =(-D2w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D2w*r*(del_r^2));                                                          % (i,j-1)
            else %Far Right Side nodes
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2w*r*del_z^2) + (2*D2w*r*del_r^2) + (sig_2w*r*del_r^2*del_z^2));  % (i,j)
                At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2w*(r-(del_r/2))*(del_z^2));                                              % (i-1,j)
                At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D2w*r*(del_r^2));                                                          % (i,j-1)
                At((j-1)*(N-1)+i,j*(N-1)+i) = (-D2w*r*(del_r^2));                                                               % (i,j+1)
            end
        elseif j <= M-2 % all of the interior Nodes
            At((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((2*D2w*r*del_z^2) + (2*D2w*r*del_r^2) + (sig_2w*r*del_r^2*del_z^2));      % (i,j)
            At((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = (-D2w*(r-(del_r/2))*(del_z^2));                                                  % (i-1,j)
            At((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = (-D2w*(r+(del_r/2))*(del_z^2));                                                  % (i+1,j)
            At((j-1)*(N-1)+i,(j-2)*(N-1)+i) = (-D2w*r*(del_r^2));                                                              % (i,j-1)
            At((j-1)*(N-1)+i,j*(N-1)+i) = (-D2w*r*(del_r^2)); 
        end
    end   
    end
end

%% B Matrices

Bf = zeros(N,M);

for i = 1:N-1;
     r = ri(i);
    for j = 1:M-2;
        if i == 1 %Centerline B
            Bf((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((nusigf_1.*(del_r.^3).*(del_z.^2))./8);
        elseif i<N-Ref; %Every Other B
            Bf((j-1)*(N-1)+i,(j-1)*(N-1)+i) = (nusigf_1.*r.*(del_r.^2).*(del_z.^2));
        else
            Bf((j-1)*(N-1)+i,(j-1)*(N-1)+i) = 0;
        end
    end
end

Bt1 = zeros(N,M);

for i = 1:N-1;
     r = ri(i);
    for j = 1:M-2;
        if i == 1 %Centerline B
            Bt1((j-1)*(N-1)+i,(j-1)*(N-1)+i) = ((nusigf_2.*(del_r.^3).*(del_z.^2))./8);
        elseif i<N-Ref; %Every Other B
            Bt1((j-1)*(N-1)+i,(j-1)*(N-1)+i) = (nusigf_2.*r.*(del_r.^2).*(del_z.^2));
        else
            Bt1((j-1)*(N-1)+i,(j-1)*(N-1)+i) = 0;
        end
    end
end

Bt2 = zeros(N,M);

for i = 1:N-1;
     r = ri(i);
    for j = 1:M-2;
        if i == 1 %Centerline B
            Bt2((j-1)*(N-1)+i,(j-1)*(N-1)+i) = (sig_r1*del_r^3*del_z^2)/8;
        elseif i<N-Ref; %Every Other B
            Bt2((j-1)*(N-1)+i,(j-1)*(N-1)+i) = (sig_r1.*r.*(del_r.^2).*(del_z.^2));
        else
            Bt2((j-1)*(N-1)+i,(j-1)*(N-1)+i) = (sigs_12.*r.*(del_r.^2).*(del_z.^2));
        end
    end
end

%% Solving for flux iteratively
% Initial guess for flux, k-eff, and source
% Define initial flux guess
fluxf = ones((N-1)*(M-2),1);
fluxt = ones((N-1)*(M-2),1);

% Define initial k and source guess
k = 1;
Sf = (1/k)*(Bt1*fluxf);
St = (Bf*fluxt);

% Iterate to find flux, k, and S
fluxdiff = 1.0;

while (fluxdiff > 0.0001)
    oldfluxf = fluxf;
    oldfluxt = fluxt;
    oldk = k;
    oldSf = Sf;
    oldSt = St;
    
    fluxf = Af\(oldSf*(1/oldk));
    fluxt = At\oldSt;
    Sf = (Bf*fluxf+Bt1*fluxt);
    St = (Bt2*fluxf);
    k = oldk*(sum(Sf)/sum(oldSf));
    
    fluxdiff = 0.0;
    
    for i = 1:(N-1)*(M-2)
        fluxdiff = fluxdiff + (fluxf(i)-oldfluxf(i)).^2;
    end
    
    fluxdiff = sqrt(fluxdiff);
    
end

% % Normalize flux to 3000MW
MW3000=(1.809e10)/(3.24*10e-11*(0.07537+0.003320)*(pi*R^2*H));
fluxcomb=fluxf+fluxt;
totc=sum(fluxcomb);
tot = sum(fluxf);
fluxf = (fluxf/totc)*MW3000;
tott = sum(fluxt);
fluxt = (fluxt/totc)*MW3000;

% Initalize Matrix for Graphing
Graph1 = zeros(N+1,M-1);
Graph2 = zeros(N+1,M-1);

for i = 1:length(fluxt)
    if mod(i,(N-1)) == 0
        Graph1(ceil(i/(N-1))+1,(N-1)) = fluxt(i);
    else
        Graph1(ceil(i/(N-1))+1,mod(i,N-1)) = fluxt(i);
    end
end

Graph2 = fliplr(Graph1);
Graph = horzcat(Graph2,Graph1);

% Plotting
bar3(Graph)
xlabel('Node Number (N)')
ylabel('Node Number (M)')
k