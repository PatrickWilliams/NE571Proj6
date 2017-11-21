clear all;
close all;
clc;

R = input('Input Radius size in Cm: ');
H = input('Input Height size in Cm: ');
N = input('Input Number of Axial Nodes: ');
M = input('Input Number of Vertical Nodes: ');

fuel_length = input('Input radius of fuel: ');
reflector_length = input('Input reflector length: ');

dr = R/(N-1);
dz = H/(M-1);

R_vec = 0:dr:R;

fuel_nodes = ceil(fuel_length/dr);
materials = zeros(M,N,4);

for i = 1:fuel_nodes
    materials(:,i,1) = 
    materials(:,i,2) =
    materials(:,i,3) =
    materials(:,i,4) =
end


D = 1/(3.62*10^-2);
E_a = .1532;
v = 1;
E_f = .1570;
%Populate LHS constant Matrix


for i = 1:N-1
    for j = 1:M-2
        
        if (i-1) == 0 && (j-1)==0
            
            
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i) = E_a*(dr^3)/8*dz^2 + 2*D*(dr^3)/8 + D*dr/2*dz^2;% (i,j)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = -D*dr/2*dz^2;   %(i+1,j)
            A((j-1)*(N-1)+i,(j)*(N-1)+i) = -D*(dr^3)/8;  %(i,j+1)
            
            B((j-1)*(N-1)+i,(j-1)*(N-1)+i) = v*E_f/8*dr^3*dz^2;
            
        elseif (i-1) == 0
            if (j == M-2)
                A((j-1)*(N-1)+i,(j-2)*(N-1) +i) = -(D/8)*dr^3; % (i,j-1)
                
                A((j-1)*(N-1)+i,(j-1)*(N-1)+i) = E_a*(dr^3)/8*dz^2 + 2*D*(dr^3)/8 + D*dr/2*dz^2;% (i,j)
                A((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = -D*dr/2*dz^2;   %(i+1,j)
                
                
                B((j-1)*(N-1)+i,(j-1)*(N-1)+i) = v*E_f/8*dr^3*dz^2;
                
            else
                A((j-1)*(N-1)+i,(j-2)*(N-1) +i) = -(D/8)*dr^3; % (i,j-1)
                A((j-1)*(N-1)+i,(j-1)*(N-1)+i) = E_a*(dr^3)/8*dz^2 + 2*D*(dr^3)/8 + D*dr/2*dz^2;% (i,j)
                A((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = -D*dr/2*dz^2;   %(i+1,j)
                A((j-1)*(N-1)+i,(j)*(N-1)+i) = -D*(dr^3)/8;  %(i,j+1)
                
                B((j-1)*(N-1)+i,(j-1)*(N-1)+i) = v*E_f/8*dr^3*dz^2;
            end
        elseif (j-1) == 0
            if (i == N-1)
                
                A((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = -D*(R_vec(i)-dr/2)*dz^2;% (i-1,j)
                A((j-1)*(N-1)+i,(j-1)*(N-1)+i) = E_a*R_vec(i)*dr^2*dz^2 + 2*D*R_vec(i)*dr^2 + 2*D*R_vec(i)*dz^2;% (i,j)
                A((j-1)*(N-1)+i,(j)*(N-1)+i) = -D*R_vec(i)*dr^2;   %(i,j+1)
                
                B((j-1)*(N-1)+i,(j-1)*(N-1)+i) = v*E_f*R_vec(i)*dr^2*dz^2;
            else
                A((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = -D*(R_vec(i)-dr/2)*dz^2;% (i-1,j)
                A((j-1)*(N-1)+i,(j-1)*(N-1)+i) = E_a*R_vec(i)*dr^2*dz^2 + 2*D*R_vec(i)*dr^2 + 2*D*R_vec(i)*dz^2;% (i,j)
                A((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = -D*(R_vec(i) + dr/2)*dz^2;   %(i+1,j)
                A((j-1)*(N-1)+i,(j)*(N-1)+i) = -D*R_vec(i)*dr^2;   %(i,j+1)
                
                B((j-1)*(N-1)+i,(j-1)*(N-1)+i) = v*E_f*R_vec(i)*dr^2*dz^2;
            end
            
        elseif (i == N-1) && (j==M-2)
            A((j-1)*(N-1)+i,(j-2)*(N-1) +i) = -D*R_vec(i)*dr^2;% (i,j-1)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = -D*(R_vec(i)-dr/2)*dz^2;% (i-1,j)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i) = E_a*R_vec(i)*dr^2*dz^2 + 2*D*R_vec(i)*dr^2 + 2*D*R_vec(i)*dz^2;% (i,j)
            
            
            B((j-1)*(N-1)+i,(j-1)*(N-1)+i) = v*E_f*R_vec(i)*dr^2*dz^2;
            
        elseif (i == N-1)
            A((j-1)*(N-1)+i,(j-2)*(N-1) +i) = -D*R_vec(i)*dr^2;% (i,j-1)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = -D*(R_vec(i)-dr/2)*dz^2;% (i-1,j)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i) = E_a*R_vec(i)*dr^2*dz^2 + 2*D*R_vec(i)*dr^2 + 2*D*R_vec(i)*dz^2;% (i,j)
            
            A((j-1)*(N-1)+i,(j)*(N-1)+i) = -D*R_vec(i)*dr^2;   %(i,j+1)
            
            B((j-1)*(N-1)+i,(j-1)*(N-1)+i) = v*E_f*R_vec(i)*dr^2*dz^2;
            
        elseif (j == M-2)
            A((j-1)*(N-1)+i,(j-2)*(N-1) +i) = -D*R_vec(i)*dr^2;% (i,j-1)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = -D*(R_vec(i)-dr/2)*dz^2;% (i-1,j)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i) = E_a*R_vec(i)*dr^2*dz^2 + 2*D*R_vec(i)*dr^2 + 2*D*R_vec(i)*dz^2;% (i,j)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = -D*(R_vec(i) + dr/2)*dz^2;   %(i+1,j)
            
            
            B((j-1)*(N-1)+i,(j-1)*(N-1)+i) = v*E_f*R_vec(i)*dr^2*dz^2;
            
        else %Interior Nodes
            
            A((j-1)*(N-1)+i,(j-2)*(N-1) +i) = -D*R_vec(i)*dr^2;% (i,j-1)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i-1) = -D*(R_vec(i)-dr/2)*dz^2;% (i-1,j)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i) = E_a*R_vec(i)*dr^2*dz^2 + 2*D*R_vec(i)*dr^2 + 2*D*R_vec(i)*dz^2;% (i,j)
            A((j-1)*(N-1)+i,(j-1)*(N-1)+i+1) = -D*(R_vec(i) + dr/2)*dz^2;   %(i+1,j)
            A((j-1)*(N-1)+i,(j)*(N-1)+i) = -D*R_vec(i)*dr^2;   %(i,j+1)
            
            B((j-1)*(N-1)+i,(j-1)*(N-1)+i) = v*E_f*R_vec(i)*dr^2*dz^2;
            
        end
    end
end

%%%%%%
% Decleare flux all ones initially
Flux=ones((N-1)*(M-2),1);
K=1.0;

S=(1/K)*B*Flux; % Calculates right-hand side

flxdiff=1.0;

inverted_a  = inv(A);

% Iterations
while(flxdiff>0.0001)
    OldFlux=Flux;
    OldK=K;
    OldS=S;
    Flux=inverted_a*OldS;
    
    S=1/OldK*B*Flux;
    K=OldK*sum(S)/sum(OldS);
    
    
    % Calculate flux difference between successive iterations
    flxdiff=0.0;
    for i=1:(N-1)*(M-2)
        flxdiff = flxdiff + (Flux(i)-OldFlux(i))^2;
    end
    flxdiff = sqrt(flxdiff);
end
% Normalize Flux (Integral of flux set to unity)
%Total = sum(Flux);
%Flux = Flux/Total;

% Seperate into groups and Add last node (zero flux condition)

R = 0:dr:R;
H = 0:dz:H;

Flux_graph = zeros(N,M);


for i = 1:length(Flux)
    if mod(i,(N-1)) == 0
        Flux_graph(ceil(i/(N-1))+1,(N-1)) = Flux(i);
    else
        Flux_graph(ceil(i/(N-1))+1,mod(i,N-1)) = Flux(i);
    end
end

x_vec = 0:N;
y_vec = 0:M;

% Plot the final plot

bar3(Flux_graph)

xlabel('Radial Nodes (N)')
ylabel('Axial Nodes (M)')
legend('Flux')


