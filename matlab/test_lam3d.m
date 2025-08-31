clear all
format compact
clc

% Carbon Epoxy AS4 3501-6
% AS4 3501-6 elasticity data from MIL-HDBK-17-3F p. 627.  59.5% fiber volume fraction
mat.E1 = 16.5e+06; % psi
mat.E2 = 1.4e+06;
mat.E3 = 1.4e+06;
mat.nu23 = 0.54;
mat.nu13 = 0.33;
mat.nu12 = 0.33;
mat.G23 = 0.45e6;
mat.G13 = 0.87e6;
mat.G12 = 0.87e6;


lay.Thickness_FEM = 0.1;
lay.Theta_deg = 0;
lay.Mat = mat;

thetas1 = [0 0 90 90 0 0];
thetas2 = [0 90 45 -45 -45 45 90 0];
thetas3 = [30 -30 -30 30];

for i=1:length( thetas1 )
    lay.Theta_deg = thetas1( i );
    layers1(i) = lay;
end

for i=1:length( thetas2 )
    lay.Theta_deg = thetas2( i );
    layers2(i) = lay;
end

for i=1:length( thetas3 )
    lay.Theta_deg = thetas3( i );
    layers3(i) = lay;
end

[Cstar, Hstar, props] = LaminateTheory3D( layers1 );
[Cstar, Hstar, props] = LaminateTheory3D( layers2 );
[Cstar, Hstar, props] = LaminateTheory3D( layers3 );

