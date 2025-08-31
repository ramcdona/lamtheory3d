function [Cstar, Hstar, props] = LaminateTheory3D(layers)
%LAMINATETHEORY3D  Effective stiffness via Chou/Carleone/Hsu (Eqs. 25–28).
%   layers is an array of structs with fields:
%     - Thickness_FEM  : scalar
%     - Theta_deg      : ply angle in degrees
%     - Mat            : struct with fields E1,E2,E3,nu12,nu13,nu23,G12,G13,G23
%
%   Returns:
%     Cstar : 6x6 effective stiffness (Voigt order: 11,22,33,23,13,12)
%     Hstar : 6x6 effective compliance (inv(Cstar))
%     props : struct with effective engineering constants Ex,Ey,Ez, vyz,vxz,vxy, Gyz,Gxz,Gxy


nlayer = numel(layers);

% Total thickness
t_sum = 0.0;
for ilay = 1:nlayer
    if ~isempty(layers(ilay))
        t_sum = t_sum + layers(ilay).Thickness_FEM;
    end
end

% Volume fractions V(l)
V = zeros(1, nlayer);
for ilay = 1:nlayer
    if ~isempty(layers(ilay))
        V(ilay) = layers(ilay).Thickness_FEM / t_sum;
    end
end

% Build rotated layer stiffnesses Cbar(l)
Cbars = zeros(6,6,nlayer);
for ilay = 1:nlayer
    lay = layers(ilay);

    mat = lay.Mat;

    Cprime = GetCprime( ...
                        mat.E1, mat.E2, mat.E3, ...
                        mat.nu12, mat.nu13, mat.nu23, ...
                        mat.G12, mat.G13, mat.G23);

    % rotate each C' into the global frame
    L = build_L_stress( lay.Theta_deg );
    T = build_T_strain( lay.Theta_deg );

    Cbar = L * Cprime / T;

    Cbars( :, :, ilay ) = Cbar;
end

% Voigt index sets
Iidx = [1 2 3 6];
Sidx = [4 5];

Cstar = zeros(6,6);

% Precompute topj and bot
topj = zeros(1,6);
bot  = 0.0;

for l = 1:nlayer
    Cl = Cbars(:,:,l);

    bot = bot + V(l) / Cl(3,3);

    for j = 1:6
        topj(j) = topj(j) + V(l) * Cl(3,j) / Cl(3,3);
    end
end

% In-plane block assembly using Iidx
for k = 1:nlayer
    Ck = Cbars(:,:,k);

    for ja = 1:4
        j = Iidx(ja);

        for ia = 1:4
            i = Iidx(ia);

            Cstar(i,j) = Cstar(i,j) + V(k) * ( ...
                         Ck(i,j) ...
                         - Ck(i,3) * Ck(3,j) / Ck(3,3) ...
                         + Ck(i,3) * topj(j) / ( Ck(3,3) * bot ) );
        end
    end
end

% Determinants of the shear 2x2 sub-blocks
detv = zeros(1,nlayer);
for k = 1:nlayer
    Ck = Cbars(:,:,k);
    detv(k) = Ck(4,4) * Ck(5,5) - Ck(4,5) * Ck(5,4);
end

% Denominator for shear part
denom = 0.0;
for k = 1:nlayer
    Ck = Cbars(:,:,k);
    for l = 1:nlayer
        Cl = Cbars(:,:,l);
        denom = denom + ( V(k)*V(l) / ( detv(k) * detv(l) ) ) * ...
            ( Ck(4,4) * Cl(5,5) - Ck(4,5) * Cl(5,4) );
    end
end

% Shear block assembly using Sidx
for ja = 1:2
    j = Sidx(ja);

    for ia = 1:2
        i = Sidx(ia);

        for k = 1:nlayer
            Ck = Cbars(:,:,k);
            Cstar(i,j) = Cstar(i,j) + ( V(k) * Ck(i,j) / detv(k) ) / denom;
        end
    end
end

% Ensure symmetry
Cstar = 0.5 * (Cstar + Cstar');

% Compliance and engineering constants
Hstar = inv(Cstar);

props = struct();
props.Ex  = 1.0 / Hstar(1,1);
props.Ey  = 1.0 / Hstar(2,2);
props.Ez  = 1.0 / Hstar(3,3);
props.vyz = -Hstar(2,3) / Hstar(2,2);
props.vxz = -Hstar(1,3) / Hstar(1,1);
props.vxy = -Hstar(1,2) / Hstar(1,1);
props.Gyz = 1.0 / Hstar(4,4);
props.Gxz = 1.0 / Hstar(5,5);
props.Gxy = 1.0 / Hstar(6,6);

% Printouts to match the table
fprintf('Ex %g\n',  props.Ex/1e6);
fprintf('Ey %g\n',  props.Ey/1e6);
fprintf('Ez %g\n',  props.Ez/1e6);

fprintf('nuyz %g\n', props.vyz);
fprintf('nuxz %g\n', props.vxz);
fprintf('nuzy %g\n', props.vxy);

fprintf('Gyz %g\n', props.Gyz/1e6);
fprintf('Gxz %g\n', props.Gxz/1e6);
fprintf('Gxy %g\n', props.Gxy/1e6);

fprintf('\n\n');
end


