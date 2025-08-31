function C = GetCprime( E1, E2, E3, nu12, nu13, nu23, G12, G13, G23 )
    % Compute the reduced stiffness matrix C' from engineering constants
    % Following Eq.(6)-(7) from Chou, Carleone & Hsu (1972).

    % Denominator V (Eq. 7)
    denom = 1.0 ...
          - nu12 * ( nu12 * (E2/E1) + 2.0 * nu23 * nu13 * (E3/E1) ) ...
          - nu13 * nu13 * (E3/E1) ...
          - nu23 * nu23 * (E3/E2);

    % C' components (Eq. 6)
    C11 = ( 1.0 - nu23 * nu23 * (E3/E2) ) * E1 / denom;
    C12 = ( nu12 + nu13 * nu23 * (E3/E2) ) * E2 / denom;
    C13 = ( nu13 + nu12 * nu23 ) * E3 / denom;
    C22 = ( 1.0 - nu13 * nu13 * (E3/E1) ) * E2 / denom;
    C23 = ( nu23 - nu12 * nu13 * (E2/E1) ) * E3 / denom;
    C33 = ( 1.0 - nu12 * nu12 * (E2/E1) ) * E3 / denom;

    C44 = G23;
    C55 = G13;
    C66 = G12;

    % Assemble symmetric stiffness matrix in Voigt notation
    C = [ C11, C12, C13,   0,   0,   0;
          C12, C22, C23,   0,   0,   0;
          C13, C23, C33,   0,   0,   0;
            0,   0,   0, C44,   0,   0;
            0,   0,   0,   0, C55,   0;
            0,   0,   0,   0,   0, C66 ];
end
