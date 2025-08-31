function L = build_L_stress( theta_deg )
    th = theta_deg * pi / 180.0;

    m = cos(th);
    n = sin(th);

    L = [ m * m,     n * n,   0, 0, 0,   2.0 * m * n;
          n * n,     m * m,   0, 0, 0,  -2.0 * m * n;
              0,         0,   1, 0, 0,            0;
              0,         0,   0, m, -n,           0;
              0,         0,   0, n,  m,           0;
         -m * n,     m * n,   0, 0, 0,  m * m - n * n ];
end
