import numpy as np

# -------------------- Rotation matrices (stress & strain transforms) --------------------

def build_L_stress(theta_deg: float) -> np.ndarray:
    """Build stress transform L(theta) per paper Eq. (13)."""
    th = np.deg2rad(theta_deg)
    m = np.cos(th)
    n = np.sin(th)

    L = np.array([
        [ m*m,   n*n,   0., 0., 0.,   2.0*m*n],
        [ n*n,   m*m,   0., 0., 0.,  -2.0*m*n],
        [ 0.,     0.,   1., 0., 0.,        0.],
        [ 0.,     0.,   0., m, -n,        0.],
        [ 0.,     0.,   0., n,  m,        0.],
        [-m*n,  m*n,   0., 0., 0., m*m - n*n]
    ], dtype=float)
    return L


def build_T_strain(theta_deg: float) -> np.ndarray:
    """Build strain transform T(theta) per paper Eq. (16)."""
    th = np.deg2rad(theta_deg)
    m = np.cos(th)
    n = np.sin(th)

    T = np.array([
        [    m*m,     n*n, 0., 0., 0.,     m*n],
        [    n*n,     m*m, 0., 0., 0.,    -m*n],
        [    0.,       0., 1., 0., 0.,       0.],
        [    0.,       0., 0., m,  n,       0.],
        [    0.,       0., 0., -n,  m,       0.],
        [-2.0*m*n, 2.0*m*n, 0., 0., 0., m*m - n*n]
    ], dtype=float)
    return T


# -------------------- C' from engineering constants (Eq. 6–7) --------------------

def GetCprime(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23) -> np.ndarray:
    """Compute reduced stiffness matrix C' from engineering constants."""
    # Denominator V (Eq. 7)
    denom = (1.0
             - nu12 * (nu12 * (E2/E1) + 2.0 * nu23 * nu13 * (E3/E1))
             - nu13 * nu13 * (E3/E1)
             - nu23 * nu23 * (E3/E2))

    # C' components (Eq. 6)
    C11 = (1.0 - nu23 * nu23 * (E3/E2)) * E1 / denom
    C12 = (nu12 + nu13 * nu23 * (E3/E2)) * E2 / denom
    C13 = (nu13 + nu12 * nu23) * E3 / denom
    C22 = (1.0 - nu13 * nu13 * (E3/E1)) * E2 / denom
    C23 = (nu23 - nu12 * nu13 * (E2/E1)) * E3 / denom
    C33 = (1.0 - nu12 * nu12 * (E2/E1)) * E3 / denom

    C44 = G23
    C55 = G13
    C66 = G12

    # Assemble symmetric stiffness matrix in Voigt notation
    C = np.array([
        [C11, C12, C13,   0.,   0.,   0.],
        [C12, C22, C23,   0.,   0.,   0.],
        [C13, C23, C33,   0.,   0.,   0.],
        [0.,   0.,   0., C44,   0.,   0.],
        [0.,   0.,   0.,   0., C55,   0.],
        [0.,   0.,   0.,   0.,   0., C66]
    ], dtype=float)
    return C


# -------------------- Effective stiffness per Chou/Carleone/Hsu (Eqs.25–28) --------------------

def LaminateTheory3D(layers):
    """
    Compute effective stiffness of a laminate.
    layers: list of dicts with keys:
        - 'thickness': ply thickness
        - 'theta': ply orientation in degrees
        - 'mat': dict with keys E1,E2,E3,nu12,nu13,nu23,G12,G13,G23
    Returns:
        Cstar: effective stiffness matrix (6x6)
        Hstar: compliance matrix (6x6)
        props: dict of effective engineering constants
    """
    nlayer = len(layers)

    # total thickness
    t_sum = sum(l['thickness'] for l in layers)

    # volume fractions
    V = np.array([l['thickness']/t_sum for l in layers])

    # rotated Cbar for each layer
    Cbars = []
    for lay in layers:
        mat = lay['mat']
        Cprime = GetCprime(mat['E1'], mat['E2'], mat['E3'],
                           mat['nu12'], mat['nu13'], mat['nu23'],
                           mat['G12'], mat['G13'], mat['G23'])

        L = build_L_stress(lay['theta'])
        T = build_T_strain(lay['theta'])
        Cbar = L @ Cprime @ np.linalg.inv(T)
        Cbars.append(Cbar)
    Cbars = np.stack(Cbars, axis=2)  # shape (6,6,nlayer)

    Iidx = [0,1,2,5]  # 0-based indices for Voigt order
    Sidx = [3,4]

    Cstar = np.zeros((6,6))

    # precompute topj and bot
    topj = np.zeros(6)
    bot = 0.0
    for l in range(nlayer):
        Cl = Cbars[:,:,l]
        bot += V[l] / Cl[2,2]
        topj += V[l] * Cl[2,:] / Cl[2,2]

    # in-plane block assembly
    for k in range(nlayer):
        Ck = Cbars[:,:,k]
        for j in Iidx:
            for i in Iidx:
                Cstar[i,j] += V[k] * (
                    Ck[i,j]
                    - Ck[i,2]*Ck[2,j]/Ck[2,2]
                    + Ck[i,2]*topj[j]/(Ck[2,2]*bot)
                )

    # shear block determinants
    detv = np.zeros(nlayer)
    for k in range(nlayer):
        Ck = Cbars[:,:,k]
        detv[k] = Ck[3,3]*Ck[4,4] - Ck[3,4]*Ck[4,3]

    # shear denom
    denom = 0.0
    for k in range(nlayer):
        Ck = Cbars[:,:,k]
        for l in range(nlayer):
            Cl = Cbars[:,:,l]
            denom += (V[k]*V[l]/(detv[k]*detv[l])) * (Ck[3,3]*Cl[4,4] - Ck[3,4]*Cl[4,3])

    # shear block assembly
    for j in Sidx:
        for i in Sidx:
            for k in range(nlayer):
                Ck = Cbars[:,:,k]
                Cstar[i,j] += (V[k]*Ck[i,j]/detv[k]) / denom

    # symmetrize
    Cstar = 0.5*(Cstar + Cstar.T)

    # compliance and engineering constants
    Hstar = np.linalg.inv(Cstar)

    props = {}
    props['Ex']  = 1.0/Hstar[0,0]
    props['Ey']  = 1.0/Hstar[1,1]
    props['Ez']  = 1.0/Hstar[2,2]
    props['vyz'] = -Hstar[1,2]/Hstar[1,1]
    props['vxz'] = -Hstar[0,2]/Hstar[0,0]
    props['vxy'] = -Hstar[0,1]/Hstar[0,0]
    props['Gyz'] = 1.0/Hstar[3,3]
    props['Gxz'] = 1.0/Hstar[4,4]
    props['Gxy'] = 1.0/Hstar[5,5]

    # print like C++ version
    print(f"Ex {props['Ex']/1e6:.6g}")
    print(f"Ey {props['Ey']/1e6:.6g}")
    print(f"Ez {props['Ez']/1e6:.6g}")
    print(f"nuyz {props['vyz']:.6g}")
    print(f"nuxz {props['vxz']:.6g}")
    print(f"nuzy {props['vxy']:.6g}")
    print(f"Gyz {props['Gyz']/1e6:.6g}")
    print(f"Gxz {props['Gxz']/1e6:.6g}")
    print(f"Gxy {props['Gxy']/1e6:.6g}\n")

    return Cstar, Hstar, props


# -------------------- Example usage --------------------
if __name__ == "__main__":
    # ---- Material definition ----
    # Carbon Epoxy AS4 3501-6
    mat = {
        "E1": 16.5e6,   # psi
        "E2": 1.4e6,
        "E3": 1.4e6,
        "nu23": 0.54,
        "nu13": 0.33,
        "nu12": 0.33,
        "G23": 0.45e6,
        "G13": 0.87e6,
        "G12": 0.87e6,
    }

    # ---- Lamina definition ----
    lay = {
        "thickness": 0.1,
        "theta": 0,
        "mat": mat,
    }

    # ---- Stacking sequences ----
    thetas1 = [0, 0, 90, 90, 0, 0]
    thetas2 = [0, 90, 45, -45, -45, 45, 90, 0]
    thetas3 = [30, -30, -30, 30]

    layers1 = []
    for th in thetas1:
        layer = lay.copy()
        layer["theta"] = th
        layers1.append(layer)

    layers2 = []
    for th in thetas2:
        layer = lay.copy()
        layer["theta"] = th
        layers2.append(layer)

    layers3 = []
    for th in thetas3:
        layer = lay.copy()
        layer["theta"] = th
        layers3.append(layer)

    # ---- Call laminate theory ----
    Cstar1, Hstar1, props1 = LaminateTheory3D(layers1)
    Cstar2, Hstar2, props2 = LaminateTheory3D(layers2)
    Cstar3, Hstar3, props3 = LaminateTheory3D(layers3)


