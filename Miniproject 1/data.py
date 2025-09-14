# UNIQUAC Parameters (Bohorquez et al., 2020)

A = dict(
    CA={"CA": 0, "EtOH": 0, "MEC": 0, "H2O": 0, "DEC": 0, "TEC": 0},
    EtOH={"CA": 0, "EtOH": 0, "MEC": 0, "H2O": 2.0046, "DEC": 0, "TEC": 0},
    MEC={"CA": 0, "EtOH": 0, "MEC": 0, "H2O": 0, "DEC": 0, "TEC": 0},
    H2O={"CA": 0, "EtOH": -2.4936, "MEC": 0, "H2O": 0, "DEC": 0, "TEC": 0},
    DEC={"CA": 0, "EtOH": 0, "MEC": 0, "H2O": 0, "DEC": 0, "TEC": 0},
    TEC={"CA": 0, "EtOH": 0, "MEC": 0, "H2O": 0, "DEC": 0, "TEC": 0}
)

B = dict(
    CA={"CA": 0, "EtOH": -139.839, "MEC": 28.751, "H2O": 92.644, "DEC": -19.605, "TEC": 90.604},
    EtOH={"CA": 54.177, "EtOH": 0, "MEC": -79.836, "H2O": -728.971, "DEC": -135.446, "TEC": 70.5},
    MEC={"CA": -33.13, "EtOH": 17.574, "MEC": 0, "H2O": 263.187, "DEC": 28.784, "TEC": 107.843},
    H2O={"CA": 53.676, "EtOH": 756.948, "MEC": -447.773, "H2O": 0, "DEC": -336.304, "TEC": 82.56},
    DEC={"CA": -0.897, "EtOH": 60.371, "MEC": -34.376, "H2O": 178.184, "DEC": 0, "TEC": 77.814},
    TEC={"CA": -172.585, "EtOH": -301.6, "MEC": -154.776, "H2O": -501.8, "DEC": -95.921, "TEC": 0}
)

r = {"CA": 5.9585, "EtOH": 2.1055, "MEC": 7.1912, "H2O": 0.92, "DEC": 8.4674, "TEC": 9.7436}
q = {"CA": 4.808, "EtOH": 1.972, "MEC": 6.380, "H2O": 1.4, "DEC": 7.424, "TEC": 8.468}

# Equilibrium Constants (Kolah et al., 2007)
K = {"R1": 6.35, "R2": 2.72, "R3": 3.78}

# Stoichiometric number of component i in reaction Rj
coeff = {"CA": {"R1": -1, "R2": 0, "R3": 0},
         "EtOH": {"R1": -1, "R2": -1, "R3": -1},
         "MEC": {"R1": 1, "R2": -1,"R3": 0},
         "H2O": {"R1": 1, "R2": 1, "R3": 1},
         "DEC": {"R1": 0, "R2": 1, "R3": -1},
         "TEC": {"R1": 0, "R2": 0, "R3": 1},
         }