import scipy
import matplotlib.pyplot as plt
import numpy as np

components = ["CA", "EtOH", "MEC", "H2O", "DEC", "TEC"]

def batch(t, x_var, temperature, catalyst_loading):
    import uniquac
    import math
    import data

    T = temperature + 273.15  # K


    # Map x_var to x[i]. That is, x["CA'] must be mapped to x_var[0]
    x = dict()
    for index, value in enumerate(components):
        x[value] = x_var[index]

    # Compute Activities of the componentns
    gamma = uniquac.gamma(components, data.A, data.B, T, data.r, data.q, x)  # Activity Coefficients
    a = {i: gamma[i] * x[i] for i in components}  # Activities

    # Compute rate constants of reactions
    reactions = ["R1", "R2", "R3"]

    k_0_self = {"R1": 8.37e6, "R2": 9.82e6, "R3": 5e6}  # Self-catalyzed reaction pre-exponential factor in 1/(molfrac * s)
    E_A_self = {"R1": 70_800, "R2": 72_000, "R3": 72_400}  # Self-catalyzed reaction activation energy in J/mol
    k_self = {i: k_0_self[i] * math.exp(-E_A_self[i] / (8.314 * T)) for i in reactions}  # Self-catalyzed reaction rate constant in 1/(molfrac * s)

    k_0_cat = {"R1": 1.199e8, "R2": 7.105e9, "R3": 9.171e8}  # Catalyzed reaction pre-exponential factor in 1/(%wt cat * s)
    E_A_cat = {"R1": 73_324, "R2": 75_797, "R3": 88_784}  # Catalyzed reaction activation energy in J/mol
    k_cat = {i: k_0_cat[i] * math.exp(-E_A_cat[i] / (8.314 * T)) for i in reactions}  # Catalyzed reaction rate constant in 1/(%wt cat * s)

    # Define catalyst loading and concentration of active sites
    w_cat = catalyst_loading  # wt% catalyst loading
    x_acid = 3 * x["CA"] + 2 * x["MEC"] + x["DEC"]

    # Compute for the rate of each reaction
    r = dict()
    r["R1"] = (w_cat * k_cat["R1"] + x_acid * k_self["R1"]) * (a["CA"] * a["EtOH"] - a["H2O"] * a["MEC"] / data.K["R1"])  # in 1/s
    r["R2"] = (w_cat * k_cat["R2"] + x_acid * k_self["R2"]) * (a["MEC"] * a["EtOH"] - a["H2O"] * a["DEC"] / data.K["R2"])  # in 1/s
    r["R3"] = (w_cat * k_cat["R3"] + x_acid * k_self["R3"]) * (a["DEC"] * a["EtOH"] - a["H2O"] * a["TEC"] / data.K["R3"])  # in 1/s

    # Mass balances for a batch reactor
    dx_CAdt = -r["R1"]
    dx_EtOHdt = - r["R1"] - r["R2"] - r["R3"]
    dx_MECdt = r["R1"] - r["R2"]
    dx_H2Odt = r["R1"] + r["R2"] + r["R3"]
    dx_DECdt = r["R2"] - r["R3"]
    dx_TECdt = r["R3"]

    # Return system of ODEs
    return [dx_CAdt, dx_EtOHdt, dx_MECdt, dx_H2Odt, dx_DECdt, dx_TECdt]

# Reaction parameters
temperature = 78 # degrees C
catalyst_loading = 1.5 # wt% catalyst
x0 = [1/16, 15/16, 0, 0, 0, 0] # Initial Mole fractions of "CA", "EtOH", "MEC", "H2O", "DEC", "TEC"

# Bound for numerical solution
end = 10_000 # seconds

# Solve ODEs
solution = scipy.integrate.solve_ivp(batch, t_span=(0, end), y0=x0, t_eval=np.linspace(0, end, 1000), args=(temperature, catalyst_loading))

# Get solution
t = solution.t
x = {"CA": solution.y[0], "EtOH": solution.y[1], "MEC": solution.y[2], "H2O": solution.y[3], "DEC": solution.y[4], "TEC": solution.y[5]}


# Calculate Substitution Degree (SD)
SD = (x["MEC"] / 3 + x["DEC"] / 2 + x["TEC"]) / (x["CA"] + x["MEC"] + x["DEC"] + x["TEC"])

# Calculate Yield of TEC
y_TEC = x["TEC"][1:] / (x0[0] - x["CA"][1:])

# Calculate Selectivity of TEC
s_TEC = x["TEC"][1:] / (x["MEC"][1:] + x["DEC"][1:])

# Plot solution
plt.plot(t, x["CA"], label="CA")
plt.plot(t, x["EtOH"], label="EtOH")
plt.plot(t, x["MEC"], label="MEC")
plt.plot(t, x["H2O"], label="H2O")
plt.plot(t, x["DEC"], label="DEC")
plt.plot(t, x["TEC"], label="TEC")
#plt.plot(t, SD, label="SD")
#plt.plot(t[1:], y_TEC, label="TEC Yield")
#plt.plot(t[1:], s_TEC, label="TEC Selectivity")


plt.legend()
plt.show()

e1 = (x["CA"][-1] - x0[0]) / (-1)
e2 = ((x["MEC"][-1] - x0[2]) - 1 * e1) / (-1)
e3 = ((x["TEC"][-1] - x0[-1])) / (1)

print(e1, e2, e3)





