import matplotlib.pyplot as plt
import numpy as np

# Parameters
a = 21
b = 100
c = 0.81

# Create an array of concentrations at 1, 2, and 3 mg/L
X = [1, 2, 3] # mg/L

# Create an array of temperatures from 0 to 50 oC
T = np.linspace(0, 50, 1000)

# Calculate k and dXdt
dXdt = dict()
for Xi in X: # Calculate dXdt for each concentration and temperature
    k = - ((T - a) ** 2) / b + c # Rate constant in 1/day
    dXdt[Xi] = -k * Xi # Rate of degradation of toxic substance, mg/L*day


# Plot results
plt.plot(T, dXdt[1], "r", label="X = 1 mg/L", )
plt.plot(T, dXdt[2], "g", label="X = 2 mg/L", )
plt.plot(T, dXdt[3], "k", label="X = 3 mg/L", )
plt.xlabel(r"Temperature [$^\circ C$]")
plt.ylabel(r"dX/dt [$mg \cdot L^{-1} \cdot day^{-1}$]")
plt.legend()
plt.savefig("t_opt.png", dpi=500) # Save the figure
plt.show() # A minimum around 21oC for all concentrations


# Find and print temp at minimum dXdt. Equivalently, this is the temp at maximum rate of degradation
min_rate = min(dXdt[1]) # Get minimum rate
T_opt = T[np.where(dXdt[1] == min_rate)] # Find index of minimum rate. Index should be the same as optimal temp.
print(T_opt) # Minimum is at 21.02 oC