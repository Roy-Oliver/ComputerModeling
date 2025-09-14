import matplotlib.pyplot as plt
import numpy as np
import scipy

# Parameters
a = 21
b = 100
c = 0.81
T_opt = 21 # oC, optimal temperature

# model for temperature at time t
def temp(t):
    delta = -np.arcsin((T_opt - 21.5) * 2 / (30 - 13)) / (2 * np.pi)  # Compute for delta, Eq. 6
    T = ((30 - 13) / 2) * np.sin(2 * np.pi * (t - delta)) + (30 - 13) / 2 + 13  # Compute for temperature, Eq. 5

    return T

# Model
def model(X, t):
    T = temp(t) # Compute for the temperature at time t
    k = -((T - a) ** 2) / b + c # Compute for the kinetic constant
    dXdt = -k * X # Compute for derivative

    return dXdt

# Solve ODE numerically
# Initial conditions
X0 = 1 # mg/L

# Integration limits
start = 0 # days
end = 7 # days
t = np.linspace(start, end, 1000)

# Get solution
X = scipy.integrate.odeint(model, X0, t)

fig, ax1 = plt.subplots()
# Graph results
ax1.plot([0, 10], [0.05, 0.05], "--", color="lightcoral") # Draw a line at 95% contaminant reduction
ax1.plot(t, X[:,0], "r", label="X") # Plot X vs t
ax1.set_xlabel(r"t [$days$]") # Label axes
ax1.set_xticks(np.arange(start, end + 1, 1)) # Set x-axis ticks
ax1.set_xlim([start, end]) # Define plotting limits for x-axis
ax1.set_ylabel(r"X [$mg \cdot L^{-1}$]") # Define x-axis label
ax1.set_ylim([0, 1]) # Define plotting limits for y-axis
ax1.set_yticks([0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]) # Set y-axis ticks

ax2 = plt.twinx()
# Plot temperature
ax2.plot(t, temp(t), linewidth=1, color="lightgrey", label="T")
ax2.set_ylabel(r"Temperature [$^\circ C$]") # Define y-axis label
ax2.set_ylim([11, 32]) # Define plotting limits for y-axis
ax2.set_yticks(np.arange(11, 32, 2)) # Set y-axis ticks
fig.legend(loc="upper right", bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)

ax1.set_zorder(ax2.get_zorder()+1)
ax1.set_frame_on(False)


plt.savefig("model_t_variable.png", dpi=500) # Save the figure
plt.show()