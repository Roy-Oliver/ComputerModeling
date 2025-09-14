import scipy
import data
import numpy as np
import matplotlib.pyplot as plt
import itertools

ratios = np.linspace(1, 20, 100) # EtOH: CA ratios to test
results = dict()
for temp in [80, 90, 100, 110, 120]: # Create a dictionary with temperatures as keys and TEC mole fraction as values
    results[temp] = [] # A blank list which will later contain TEC mole fraction for every EtOH: CA ratio

for temp, ratio in itertools.product(results.keys(), ratios): # Loop through all combinations of ratios and temperatures
    # Input Parameters
    n0 = {"CA": 1, "EtOH": ratio, "MEC": 0, "H2O": 0, "DEC": 0, "TEC": 0} # Initial Moles of "CA", "EtOH", "MEC", "H2O", "DEC", "TEC"
    T = temp + 273.15 # reaction temperature in K

    # Names of reactions and components
    reactions = ["R1", "R2", "R3"] # R1: prod of MEC, R2: prod of DEC, R3: prod of TEC
    components = ["CA", "EtOH", "MEC", "H2O", "DEC", "TEC"]

    def func(X):
        import uniquac

        # Unpack the unknown extents (xi) into a dictionary
        xi = dict()
        xi["R1"], xi["R2"], xi["R3"] = X

        # Calculate for the equilibrium moles and mole fractions of components using the extent
        n = {component: n0[component] + sum(data.coeff[component][reaction] * xi[reaction] for reaction in reactions) for component in components} # Uses stoichiometric numbers from "data" file
        x = {component: n[component] / sum(n[comp] for comp in components) for component in components}

        # Compute activities of the components using the UNIQUAC Model
        gamma = uniquac.gamma(components, data.A, data.B, T, data.r, data.q, x)  # Import UNIQUAC parameters from "data" file and calculate activity coefficients
        a = {component: gamma[component] * x[component] for component in components}  # Calculate activities

        # Solve expressions for equilibrium
        eq_R1 = a["MEC"] * a["H2O"] - data.K["R1"] * a["CA"] * a["EtOH"]
        eq_R2 = a["DEC"] * a["H2O"] - data.K["R2"] * a["MEC"] * a["EtOH"]
        eq_R3 = a["TEC"] * a["H2O"] - data.K["R3"] * a["DEC"] * a["EtOH"]

        return [eq_R1, eq_R2, eq_R3]

    # Solve the equilibrium equations using the extent of reaction e
    xi = dict()
    guess = [min(n0.values()), min(n0.values()), min(n0.values())] # Take as guess the minimum initial moles among the components
    xi["R1"], xi["R2"], xi["R3"] = scipy.optimize.fsolve(func, x0=guess)

    # Calculate for the equilibrium moles and mole fractions of components using the extent
    n = {component: n0[component] + sum(data.coeff[component][reaction] * xi[reaction] for reaction in reactions) for component in components}
    x = {component: n[component] / sum(n[comp] for comp in components) for component in components}

    # Assign TEC mole fraction to "results" variable
    results[temp].append(x["TEC"])

# Graph results
for temp, color in zip(results.keys(), ["skyblue", "royalblue", "black", "lightcoral", "red"]):
    plt.plot(ratios, results[temp], color, label=f"{temp} C", linewidth=1)
plt.legend() # Attach a legend
plt.xlim([1, 15]) # Set x-axis limits
plt.xticks(np.arange(1, 15, step=1)) # Set ticks on the x-axis
plt.xlabel("Initial EtOH:CA Molar Ratio") # Set xlabel
plt.ylim([0.04, None]) # Set y-axis limits
plt.ylabel("Mole Fraction TEC at Equilibrium") # Set ylabel
plt.savefig("graph.png", bbox_inches='tight', dpi=300)
plt.show()
