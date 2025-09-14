import numpy as np
import scipy
import matplotlib.pyplot as plt
import uniquac

K = {"R1": 6.35, "R2": 2.72, "R3": 3.78}
ratios = np.linspace(1, 20, 100)
tec = []

for ratio in ratios:
    n0 = {"EtOH": ratio, "CA": 1, "MEC": 0, "H2O": 0, "DEC": 0, "TEC": 0}

    def func(X):

        xi1, xi2, xi3 = X

        n = dict()
        n["EtOH"] = n0["EtOH"] - xi1 - xi2 - xi3
        n["CA"] = n0["CA"] - xi1
        n["MEC"] = n0["MEC"] + xi1 - xi2
        n["H2O"] = n0["H2O"] + xi1 + xi2 + xi3
        n["DEC"] = n0["DEC"] + xi2 - xi3
        n["TEC"] = n0["TEC"] + xi3

        x = dict()
        x["EtOH"] = n["EtOH"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])
        x["CA"] = n["CA"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])
        x["MEC"] = n["MEC"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])
        x["H2O"] = n["H2O"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])
        x["DEC"] = n["DEC"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])
        x["TEC"] = n["TEC"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])

        eq_R1 = x["MEC"] * x["H2O"] - K["R1"] * x["CA"] * x["EtOH"]
        eq_R2 = x["DEC"] * x["H2O"] - K["R2"] * x["MEC"] * x["EtOH"]
        eq_R3 = x["TEC"] * x["H2O"] - K["R3"] * x["DEC"] * x["EtOH"]

        return [eq_R1, eq_R2, eq_R3]

    guess = [0, 0, 0]
    xi1, xi2, xi3 = scipy.optimize.fsolve(func, x0=guess)

    n = dict()
    n["EtOH"] = n0["EtOH"] - xi1 - xi2 - xi3
    n["CA"] = n0["CA"] - xi1
    n["MEC"] = n0["MEC"] + xi1 - xi2
    n["H2O"] = n0["H2O"] + xi1 + xi2 + xi3
    n["DEC"] = n0["DEC"] + xi2 - xi3
    n["TEC"] = n0["TEC"] + xi3

    x = dict()
    x["EtOH"] = n["EtOH"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])
    x["CA"] = n["CA"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])
    x["MEC"] = n["MEC"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])
    x["H2O"] = n["H2O"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])
    x["DEC"] = n["DEC"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])
    x["TEC"] = n["TEC"] / (n["EtOH"] + n["CA"] + n["MEC"] + n["H2O"] + n["DEC"] + n["TEC"])

    tec.append(x["TEC"])

plt.plot(ratios, tec)
plt.show()

