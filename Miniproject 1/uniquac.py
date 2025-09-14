def gamma(components, A, B, T, r, q, x):
    """
    Computes for the activity coefficient, gamma, based on the UNIQUAC Model.

    components: List of components in the system
    A, B: Parameters for the binary interaction parameter tau[i][j] = exp(A[i][j] + B[i][j]/T)
    T: Temperature in Kelvin.
    r: Parameter for the UNIQUAC Model, "relative molecular volume".
    q: Parameter for the UNIQUAC Model, "relative molecular surface area".
    x: Mole fraction of components.
    """

    import math

    # Model for the binary interaction parameter tau as shown in Bohorquez et al. (2020): tau[i][j] = exp(A[i][j] + B[i][j]/T)
    tau = dict()
    for i in components:
        dummy_dict = dict()
        for j in components:
            dummy_dict[j] = math.exp(A[i][j] + B[i][j] / T)
        tau[i] = dummy_dict

    # UNIQUAC Model formulation for activity coefficient (gamma) as shown in Perry's Chemical Engineers' Handbook (8th ed) p. 4-25
    theta = dict()
    for i in components:
        theta[i] = x[i] * q[i] / sum(x[j] * q[j] for j in components) # Eq. 4-285

    J = dict()
    for i in components:
        J[i] = r[i] / sum(r[j] * x[j] for j in components) # Eq 4-289

    L = dict()
    for i in components:
        L[i] = q[i] / sum(q[j] * x[j] for j in components) # Eq 4-290

    s = dict()
    for i in components:
        s[i] = sum(theta[l] * tau[l][i] for l in components) # Eq 4-291

    yR = dict()
    for i in components:
        ln_yR = q[i] * (1 - math.log(s[i]) - sum(theta[j] * tau[i][j] / s[j] for j in components)) # Eq 4-288
        yR[i] = math.exp(ln_yR)

    yC = dict()
    for i in components:
        ln_yC = 1 - J[i] + math.log(J[i]) - 5 * q[i] * (1 - J[i] / L[i] + math.log(J[i] / L[i])) # Eq 4-287
        yC[i] = math.exp(ln_yC)

    gamma = dict()
    for i in components:
        ln_gamma = math.log(yC[i]) + math.log(yR[i]) # Eq 4-286
        gamma[i] = math.exp(ln_gamma)

    return gamma