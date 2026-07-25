import numpy as np

# =====================================
# Input parameters from Quantum ESPRESSO
# =====================================

lambda_ep =  0.4260         # Electron-phonon coupling constant
omega_ln =  242.7609      # Logarithmic average frequency (K)

# =====================================
# McMillan equation
# =====================================

def Tc_McMillan(lam, omega_ln, mu):

    denominator = lam - mu * (1.0 + 0.62 * lam)

    if denominator <= 0:
        return 0.0

    exponent = -1.04 * (1.0 + lam) / denominator

    Tc = (omega_ln / 1.2) * np.exp(exponent)

    return Tc

# =====================================
# Generate Tc vs mu*
# =====================================

outfile = "Tc_vs_mu.dat"

with open(outfile, "w") as f:

    f.write("# mu*      Tc(K)\n")

    for mu in np.arange(0.10, 0.201, 0.01):

        Tc = Tc_McMillan(lambda_ep, omega_ln, mu)

        f.write("{:6.2f}   {:12.6f}\n".format(mu, Tc))

print("Output written to:", outfile)
