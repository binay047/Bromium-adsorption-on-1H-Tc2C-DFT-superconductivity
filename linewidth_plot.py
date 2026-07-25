import numpy as np

# Convert from GHz to meV
ghz2mev = 4.13567587265e-3

# Load elph.gamma.5.gnu file
data = np.loadtxt('elph.gamma.5.gnu')
k = np.unique(data[:, 0])
gamma = np.reshape(data[:, 1], (-1, len(k))) * ghz2mev  # convert once, in meV

# High-symmetry point positions (for reference / annotation)
gG1 = k[0]
M   = k[366]
K   = k[577]
gG2 = k[999]

# Write xmgrace-compatible .dat file: one block per band, blank line between blocks
with open('gamma_linewidth.dat', 'w') as f:
    f.write(f"# High-symmetry points: Gamma={gG1:.6f}  M={M:.6f}  K={K:.6f}  Gamma={gG2:.6f}\n")
    f.write("# Columns: k-path   linewidth(meV)\n")
    for i in range(gamma.shape[0]):
        f.write(f"# Band {i+1}\n")
        for kx, gy in zip(k, gamma[i, :]):
            f.write(f"{kx:.8f}  {gy:.8f}\n")
        f.write("\n")  # blank line = new set boundary in xmgrace

print("Wrote gamma_linewidth.dat")
print(f"High-symmetry points -> Gamma: {gG1:.4f}, M: {M:.4f}, K: {K:.4f}, Gamma: {gG2:.4f}")
