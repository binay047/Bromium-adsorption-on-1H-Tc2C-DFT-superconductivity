import numpy as np
file_path = "/home/binay/Desktop/DFT/Superconductivity/Br_1H_Tc_2C/a2F.dos10" # Here, you need to paste the path to a2F.dos10 and a2F.dos10 is chosen using the converged value of Tc calculated from lambda.in
freq, a2F = [], []
with open(file_path, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) > 0 and cols[0] == "lambda":
            break
        if len(cols) >= 2:
            try:
                freq.append(float(cols[0]))
                a2F.append(float(cols[1]))
            except:
                pass

freq = np.array(freq)
a2F = np.array(a2F)
mask = freq > 0
freq = freq[mask]
a2F = a2F[mask]
idx = np.argsort(freq)
freq = freq[idx]
a2F = a2F[idx]

RY_TO_CM = 109737.31568539
freq_cm = freq * RY_TO_CM

lam = np.zeros_like(freq)
for i in range(1, len(freq)):
    dw = freq[i] - freq[i-1]
    y1 = 2.0 * a2F[i-1] / max(freq[i-1], 1e-12)
    y2 = 2.0 * a2F[i]   / max(freq[i], 1e-12)
    lam[i] = lam[i-1] + 0.5 * (y1 + y2) * dw

print("Total lambda =", lam[-1])

np.savetxt("a2F_combined.dat",
           np.column_stack([freq_cm, a2F, lam]),
           header="freq_cm-1  a2F  lambda",
           comments="# ")
