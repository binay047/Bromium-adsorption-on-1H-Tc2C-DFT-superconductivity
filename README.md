# Br-adsorbed Tc₂C superconductivity workflow (Quantum ESPRESSO)
This guide explains the complete workflow used to calculate the superconducting critical temperature (Tc) of a bromine (Br) adsorbed Tc₂C monolayer.

## 1. Build adsorption structures
Create three structures by placing the Br atom at different adsorption sites:
* **Top of Tc** → Br has the same x and y coordinates as a Tc atom.
* **Top of C** → Br has the same x and y coordinates as the C atom.
* **Bridge site** → Br is placed midway between two Tc atoms.
Run a quick SCF calculation for each structure and compare the total energies.
pw.x < scf_topTc.in > scf_topTc.out
pw.x < scf_topC.in  > scf_topC.out
pw.x < scf_bridge.in > scf_bridge.out
Choose the structure with the **lowest total energy**, and also tot and abs magnetisation need to be zero.

## 2. Optimise the structure (vc-relax)
Run variable-cell relaxation:
pw.x < vc_relax.in > vc_relax.out
Copy the following from "vc_relax. out":
* "CELL_PARAMETERS"
* "ATOMIC_POSITIONS"
Note: BFGS must be converged, and to find relaxed  cell_parameters and atomic_positions, press Ctrl + F inside vc_relax.out and type 'final bfgs'
Paste them into:
* "scf.in"
* "dense_scf.in"
Use these optimised positions for all remaining calculations.

## 3. If negative phonon frequencies appear
After phonon calculations, check the *.frq.gp files.
* If all frequencies are **positive**, continue.
* If any frequency is **negative**, apply a small strain and relax again.

### Compressive strain
* −1%, −2%, or −3%
### Tensile strain
* +1%, +2%, or +3%
Note: apply uniaxial strain, leaving the last row of cell_parameters untouched 
For example: to apply 1% tensile strain, you need to multiply the 1st and 2nd rows of cell_parameters by 1.01, and for compressive 1% by 0.99
Run relaxation again:
pw.x < relax.in > relax. out
Update the atomic positions from relax.out and use them in every later step.

## 4. SCF calculation

Use the **same number of MPI cores** for all calculations below.

Example with 8 cores:
mpirun -np 8 pw.x < scf.in > scf.out

## 5. Dense electronic calculation

## 6. D3 Hessian calculation (only if D3 correction is used)
If scf. in contains:
vdw_corr = 'grimme-d3'
dftd3_threebody = .false.
run:
mpirun -np 8 d3hess.x < d3hess.in > d3hess.out
Otherwise, skip this step.

## 7. Phonon calculation
mpirun -np 8 ph.x < ph.in > ph.out

## 8. Convert dynamical matrices
q2r.x < q2r.in > q2r.out
**Purpose:** Convert phonon data to real-space force constants.

## 9. Phonon band structure
matdyn.x < matdyn.in > matdyn.out
plotband.x < plotband.in > plotband.out
**Purpose:** Generate phonon dispersion data for plotting.
Check the generated *.frq.gp files for negative frequencies.

## 10. Phonon density of states (PhDOS)
matdyn.x < phdos.in > phdos.out
Extract total and atomic contributions:
awk '{print $1,$2}' matdyn.phdos > total.dat
awk '{print $1,$3}' matdyn.phdos > Br.dat
awk '{print $1,$5}' matdyn.phdos > C.dat
awk '{print $1,$6}' matdyn.phdos > Tc.dat
**Purpose:** Obtain total and atom-projected phonon DOS.

## 11. Prepare electron–phonon coupling input for lambda.in
Open ph. out and find the q-points listed for dyn.*.
Example:
16 0.02 0
4
0.000000000   0.000000000   0.000000000   1
0.000000000   0.288675135   0.000000000   6
0.000000000  -0.577350269   0.000000000   3
0.250000000   0.433012702   0.000000000   6
elph_dir/elph.inp_lambda.1
elph_dir/elph.inp_lambda.2
elph_dir/elph.inp_lambda.3
elph_dir/elph.inp_lambda.4
0.10
* 16 → total number of q-points.
* 0.02 → degauss value used in ph. in.
* 4 → number of irreducible q-points.
* The next four lines are q-point coordinates and weights.
* elph.inp_lambda.* → electron–phonon data files.
* 0.10 → Coulomb pseudopotential (μ*).
lambda.x <lambda.in> lambda.out

## 12. Generate α²F(ω) and λ(ω)

Run:
python3 extract_a2F_lambda.py
plot. dat using xmgrace



