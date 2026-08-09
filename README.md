# Br-adsorbed Tc₂C superconductivity workflow (Quantum ESPRESSO)

This guide explains the complete workflow used to calculate the superconducting critical temperature (Tc) of a bromine (Br) adsorbed Tc₂C monolayer.

---

## 1. Build adsorption structures

### Theory
Adatom adsorption on a 2D monolayer can occupy several inequivalent high-symmetry sites — directly above a metal atom ("top"), above the non-metal atom, or in a "bridge" position between two neighboring atoms — and the energetically preferred site is not known a priori; it must be determined by comparing total energies. The most stable configuration is the global minimum of

$$
E_{\text{ads}} = E_{\text{total}}(\text{site})
$$

across the candidate sites, evaluated at fixed (unrelaxed or lightly relaxed) geometry via a quick SCF pass. The magnetization check (total and absolute magnetization $\to 0$) is a physical consistency requirement here: Tc₂C-based systems are expected to be nonmagnetic in the correct ground state for this study, so a candidate structure converging to a finite magnetic moment signals either the wrong adsorption site or an unphysical electronic configuration, not just a higher-energy one.

### Procedure
Create three structures by placing the Br atom at different adsorption sites:
* **Top of Tc** → Br has the same x and y coordinates as a Tc atom.
* **Top of C** → Br has the same x and y coordinates as the C atom.
* **Bridge site** → Br is placed midway between two Tc atoms.
* Run a quick SCF calculation for each structure and compare the total energies.
* pw.x < scf_topTc.in > scf_topTc.out
* pw.x < scf_topC.in  > scf_topC.out
* pw.x < scf_bridge.in > scf_bridge.out
* Choose the structure with the **lowest total energy**, and also total and absolute magnetisation need to be zero.

---

## 2. Convergence Test

### Theory
As with any plane-wave DFT calculation, results are only meaningful once three numerical parameters — the plane-wave cutoff, the k-point sampling density, and the in-plane lattice constant — are converged, since all downstream quantities (density, energy, phonons, electron-phonon coupling) inherit any residual error from these choices. Each is optimized sequentially and fed forward into the next: a converged $E_{\text{cut}}$ fixes the basis for the k-point scan, a converged k-mesh fixes the sampling for the lattice/volume scan, and the resulting equilibrium lattice constant $a_0$ (from fitting $E_{\text{tot}}(a)$ to an equation of state via `ev.x`) is then used in a variable-cell relaxation (`vc-relax`) to obtain the fully relaxed cell geometry and atomic positions used in every subsequent step.

### 2.A. Optimisation of planewave (ecutwfc)
* chmod +x ecut.sh
* ./ecut.sh
* Take an ecut that is converged and use that in the input file for optimisation of k-points

### 2.B. Optimisation of k-points
* chmod +x kpoi.t.sh
* ./kpoint.sh
* Take a converged and use that in the input file for optimisation of lattice

### 2.C. Optimisation of lattice
* chmod +x lattice.sh
* ./lattice.sh
* In the same lattice directory, in the terminal enter ev.x
* ang
* noncubic
* 4
* input file name: lattice.dat
* output file name: bin
* take this a0 from bin into vc_relax.in

### 2.D. Do vc-relax
* Run variable-cell relaxation:
* pw.x < vc_relax.in > vc_relax.out
* Copy the following from "vc_relax. out":
* "CELL_PARAMETERS"
* "ATOMIC_POSITIONS"
* Note: BFGS must be converged, and to find relaxed  cell_parameters and atomic_positions, press Ctrl + F inside vc_relax.out and type 'final bfgs'
* Paste them into:
* "scf. in" that is inside pseudo.sh
* Use these optimised positions for all remaining calculations.

### 2.E. Optimisation of pseudopotential
* chmod

### 2.F. Optimisation of degauss value

#### Theory
For metallic (or narrow-gap) 2D systems, partial occupations near $E_F$ are handled via smearing, and the choice of smearing scheme and width (`degauss`) can shift total energies and derived quantities like the electron-phonon coupling. This step scans several smearing types (Fermi-Dirac "fd", Gaussian "gauss", Methfessel-Paxton "mp", Marzari-Vanderbilt "mv") and widths, so a value can be chosen where the total energy is stable with respect to further reduction of `degauss` — since the McMillan/Eliashberg electron-phonon calculations in later steps are directly sensitive to this choice through the double-delta-function sum at $E_F$.

* Here, run the following two lines of code
* chmod +x degauss.sh
* ./degauss.sh
* After that, the smearing_results.dat file will be created, and use the awk command
* awk '$1=="fd" {print $2, $3}' smearing_results.dat > fd.dat
* awk '$1=="gauss" {print $2, $3}' smearing_results.dat > gauss.dat
* awk '$1=="mp" {print $2, $3}' smearing_results.dat > mp.dat
* awk '$1=="mv" {print $2, $3}' smearing_results.dat > mv.dat

---

## 3. SCF Calculation

### Theory
The self-consistent field calculation solves the Kohn-Sham equations to self-consistency on the relaxed structure, producing the converged ground-state charge density $n(\mathbf{r})$ that all subsequent steps (phonons, electron-phonon coupling, DOS) depend on. Using identical MPI core counts across the workflow matters specifically for QE's recover/restart mode, since parallelization-dependent data (k-point and PW pool distribution) written to the scratch directory must match on restart.

* Use the **same number of MPI cores** for all calculations below if you are doing calculations using recover mode.
* Example with 8 cores:
* mpirun -np 8 pw.x < scf.in > scf.out
* mpirun -np 8 pw.x <dense.in> dense.out

---

## 4. Dense Electronic Calculation

### Theory
A non-self-consistent calculation on a denser k-mesh than the SCF run is required for an accurate evaluation of the double-delta-function surface integral at the Fermi level, which underlies the electron-phonon coupling strength $\lambda_{\mathbf{q}\nu}$ computed in later steps — a sparse mesh under-resolves the Fermi surface and produces noisy or inaccurate $\lambda$ values.

---

## 5. D3 Hessian Calculation (only if D3 correction is used)

### Theory
When Grimme's DFT-D3 dispersion correction is included in the SCF Hamiltonian, its contribution to the dynamical matrix (Hessian of the total energy, including the D3 pairwise-correction term) must be computed separately and added to the DFPT phonon calculation, since standard `ph.x` DFPT does not natively differentiate the semi-empirical D3 energy term. Omitting this step when D3 is active would produce a dynamical matrix inconsistent with the SCF Hamiltonian actually used.

* If scf. in contains:
* vdw_corr = 'grimme-d3'
* dftd3_threebody = .false.
* run:
* mpirun -np 8 d3hess.x < d3hess.in > d3hess.out
* Otherwise, skip this step.

---

## 6. Phonon Calculation

### Theory
Using Density Functional Perturbation Theory, `ph.x` computes the dynamical matrix $D(\mathbf{q})$ at each q-point in the chosen mesh directly from the linear response of the self-consistent density to atomic displacements, without the need for finite supercells. This is the foundational step for everything electron-phonon and superconductivity related in this workflow: the phonon frequencies, eigenvectors, and (via the `elph` machinery invoked here) the electron-phonon matrix elements are all obtained from this single DFPT run.

* mpirun -np 8 ph.x < ph.in > ph.out

---

## 7. Convert Dynamical Matrices

### Theory
`q2r.x` inverse-Fourier-transforms the dynamical matrices $D(\mathbf{q})$, known only on the coarse DFPT q-mesh, into real-space interatomic force constants $C(\mathbf{R})$. This real-space representation can then be interpolated back onto an arbitrarily dense q-path or q-mesh (Step 8, 10), which is far cheaper than running DFPT directly on a dense mesh.

* q2r.x < q2r.in > q2r.out

**Purpose:** Convert phonon data to real-space force constants.

---

## 8. Phonon Band Structure

### Theory
`matdyn.x` Fourier-interpolates the real-space force constants back onto a dense path of q-points through the Brillouin zone to produce the phonon dispersion $\omega_s(\mathbf{q})$, and `plotband.x` formats this for plotting. Checking for negative (imaginary) frequencies anywhere along the path is the standard test of dynamical stability — an imaginary branch indicates the relaxed structure sits at a saddle point of the potential energy surface rather than a true local minimum, which invalidates any electron-phonon/Tc result computed on that structure.

* matdyn.x < matdyn.in > matdyn.out
* plotband.x < plotband.in > plotband.out
* Check the generated *.frq.gp files for negative frequencies.

**Purpose:** Generate phonon dispersion data for plotting.

---

## 9. If Negative Phonon Frequencies Appear

### Theory
An imaginary phonon branch in a 2D monolayer very often originates from residual in-plane stress left over from an imperfect relaxation (common in 2D systems due to the interplay of the vacuum region and the fixed out-of-plane lattice vector). Applying a small in-plane (uniaxial) strain and re-relaxing perturbs the structure enough to let BFGS find a genuinely stable local minimum, which is why the strain is applied only to the in-plane lattice vectors (rows 1 and 2) while the vacuum-containing row is left untouched — straining the vacuum direction has no physical meaning for a 2D sheet.

After phonon calculations, check the *.frq.gp files.
* If all frequencies are **positive**, continue.
* If any frequency is **negative**, apply a small strain and relax again.

### Compressive strain
* −1%, −2%, or −3%

### Tensile strain
* +1%, +2%, or +3%
* Note: apply uniaxial strain, leaving the last row of cell_parameters untouched
* For example: to apply 1% tensile strain, you need to multiply the 1st and 2nd rows of cell_parameters by 1.01, and for compressive 1% by 0.99
* Run relaxation again:
* pw.x < relax.in > relax. out
* Update the atomic positions from relax.out and use them in every later step from step 3 to step 9.

---

## 10. Phonon Density of States (PhDOS)

### Theory
The phonon DOS $g(\omega)$, obtained the same way as the dispersion but integrated over a dense q-mesh rather than a path, is required as the vibrational input to the McMillan-Allen-Dynes formula for $T_c$ (via the Eliashberg spectral function in later steps) and lets low- versus high-frequency contributions be attributed to specific atomic species. Splitting into atomic contributions (Br, C, Tc) is done by summing the appropriate atom-resolved columns of the `matdyn.phdos` output, matching each atom index to its species in the unit cell.

* matdyn.x < phdos.in > phdos.out
* Extract total and atomic contributions:
* awk '{print $1,$2}' matdyn.phdos > total.dat
* awk '{print $1, $3+$4}' matdyn.phdos > Br.dat   # atoms 1+2 = Br
* awk '{print $1, $5}'     matdyn.phdos > C.dat    # atom 3 = C
* awk '{print $1, $6+$7}' matdyn.phdos > Tc.dat   # atoms 4+5 = Tc

**Purpose:** Obtain total and atom-projected phonon DOS.

---

## 11. Prepare Electron–Phonon Coupling Input for lambda.in

### Theory
The mode- and q-point-resolved electron-phonon coupling strengths $\lambda_{\mathbf{q}\nu}$, computed by `ph.x` during Step 6 and written to the `elph_dir/elph.inp_lambda.*` files, must be assembled with their q-point weights (from the irreducible q-point list QE prints in `ph.out`) into a single `lambda.in` file. `lambda.x` then performs the Brillouin-zone sum

$$
\lambda = \sum_{\mathbf{q}\nu} \lambda_{\mathbf{q}\nu}\, w_{\mathbf{q}}
$$

to obtain the total isotropic electron-phonon coupling constant $\lambda$, and combines it with the Coulomb pseudopotential $\mu^*$ (empirical parameter capturing residual Coulomb repulsion between Cooper-pair electrons, screened by the phonon-mediated retardation) in the McMillan-Allen-Dynes equation to estimate $T_c$.

* Open ph. out and find the q-points listed for dyn.*.
* Example:
* 16 0.02 0
* 4
* 0.000000000   0.000000000   0.000000000   1
* 0.000000000   0.288675135   0.000000000   6
* 0.000000000  -0.577350269   0.000000000   3
* 0.250000000   0.433012702   0.000000000   6
* elph_dir/elph.inp_lambda.1
* elph_dir/elph.inp_lambda.2
* elph_dir/elph.inp_lambda.3
* elph_dir/elph.inp_lambda.4
* 0.10
* 16 → total number of q-points.
* 0.02 → degauss value used in ph. in.
* 4 → number of irreducible q-points.
* The next four lines are q-point coordinates and weights.
* elph.inp_lambda.* → electron–phonon data files.
* 0.10 → Coulomb pseudopotential (μ*).
* lambda.x <lambda.in> lambda.out
* Note: if there is any minus sign in elph. inp_lambda.* lines before Gauss broadening, then remove the negative sign

---

## 12. Generate α²F(ω) and λ(ω)

### Theory
The Eliashberg spectral function $\alpha^2F(\omega)$ encodes the frequency-resolved strength of electron-phonon coupling and is the central input to the Allen-Dynes formula for $T_c$:

$$
\lambda = 2\int_0^\infty \frac{\alpha^2F(\omega)}{\omega}\, d\omega
$$

The `a2F.dos*` files QE produces correspond to different smearing widths (`degauss` values) used in the double-delta sum; the specific file chosen (`a2F.dos10` here) is the one whose resulting $\lambda$/$T_c$ was found converged with respect to smearing width in Step 2.F, so that the reported $T_c$ is not an artifact of an under- or over-broadened Fermi-surface sum.

* Run:
* Inside extract_a2F_lambda.py, in file_path, you need to paste the path to a2F.dos10, and a2F.dos10 is chosen using the converged value of Tc calculated from lambda.in
* python3 extract_a2F_lambda.py
* plot. dat using xmgrace

---

## 13. Plot Phonon Linewidth

### Theory
The phonon linewidth $\gamma_{\mathbf{q}\nu}$ is directly proportional to the mode-resolved electron-phonon coupling $\lambda_{\mathbf{q}\nu}\,\omega_{\mathbf{q}\nu}$ — phonon modes that couple strongly to the electronic states at $E_F$ are damped (broadened) more strongly by their interaction with the electron gas. Plotting linewidth against the phonon dispersion identifies which specific vibrational branches and q-points drive the total coupling $\lambda$, which is useful for physically interpreting which bonds/motions are responsible for the superconducting pairing.

* plotband.x <linewidth.in> linewidth.out
* and run python3 linewidth_plot.py
* plot elph.gamma.5.gnu using xmgrace

---

## 13. Fermi Surface Plot

### Theory
The Fermi surface — the constant-energy surface $\varepsilon_n(\mathbf{k}) = E_F$ in reciprocal space — determines which electronic states participate in the double-delta-function sum underlying $\lambda_{\mathbf{q}\nu}$; its shape (nesting features, sheet multiplicity, curvature) can directly explain why certain phonon q-vectors couple more strongly than others in Step 13's linewidth plot. `plot_fermi_surface.py` visualizes this surface from the `.bxsf` file QE writes containing $\varepsilon_n(\mathbf{k})$ on the full k-mesh.

* python3 plot_fermi_surface.py aiida_fs.bxsf -o fermi_surface_Tc2CBr2.png
