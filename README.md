# Br_2Tc_2C-superconductivity
! Here, I have calculated Tc of Bromium-adsorbed Tc2C monolayers. 
! At first, run scf. in keeping Bromine above Tc, above Carbon, and in the middle of two Tc; whichever case has the lowest energy, take that configuration of atomic position and proceed to vc-relax calculation and reach up to phonon calculation.
! If you get any negative frequency in xxxx.frq.gp files, return to the initial calculation and apply either compressive or tensile strain (-3% to +3%), and do a relax calculation
pw.x < relax.in > relax. out
Now, update the final atomic positions in the SCF calculations from relax.out and keep them the same throughout the entire calculation. 
!For below (scf.in, dense. in, d3hess.in, ph.in), please use the same number of cores (n) in mpirun -np n; otherwise, you will encounter a crash. 
pw.x < scf.in > scf.out 
pw.x < dense.in > dense.out
d3hess.x < d3hess.in > d3hess.out
ph.x <ph.in> ph.out
q2r.x < q2r.in > q2r.out 
matdyn.x < matdyn.in > matdyn.out 
plotband.x < plotband.in > plotband.out 
matdyn.x < phdos.in > phdos.out

# Now, extract phdos from matdyn.modes
awk '{print $1,$2}' matdyn.phdos > total.dat
awk '{print $1,$3}' matdyn.phdos > Br.dat
awk '{print $1,$5}' matdyn.phdos > C.dat
awk '{print $1,$6}' matdyn.phdos > Tc.dat

# For extracting Critical temperature
! from ph. out, find all initial q from dyn.1, dyn.2, dyn.3 and dyn.5 and below (1,6,3,6) are no of q's in the given q of the respective dyn
! 16 is written from (1+6+3+6), 0.02 degauss used in ph. in, elph_dir/elph.inp_lamdaa.* is listed in elph_dir folder and 0.10 is same for all calculations
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

# Now plotting Eliasberg spectral function and electron-coupling constant
python3 extract_a2F_lambda.py
! For plotting it using xmgrace, you can watch here: www.youtube.com/@BinayLimbu-bk7yw

