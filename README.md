# Br_2Tc_2C-superconductivity
Here, I have calculated Tc of Bromium-adsorbed monolayers. 
pw.x < relax.in > relax. out
Now, update final atomic positions in SCF calculations from relax.out and keep them the same throughout the whole calculations. 
pw.x < scf.in > scf.out 
pw.x < dense.in > dense.out
d3hess.x < d3hess.in > d3hess.out
q2r.x < q2r.in > q2r.out 
matdyn.x < matdyn.in > matdyn.out 
plotband.x < plotband.in > plotband.out 
matdyn.x < phdos.in > phdos.out

# Now, extract phdos from matdyn.modes
awk '{print $1,$2}' matdyn.phdos > total.dat
awk '{print $1,$3}' matdyn.phdos > Br.dat
awk '{print $1,$5}' matdyn.phdos > C.dat
awk '{print $1,$6}' matdyn.phdos > Tc.dat
