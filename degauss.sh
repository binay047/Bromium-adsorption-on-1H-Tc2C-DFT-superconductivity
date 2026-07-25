#!/bin/bash

# Smearing methods
for sf in fd gauss mp mv
do

# Smearing values (Ry)
for se in 0.005 0.010 0.015 0.020 0.025 0.030 0.035 0.040
do

cat > ${sf}_${se}.in << EOF
&CONTROL
  calculation = 'scf'
  etot_conv_thr =   5.0000000000d-06
  forc_conv_thr =   1.0000000000d-05
  outdir = './out/'
  prefix = 'aiida'
  pseudo_dir = './pseudo/'
  verbosity = 'high'
/
&SYSTEM
  ecutrho =  600
  ecutwfc =  60
  ibrav = 0
  nat = 5
  nbnd = 40
  nosym = .false.
  ntyp = 3
  smearing='$sf'
  degauss=$se
  occupations = 'smearing'
  vdw_corr                  = "grimme-d3"
  dftd3_threebody           = .false.
  assume_isolated = '2D'
/
&ELECTRONS
  conv_thr =   1.0000000000d-10
  electron_maxstep = 300
  mixing_beta =   4.0000000000d-01
/
ATOMIC_SPECIES
Br     79.904   Br.pbe-dn-rrkjus_psl.1.0.0.UPF
C      12.0107  C.pbe-n-rrkjus_psl.1.0.0.UPF
Tc     98.0     Tc.pbe-spn-rrkjus_psl.0.3.0.UPF
ATOMIC_POSITIONS (crystal)
Br               0.3333332865        0.6666665730        0.6154058876
Br               0.3333332865        0.6666665730        0.3845945370
C                0.0000000000        0.0000000000        0.4999999987
Tc               0.3333332865        0.6666665730        0.4594255155
Tc               0.3333332865        0.6666665730        0.5405744514
K_POINTS automatic
24 24 1 0 0 0
CELL_PARAMETERS angstrom
      2.9426199818       0.0000000000       0.0000000000  
     -1.4713099909       2.5483836579       0.0000000000  
      0.0000000000       0.0000000000      32.7809600830  

EOF

echo "Running ${sf}  degauss=${se}"

mpirun -np 2 pw.x < ${sf}_${se}.in > ${sf}_${se}.out

energy=$(grep "!" ${sf}_${se}.out | tail -1 | awk '{print $5}')

printf "%-6s %8s %20s\n" "$sf" "$se" "$energy" >> smearing_results.dat

done
done

echo "Finished!"
