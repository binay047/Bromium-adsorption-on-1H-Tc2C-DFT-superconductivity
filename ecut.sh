#bin/sh/     
                                                                                                                                               
NAME="Ecutwfc"
for ecut in 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 ; do
cat > "${NAME}_${ecut}.in" <<EOF 

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
  degauss =   2.0000000000d-02
  ecutrho =   $((10*ecut))
  ecutwfc =   $ecut
  ibrav = 0
  nat = 5
  nbnd = 40
  nosym = .false.
  ntyp = 3
  occupations = 'smearing'
  smearing                  = "methfessel-paxton"
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
  mpirun -np 2 pw.x -inp "${NAME}_${ecut}.in" | tee "${NAME}_${ecut}.out"   
#   pw.x <"${NAME}_${ecut}.in"> "${NAME}_${ecut}.out"
    echo "${NAME}_${ecut}"
    grep "!" "${NAME}_${ecut}.out"

 # Write cut-off and total energies in calcecut.dat.                       
                                                                                  
    awk '/!/ {printf "%d %s\n", ('$ecut'), $5}' "${NAME}_${ecut}.out" >>"ecutwfc.dat"  
done

