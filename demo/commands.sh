#
# step 1: generate the solute conformation
cd step-1.solute-confomation
# if you have GROMACS 5 or above
  gmx editconf -f methylcyclopentane.pdb -o solute.gro -box 3 3 3 -center 1.5 1.5 1.5
  gmx editconf -f solute.gro -o solute.pdb
# else
  echo 'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1' > solute.pdb
  cat methylcyclopentane.pdb | egrep '(ATOM|HETATM)' | sed s/HETATM/ATOM\ \ /g >> solute.pdb
# end if
cd ..

# step 2: generate the force field with AMBER Tools and acpype 
cp step-1.solute-confomation/solute.pdb step-2.solute-ff
cd step-2.solute-ff
# if AMBER Tools are installed installed in $HOME/bin/amber18
  export AMBERHOME=$HOME/bin/amber18; export PATH=$PATH:$AMBERHOME/bin; export LD_LIBRARY_PATH=$AMBERHOME/lib:$AMBERHOME/lib64
  antechamber -i solute.pdb -fi pdb -o solute.mol2 -fo mol2 -s 2 -c bcc
# if acpype is installed in $HOME/bin/acpype-master
  $HOME/bin/acpype-master/acpype.py -i solute.mol2 -c user
# finally generate solute.ff
  gmxtop2solute -top solute.acpype/solute_GMX.top -o solute.ff
cd ..

# step 3: run rismhi3d
cp step-1.solute-confomation/solute.pdb step-3.run-rismhi3d
cp step-2.solute-ff/solute.ff step-3.run-rismhi3d
export IETLIB=`pwd`/solvent
cd step-3.run-rismhi3d
