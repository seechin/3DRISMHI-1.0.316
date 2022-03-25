# 3DRISMHI

The original release of 3DRISM-HI

3DRISMHI 1.0.313, (c) 2022 Cao Siqin

# Licence

3DRISM-HI is free software for non-commercial use. You can use, modify or redistribute for non-commercial purposes under terms of the GNU General Public (GPL) License version 3 as published by the Free Software Foundation. A GNU GPL3 Licence is attached to this program.

3DRISM-HI uses FFTW 3 for fast Fourier transformations [1].

# Introduction

The 3DRISMHI is an open-source software to calculate the atomic distributions of solvents around given solutes based on integration equation theory (IET) or liquid and hydrophobicity induced density inhomogeneity (HI) theory [2-4]. With the structure of solutes, the force field parameters of solute, the solvent-solvent correlations (provided in the package or generated from MD simulations), the molecular structure of solvents (provided in the package or generated from trajectories) and some experimental parameters (provided in the package), the 3DRISMHI can efficiently calculate the solvent density distrubution.

3DRISMHI is compatible to GROMACS and AMBER.

The 3DRISMHI can run on Linux, Mac OS, Windows Subsystem Linux and Windows Cygwin. FFTW version 3 is required. Additional libraries can be used: libz (for compression of output data), GROMACS (for reading XTC trajectories). The graphic interface of 3DRISMHI is provided by CWBSol.

The 3DRISMHI can use multi cores (on the same computer) to speed up the calculation. 3DRISMHI doesn't provide MPI or GPU paralleling. The memory required to perform a 3DRISMHI calculation can be huge, so a large RAM is expected if the grid number is huge.

# References

[1] Matteo Frigo and Steven G. Johnson, “The design and implementation of FFTW3,” Proc. IEEE 93 (2), 216–231 (2005)

[2] Siqin Cao, Fu Kit Sheong, and Xuhui Huang, “Reference interaction site model with hydrophobicity induced density inhomogeneity: An analytical theory to compute solvation properties of large hydrophobic solutes in the mixture of polyatomic solvent molecules”, Journal of Chemical Physics 143, 054110 (2015)

[3] Siqin Cao, Lizhe Zhu and Xuhui Huang, 3DRISM-HI-D2MSA: an improved analytic theory to compute solvent structure around hydrophobic solutes with proper treatment of solute–solvent electrostatic interactions, Mol. Phys. 116, 1003 (2017)

[4] Siqin Cao, Kirill Konovalov, Ilona Christy Unarta and Xuhui Huang, Recent Developments in Integral Equation Theory for Solvation to Treat Density Inhomogeneity at Solute-Solvent Interface, Advanced Theory and Simulations (2019) DOI: 10.1002/adts.201900049
