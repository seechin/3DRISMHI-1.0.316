# 3DRISMHI

The original release of 3DRISM-HI

Alpha version 0.246.1472, (c) Cao Siqin, 2015-2019

# Introduction

The 3DRISMHI is an open-source software package to calculate the atomic distributions of solvents around given solutes based on integration equation theory (IET) or liquid and hydrophobicity induced density inhomogeneity (HI) theory [1–3]. Based on the structure of solutes, the force field parameters of solute, the solvent-solvent correlations (provided in the package or generated from MD simulations), the molecular structure of solvents (provided in the package or generated from structures) and some experimental parameters (provided in the package or modify youself), the 3DRISMHI can efficiently calculate the solvent density distrubution.

The 3DRISMHI can use multi cores (of one computer) to speed up the calculation. The memory required to perform a 3DRISMHI calculation is probably huge, so by default the allocated memory will not exceed the physical RAM installed on your computer (you can bypass this check via -ignore-memory-capacity). 3DRISMHI doesn’t provide MPI or GPU paralleling due to highly and frequently exchange of memories.

The 3DRISMHI is an independent software that can run on Linux, Mac OS, Windows Subsystem Linux and Windows Cygwin. The only required package is FFTW version 3. More features will be avaiable if compiled with the following packages: libz, GROMACS. The 3DRISMHI doesn’t have any graphical interface, and need to run in terminals or computer clusters.

# References

[1] Siqin Cao, Fu Kit Sheong, and Xuhui Huang, “Reference interaction site model with hydrophobicity induced density inhomogeneity: An analytical theory to compute solvation properties of large hydrophobic solutes in the mixture of polyatomic solvent molecules”, Journal of Chemical Physics 143, 054110 (2015)

[2] Siqin Cao, Lizhe Zhu and Xuhui Huang, 3DRISM-HI-D2MSA: an improved analytic theory to compute solvent structure around hydrophobic solutes with proper treatment of solute–solvent electrostatic interactions, Mol. Phys. 116, 1003 (2017)

[3] Siqin Cao, Kirill Konovalov, Ilona Christy Unarta and Xuhui Huang, Recent Developments in Integral Equation Theory for Solvation to Treat Density Inhomogeneity at Solute-Solvent Interface, Advanced Theory and Simulations (2019) DOI: 10.1002/adts.201900049
