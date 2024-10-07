# Polarizability tensor calculation with LAMMPS

## Goal
The goal of this code is to evaluate the polarizability tensor $\alpha$ of one substitutional solute atom.

## Polarizability tensor

The polarizability tensor can be determined by measuring the change of the elastic dipole when an external strain is applied and performing a linear regression. The elastic can be determined knowing the volume V of the box and the average residual stresses induced by the solute atom over the box volume :
$$P_{ij}^0 = -V \langle \sigma_{ij}\rangle$$

$$P_{ij}=P_{ij}^0 + \alpha_{ijkl}\varepsilon_{kl}^{ext}$$

Where $P_{ij}^0$ is the permanent elastic dipole when no strain is applied.

## How to use

1. Change the elastic constant of the solvant and the potential-related lines in the "Polarizability.lammps" code and launch the calculation.
2. Many text files are generated, corresponding to the $P_{ij}$ values when strains are applied.
3. Launch the "Interpolation-post-treatment.py" file with Python. It should be in the same folder than the results of the previous LAMMPS calculcation. It will compute all components of the polarizability tensor in Voigt notation and write them in a .txt file, or in the "alpha" variable for more digits. Curves of the linear regression can be found in a sub-folder "figures" created.

