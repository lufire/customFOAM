# customFOAM
Custom solvers and libraries for OpenFOAM (Version 7.x)

Prerequisites:
- OpenFOAM 7.0 
- Optional: corresponding CFDEM library for CFDEM solvers (cfdemElectrochemicalFoam, explicitCfdemSolverPiso)

Installation instructions:
* Download or clone repository
* In 'libs/specific_libs' directory:
  ```sh
  wclean all
  wmake
  ```
* In 'solvers/specific_solver' directory:
  ```sh
  wclean all
  wmake
  ```
