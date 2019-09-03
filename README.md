# customFOAM
Custom solvers and libraries for OpenFOAM (Version 2.3.x)

Prerequisites:
- OpenFOAM 2.3.x
- optional: corresponding CFDEM library for particle electrode (cfdemElectrochemicalFoam, newCfdemElectrochemicalFoam)

Installation instructions:
* Download or clone repository
* In 'libs' directory:
  ```sh
  wclean all
  wmake all
  ```
* In 'solvers/specific_solver' directory:
  ```sh
  wclean all
  wmake
  ```
