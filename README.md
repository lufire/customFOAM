# customFOAM
Custom solvers and libraries for OpenFOAM (Version 2.3.x)

Prerequisits:
- OpenFOAM 2.3.x
- optional: corresponding CFDEM library for particle electrode (cfdemElectrochemicalFoam, newCfdemElectrochemicalFoam)

Installation instructions:
* Download or clone repository
* In 'libs' directory, do:
  ```sh
  wclean all
  wmake all
  ```
* Within 'solvers/specific_solver' directory, do:
  ```sh
  wclean all
  wmake
  ```
