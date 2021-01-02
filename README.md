# **Laser cooling simulation**

## **About the project**

This project is a python based numercial simulation program to study solid phase laser cooling. This program employs foward-time center-space finite element approach (explicit Runge Kutta method) to solve charge carrier rate equations as well as heat diffusion equation. The former one includes photoinduced charge carrier generation, three types of charge carrier recombinations (trap states assisted first order recombination, bimolecular radiative recombination, Auger recombination),  charge carrier drift-diffusion process. The latter one includes charge carrier recombination induced cooling/heating as well as heat diffusion process. As an illustration, the simulation laser cooling process within a GaAs disk (radius = 9 um, thickness = 0.1um) is shown here and results can be found in https://aip.scitation.org/doi/10.1063/1.5049376. More laser cooling material and configuration can be found in https://aip.scitation.org/doi/10.1063/1.5049376 and https://www.nature.com/articles/s41427-019-0156-4. 


## **Prerequisites**

* Material parameters required in the simualation, more detail can be found in "config.py"
* Python 2.7. Other required packages can be found in "requirements.txt". 


## **Usage**

1. Modify the parameters "inconfig.py" accoding to the material and dimension configuration.
2. Execute "simulation.py". If "showprogress" is set to be "True", transient charge carrier density distribution, net charge density distribution, electric field distribution  will show up during the simuation. If "save" is set to be "True", data will be saved automatically.

## **Issues**

1. Many material parameters need to be known ahead. Some parameters changes for different samples.
2. It's complicated to change the shape of cooling material or add other components. 


## **Liscence**

Distributed under the MIT License. See `LICENSE` for more information.



