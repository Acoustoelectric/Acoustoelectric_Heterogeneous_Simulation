# ATI - Acoustoelectric Temporal Interference

This repository contains the smash files and python scripts necessary to create a simulation of Acoustoelectric Temporal Interference neuromodulation. 

The basic concept is to utilize the acoustoelectric effect, whereby an acoustic wave causes a localized change in conductivity in a medium to create an oscillating electric field, which can then be used in TI stimulation. The advantage of this is the smaller focal point of the acoustic transducer, enabling a more localized electrical stimulation. 


## Model Description
For the feasibility test the following model will be used: 
<p align="center">
	<img src="images/experiment_setup.jpg" height="300">
</p>


## Acoustoelectric Sim4Life Simulation

Jean Rintoul
Questions: jeanrintoul@imperial.ac.uk

## Structure: 
The Acoustoelectric Simulation Consists of a 4 simulation pipeline. 

1. The acoustic and lead field simulation should be run. 
2. Run acoustoelectric.py through the S4L scripter. Note that this script calls out to an external python 3 process to compute the gradient of pressure correctly. If you do not have python 3 installed on your path you will need to do this to run the script successfully. All interim quantities are stored out to NPZ files separately, to enable easy checking at every step. 
When the script has finished running, 'success!' will be written in the consol in s4l, and in the analysis tab you will see the two source terms. One source term is called 'real part diffusion', the other is 'imag part diffusion'. 
3. To run the diffusion, we compute the result for real and imaginary parts separately. In the simulation tab, in 'Image_Diffusion' simulation, import an analysis source. When you do this, you should see 'imag part diffusion' in the option for sources to import. Similarly for the 'Real_Diffusion'. Once you have ensured that both diffusion sims have the correct source term, run them.
4. Once both diffusion simulation have completed, from the scripter run 'join_ae_field'. This also calls out to an external python 3 process to compute a gradient term, and will also re-import the resultant phi and AE field into the analysis. 
5. Switch to the analysis tab and add a slice viewer to the Acoustoelectric Voltage field, and the Acoustoelectric E field(note this is the magnitude sqrt(E_x^2+ E_y^2+E_Z^2) and note a vector field).

Other files for verification and residual analysis: 
- poisson verification.py. This file takes in the resultant phi, calculates the laplacian and compares it with the source term. 
- acoustic_accuracy_test.py- this file was used to assess the accuracy of the acoustic simulation in comparison to the output of the analytic solution focused.m. The acoustic_accuracy_test_exporter can be run from within s4l once the acoustic simulation is complete, to export the results ready to plot. 
