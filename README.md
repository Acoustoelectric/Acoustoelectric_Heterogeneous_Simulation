# ATI - Acoustoelectric Temporal Interference

This repository will contain the smash files from SIM4Life and python code to create a simulation of Acoustoelectric Temporal Interference. 

The basic concept is to utilize the acoustoelectric effect, whereby an acoustic wave causes a localized change in conductivity in a medium to create an oscillating electric field, which can then be used in TI stimulation. The advantage of this is the smaller focal point of the acoustic transducer, enabling a more localized stimulation focal point. 

## Intall environment

Jean is running a Windows 10 machine(through VMWare) with Sim4LifeLight 5.2.1.1375 installed. The lab server based full version is available later when more complex models are used. 

## Model Description
For the feasibility test the following model will be used: 
<p align="center">
	<img src="images/experiment_setup.jpg" height="300">
</p>

The images folder in this repository contains the original proposal, as well as Esra's initial mathematical outline of steps to set up in the the processing pipeline.

## Goal
We desire a python processing pipeline which runs from acoustic simulation through to TI simulation. 

## Specifications: 
### Step 1: Acoustic Simulation
Currently, I've got the SEFT tutorial at transducer frequency of 55kHz, with a sinusoidal pulse. Is it possible to run through a unipolar pulse instead? So it's a 55kHz transducer modulating at 2kHz with a unipolar pulse? This would be a common set up for tumor ablation-though we are using it for a different purpose. 

Example Unipolar waveform = 0.5*(sin(2\*pi\*fs\*t)+|sin(2\*pi\*fs\*t) | )

#### 








