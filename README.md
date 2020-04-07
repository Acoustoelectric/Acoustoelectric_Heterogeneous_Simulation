# ATI - Acoustoelectric Temporal Interference

This repository will contain the smash files from SIM4Life and python code to create a simulation of Acoustoelectric Temporal Interference. 

## Intall environment

Jean is running a Windows 10 machine(through VMWare) with Sim4LifeLight 5.2.1.1375 installed. She also has access to the full version on the lab server if need be, but it's easier to develop a basic simulation locally on her computer so she's working on that first. Once the model becomes more complex, we expect to move to full version. Hence keeping the cell count under the limit is useful currently, if only to speed up development time. 

## Model Description
For the feasibility test the following model will be used: 
<p align="center">
	<img src="images/experiment_setup.jpg" height="300">
</p>

The images folder in this repository contains the original proposal, as well as Esra's initial mathematical outline of steps to set up in the the processing pipeline.

## Pipeline
We desire a python processing pipeline which runs from acoustic simulation through to TI simulation. Esra Neufeld has provided a mathematical outline with processing steps involved, which we are now trying to code in Sim4Life/Python. 

## Specifications: 
### Step 1: Acoustic Simulation
Unipolar Waveform, 1Mhz transducer with a unipolar(initially half a sine wave) pulse which is pulsed through at a 2kHz frequency. 

Unipolar waveform = 0.5*(sin(2*pi*fs*)+|sin(2*pi*fs*t) | )

#### 








