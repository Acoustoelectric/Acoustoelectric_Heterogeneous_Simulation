# 
# Unipolar Equation Modelling
#
import numpy as np 
import matplotlib.pyplot as plot

f 			= 1000000 # 1MHz. 
phi 		= 0 # phase shift
A   		= 0.5
fs 			= 1e+09 # sampling frequency.
duration 	= 0.00001
time = np.arange(0, duration, 1/fs); # time array
#amp = 1/2(sin(x)+|sin(x)|)
amplitude=0.5*(abs(np.sin(2*np.pi*np.arange(fs*duration)*f/fs+phi))+np.sin(2*np.pi*np.arange(fs*duration)*f/fs+phi))

#
# This is what we cut and paste into the equation area for the SIM4Life simulation
# cutandpaste = 0.5*sin(2*pi*0.75E9*_t) + 0.5*abs(sin(2*pi*0.75E9*_t))
# 

print(len(amplitude))
# Plot a sine wave using time and amplitude obtained for the sine wave
plot.plot(time, amplitude)
# Give a title for the sine wave plot
plot.title('Unipolar Pulse')
# Give x axis label for the sine wave plot
plot.xlabel('Time')
plot.ylabel('Amplitude')
#plot.grid()
plot.grid(True, which='both')
plot.axhline(y=0, color='k')
plot.show()

