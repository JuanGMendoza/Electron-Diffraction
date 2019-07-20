import numpy as np
import matplotlib.pyplot as plt

RADIUS  = .0646 #Meters
C = 3.0 * np.power(10,8)
ENERGY_ELECTRON = .511 * np.power(10,6) #eV
MASS_ELECTRON = 9.109 *(10 ** -31)
PLANCK = 6.63 * (10 ** -34)

def N1():
    
    voltages = (5, 4.5, 4, 3.5, 3, 2.5)
    arcLength = [(0),(0),(0),(0),(0),(0)] #all in mm
    alpha = [0,0,0,0,0,0]
    theta = [0,0,0,0,0,0]
    momentum = [0,0,0,0,0,0]
    sinTheta = [0,0,0,0,0,0]
    xValues = [0,0,0,0,0,0]
    errorY = [0,0,0,0,0,0]
    
    arcLength[0] = (10.04, 9.50, 9.02, 9.53) #5Kv

    arcLength[1] = (10.57, 10.37, 10.41, 10.86) #4.5Kv

    arcLength[2] = (11.56, 11.56, 11.49, 11.49) #4Kv

    arcLength[3] = (12.77, 12.10, 12.57, 12.40) #3.5 Kv

    arcLength[4] = (13.17, 12.84, 13.30, 13.37) #3 Kv

    arcLength[5] = (14.84, 13.85, 14.61, 14.45) #2.5Kv

    
   
    
    
    for i in range (0,6):
        
        arcLength[i] = (np.mean(arcLength[i])) * .001 #Convert to Meters
        
        alpha[i] = arcLength[i]/RADIUS
        
        theta[i] = np.arctan(np.sin(alpha[i])/(1+np.cos(alpha[i])))

        sinTheta[i] = np.sin(theta[i])

        momentum[i] = Voltage_to_Momentum(voltages[i])
        xValues[i] = PLANCK / (2*momentum[i])
        errorY[i] = error_prop(arcLength[i],theta[i])


    plt.errorbar(xValues, sinTheta, yerr = errorY, fmt = 'o', label = 'Data')

    slope, cov = np.polyfit(xValues, sinTheta,1, cov = True)
    print(cov)
    newSqueleton = np.poly1d(slope)
    newFunc = newSqueleton(xValues)
    plt.xlabel('h/2p')
    plt.ylabel('sin(θ)')
    plt.title('Slope = ' + str(round(slope[0] * 10,1)) + '      d = ' + str(round(1/(slope[0] * 10),12)))
    plt.plot(xValues,newFunc, label = 'Linear Fit')
    plt.legend()
    plt.show()

def N2():
    voltages = (5, 4.5, 4, 3.5, 3, 2.5)
    arcLength = [(0),(0),(0),(0),(0),(0)] #all in mm
    alpha = [0,0,0,0,0,0]
    theta = [0,0,0,0,0,0]
    momentum = [0,0,0,0,0,0]
    sinTheta = [0,0,0,0,0,0]
    xValues = [0,0,0,0,0,0]
    errorY = [0,0,0,0,0,0]
    
    arcLength[0] = (18.33, 18.33, 18.95, 16.9) #5Kv

    arcLength[1] = (18.69, 18.60, 19.05, 19.04) #4.5Kv

    arcLength[2] = (20.52, 19.42, 20.51, 19.85) #4Kv

    arcLength[3] = (21.81, 19.72, 21.3, 22.16) #3.5 Kv

    arcLength[4] = (23.53, 22.64, 23.29, 22.32) #3 Kv

    arcLength[5] = (25.28, 23.87, 25.7, 24.88) #2.5Kv

   
    
    
    for i in range (0,6):
        
        arcLength[i] = (np.mean(arcLength[i])) * .001 #Convert to Meters
        
        alpha[i] = arcLength[i]/RADIUS
        
        theta[i] = np.arctan(np.sin(alpha[i])/(1+np.cos(alpha[i])))

        sinTheta[i] = np.sin(theta[i])

        errorY[i] = error_prop(arcLength[i],theta[i])
        
        momentum[i] = Voltage_to_Momentum(voltages[i])
        xValues[i] = PLANCK / (momentum[i])

        

    
    plt.errorbar(xValues, sinTheta, yerr = errorY, fmt = 'o', label = 'Data')

    slope = np.polyfit(xValues, sinTheta,1)
    newSqueleton = np.poly1d(slope)
    newFunc = newSqueleton(xValues)
    plt.xlabel('h/p')
    plt.ylabel('sin(θ)')
    plt.title('Slope = ' + str(round(slope[0] * 10,1)) + '      d = ' + str(round(1/(slope[0] * 10),12)))
    plt.plot(xValues,newFunc, label = 'Linear Fit')
    plt.legend()
    plt.show()
    
def Voltage_to_Momentum(Voltage):
    K_Energy = Voltage #eV

    fraction =   ENERGY_ELECTRON / (K_Energy + ENERGY_ELECTRON)
    velocity = np.sqrt(1.0 - np.power(fraction,2)) * C
    gamma = 1.0 / (1.0 - np.power((velocity/C),2))
    return  velocity * MASS_ELECTRON * gamma

def error_prop(s, theta):

    dS = 1.49 #mm
    dR = .5 #mm
    r = 64.6 #mm
    
    dS_R = (s/r) * np.sqrt(np.power(dR/r,2) + np.power(dS/s,2))
    dSin = np.cos(s/r) * dS_R
    dCos = -np.sin(s/r) * dS_R
    dTheta = np.abs(np.sin(s/r)/(1+np.cos(s/r))) * np.sqrt(np.power(dSin/np.sin(s/r),2) + np.power(dCos/np.cos(s/r),2))
    dArctan = (1/(1 + (np.sin(s/r)/(1 + np.cos(s/r))))) * dTheta
    dSinFinal = np.cos(theta) * dArctan
    return dSinFinal
N2()

#round(4.33, 1)





