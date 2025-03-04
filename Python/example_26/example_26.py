# ####################################################################### #
#                                                                         #
#  EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   222222 66    #
#  EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       22 22  66    #
#  EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE      22   66666 #
#  EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        22    66 66 #
#  EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   222222 66666 #
#                                                                         #
# ####################################################################### #

# Please check the following papers
#  Massoud, E.C., J. Huisman, E. Benincà, M.C. Dietze, W. Bouten, and J.A.
#      Vrugt (2018), Probing the limits of predictability: data 
#      assimilation of chaotic dynamics in complex food webs, Ecology
#      Letters, 21 (1), pp 93-103, https://doi.org/10.1111/ele.12876
#  Vandermeer, J. (1993), Loose coupling of predator-prey cycles: 
#      entrainment, chaos, and intermittency in the classic MacArthur 
#      consumer-resource equations, American Naturalist, 141, pp. 687-716
#  Vandermeer, J. (2004), Coupled oscillations in food webs: balancing 
#      competition and mutualism in simple ecological models. American 
#      Naturalist, 163, pp. 857-867 
#  Vandermeer, J. (2006), Oscillating populations and biodiversity 
#      maintenance, Bioscience, 56, pp. 967-975

import numpy as np
import scipy.io as sio
from scipy.integrate import odeint
import os
import random
from multiprocessing import Pool
from VDM_model import VDM_model

# Get the current working directory
current_directory = os.getcwd()
# Go up one directory
parent_directory = os.path.abspath(os.path.join(current_directory, '..'))
# add this to path
sys.path.append(parent_directory)
# Add another directory
misc_directory = os.path.abspath(os.path.join(parent_directory, 'miscellaneous'))
# add this to path
sys.path.append(misc_directory)

from HDMR_EXT import HDMR_EXT
from HDMR_EXT_functions LH_sampling

model_version = 1  	# Model version (1 or 2)
N = 10 			# Number of simulations

if model_version == 1:
    d = 8
    X_min = [0.0001, 0.0001, 0.01, 0.01, 0.1, 0.50, 0.001, 1]
    X_max = [1, 2, 2.50, 2.50, 2.5, 2.50, 0.700, 3]
    label_par = ['$\\beta$', '$\\alpha$', '$r_{1}$', '$r_{2}$', '$g$', '$K$', '$m$', '$H$']
elif model_version == 2:
    d = 11
    X_min = [0.0001, 0.0001, 0.0001, 0.01, 0.01, 0.1, 0.50, 0.50, 0.010, 0.010, 1]
    X_max = [1, 2, 2, 2.50, 2.50, 2.5, 2.50, 2.50, 0.700, 0.700, 3]
    label_par = ['$\\beta$', '$\\alpha_{1}$', '$\\alpha_{2}$', '$r_{1}$', '$r_{2}$', '$g$', '$K_{1}$', '$K_{2}$', '$m_{1}$', '$m_{2}$', '$H$']

# Generate N parameter vectors using Latin Hypercube Sampling
X = LH_sampling(np.array(X_min), np.array(X_max), N)

# Initial conditions (prey and predator populations)
u0 = [1.0096, 1.6658, 0.5705, 1.0006] if model_version == 1 else [1.0096, 1.6658, 0.5705, 1.0006, 1.0096, 1.6658, 0.5705, 1.0006]

# ODE options
ode_options = {'atol': 1e-5, 'rtol': 1e-5}  	# Absolute and relative tolerance
dt = 3.35
t_max = 2656.55
n = int(t_max / dt + 1)  			# Number of time steps
K = 4  						# Number of species

# Define the model version
model_version = 1  # This needs to be assigned based on your version

# Check whether the file with model simulations exists
file_name = f'XY_{model_version}.mat'

if not os.path.isfile(file_name):  # File does not exist
    # Simulate the model if data doesn't exist
    X, Y, SS = VDM_model(X, n, d, u0, ode_options, dt, t_max, model_version)
    # Save X, Y, and SS to a .mat file for future use
    # Now compute the residual with observed data
    # Read the Excel data
    xls_data = pd.read_excel('data_Beninca_et_al_Nature&Ecology_Letters.xlsx', sheet_name='sheet1', usecols="A:E", nrows=794)
    MeasData = xls_data.iloc[:, 1:5].values.T  # Measured data (4 species)
    # Initialize matrix of residuals
    Res = np.full((N, n, K), np.nan)   
    # Loop to determine residuals
    for k in range(K):  # Loop over species
        for i in range(N):  # Loop over parameter vectors
            Res[i, :, k] = Y[i, :, k] - MeasData[k, :]  # Residual calculation
    
    # Save X, Y, SS, and Res with model version
    savemat(f'XY_{model_version}.mat', {'X': X, 'Y': Y, 'SS': SS, 'Res': Res})

else:
    # Load the existing data
    data = sio.loadmat(f'XY_{model_version}.mat')
    X = data['X']
    Y = data['Y']
    SS = data['SS']
    Res = data['Res']

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Now run HDMR_EXT code for each simulation time
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# Define options for HDMR_EXT
options = {
    'graphics': 0,    # Graphics flag [screen output=1]
    'basis': 1,       # choice of orthonormal basis [orthonormal=1, legendre=2, hermite=3]
    'maxorder': 2,    # maximum order for HDMR_EXT expansion
    'm': 3,           # Polynomial degree
    'K': 1,           # Number of bootstrap trials
    'R': N,           # Number of samples for each bootstrap trial
    'alfa': 0.01,     # Significance level for F-test [component function significance]
    'method': 1,      # Selection method [forward=1, backward=2]
}

# Get the shape of Res (n, K)
_, n, K = Res.shape

# Loop over each species
if __name__ == '__main__':
	for k in range(K):
    		for t in range(1, 10):  # Limited to 10 time points
		    	# Call the HDMR_EXT toolbox
		        results, Ss, Fx, Em, XY = HDMR_EXT(X, Res[:, t, k], options)  # Need to downselect Y!
        
			# Store the results, if needed (e.g., saving or further processing)
        		# results[:, :, t, k] = results
		        # Further processing can be done here

