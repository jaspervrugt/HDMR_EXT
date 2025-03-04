# ####################################################################### #
#                                                                         #
# EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   222222 44     #
# EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       22 22  44  44 #
# EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE      22   44  44 #
# EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        22    444444 #
# EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   222222     44 #
#                                                                         #
# ####################################################################### #

# Example of the following paper
#  Cariboni, J., D. Gatelli, R. Liska, and A. Saltelli (2007), The role of 
#      sensitivity analysis in ecological modelling, Ecological Modelling, 
#      203, pp. 167–182.
#      http://izt.ciens.ucv.ve/ecologia/Archivos/ECO_POB%202007/ ...
#          ECOPO2_2007/Cariboni%20et%20al%202007.pdf
#  Massoud, E.C., J. Huisman, E. Benincà, M.C. Dietze, W. Bouten, and J.A.
#      Vrugt (2018), Probing the limits of predictability: data 
#      assimilation of chaotic dynamics in complex food webs, Ecology
#      Letters, 21 (1), pp 93-103, https://doi.org/10.1111/ele.12876

import sys
import os
import numpy as np
import scipy.io as sio
from scipy.integrate import odeint

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
from evaluate_PP import evaluate_PP
from HDMR_EXT_functions import LH_sampling

d = 4  		# number of parameters
N = 100  	# number of simulations
label_par = ['$r$', '$\\alpha$', '$m$', '$\\theta$']

X_min = np.array([0.8, 0.2, 0.6, 0.05])  # min values of parameters
X_max = np.array([1.8, 1.0, 1.0, 0.15])  # max values of parameters

# Generate N parameter vectors
X = LH_sampling(X_min, X_max, N)

# Initial conditions for prey and predator populations
u0 = [8, 2]

# ODE settings
ode_options = {'atol': 1e-5, 'rtol': 1e-5}  	# Absolute and relative tolerance
dt = 1  					# Time step
t_max = 60  					# Maximum simulation time
n = int(t_max / dt) + 1  			# Number of time points

# Check if precomputed simulation data exists
if os.path.exists('XY.mat'):
    # Load pre-simulated data
    data = sio.loadmat('XY.mat')
    X = data['X']
    Y = data['Y']
else:
    # Simulate the predator-prey model if data doesn't exist
    Y = evaluate_PP(X, n, d, u0, ode_options, dt, t_max)
    # Save X and Y to a .mat file for future use
    sio.savemat('XY.mat', {'X': X, 'Y': Y})

# HDMR_EXT options
options = {
    'graphics': 1,    # Graphics flag [screen output=1]
    'basis': 1,       # choice of orthonormal basis [orthonormal=1, legendre=2, hermite=3]
    'maxorder': 2,    # maximum order for HDMR_EXT expansion
    'm': 3,           # Polynomial degree
    'K': 1,           # Number of bootstrap trials
    'R': N,           # Number of samples for each bootstrap trial
    'alfa': 0.01,     # Significance level for F-test [component function significance]
    'method': 1,      # Selection method [forward=1, backward=2]
}

# Run HDMR_EXT toolbox for each of the K = 2 species
results = {}

if __name__ == '__main__':
    	# Call the HDR_EXT toolbox
	for k in range(2):  # Loop over each species
    		for t in range(10):  # Loop over the first 10 time points
        		Ss, Fx, Em, XY = HDMR_EXT(X, Y[:, t, k], options)
        		results[(k, t)] = {
            			'Ss': Ss,
		            	'Fx': Fx,
			        'Em': Em,
			        'XY': XY}

# You can now analyze and plot the results stored in `results`
# If you look at t = 9 and t = 61 then you get table 1 as in referenced 
# paper above. Results match quite closely. So for a paper we can start 
# with this case study and then followed by the two-predator-two-prey model
