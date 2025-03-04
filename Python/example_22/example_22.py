# ####################################################################### #
#                                                                         #
# EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   222222 222222 #
# EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       22 22  22 22  #
# EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE      22     22   #
# EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        22     22    #
# EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   222222 222222 #
#                                                                         #
# ####################################################################### #

# Example of the following paper
#  Gao, Y., A. Sahin, and J.A. Vrugt (2023), Probabilistic Sensitivity 
#      Analysis With Dependent Variables: Covariance-Based Decomposition 
#      of Hydrologic Models, Water Resources Research, 59 (4), 
#      e2022WR0328346, https://doi.org/10.1029/2022WR032834

import sys
import os

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

import numpy as np
from scipy.stats import qmc
from concurrent.futures import ProcessPoolExecutor
from HDMR_EXT import HDMR_EXT
from HDMR_EXT_functions import LH_sampling
from rainfall_runoff import rainfall_runoff

# Dimensionality of the model
d = 7
N = 1000  # Number of samples used by HDMR_EXT

# Parname:      Imax Smax Qsmax  alE alF Kfast Kslow
Par_info = {
    'min': np.array([0.5, 10, 0, 1e-6, -10, 0, 0]),  # minimum values
    'max': np.array([10, 1000, 100, 100, 10, 10, 150])  # maximum values
}
label_par = ['$I_{\rm max}$', '$S_{\rm max}$', '$Q_{\rm s,max}$', 
             '$\alpha_{\rm E}$', '$\alpha_{\rm F}$', '$K_{\rm f}$', '$K_{\rm s}$']


## FIX
# Dummy function for rainfall_runoff (replace with actual model)
def rainfall_runoff(x):
    # Placeholder for the actual model function
    return np.sum(x)  # Example, simply summing the parameters

# Generate Latin Hypercube samples for the parameters
X = LH_sampling(Par_info['min'], Par_info['max'], N)

# Initialize Y matrix (based on rainfall_runoff output)
Y = np.zeros((1, N))  # Assuming one output per run for simplicity

# Run the model once to determine n (the output size)
Y[:, 0] = rainfall_runoff(X[0, :d])  # Initial value

# Parallelize the remaining model evaluations
with ProcessPoolExecutor() as executor:
    results = list(executor.map(lambda x: rainfall_runoff(x), X[1:, :d]))

# Fill the Y matrix with results
Y[:, 1:] = np.array(results).T  # Transpose to match the shape

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

# Initialize SA as a 3D numpy array (9x15xN), where N is defined elsewhere
ny = Y.shape[0]  # Get the number of rows
N = Y.shape[1]   # Get the number of columns

results = np.empty((9, 15, ny))  # Initialize the results array (adjust size as needed)

if __name__ == '__main__':
    	# Call the HDMR_EXT toolbox
	for t in range(n):
	    	results[:,:,t], Ss, Fx, Em, XY = HDMR(X, Y[t, :N], options)
