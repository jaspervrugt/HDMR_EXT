# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      777777   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE              77   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE          77    #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE            77     #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE       77      #
#                                                                         #
# ####################################################################### #

import sys
import os
import numpy as np

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

# Parameters
d = 10    # Number of parameters
N = 5000  # Number of samples

# Generate samples from U(0, 1) for X
X = np.random.rand(N, d)

# a-values of Sobol's g-function (you can modify this as needed)
a = 100 * np.random.rand(1, d)

# Sobol g-function computation

# Compute y-values based on Sobol's g-function formula
X_expanded = np.tile(X, (1, 1))  # Replicate X along the axis
a_expanded = np.tile(a, (N, 1))  # Replicate 'a' along the samples axis

y = np.sum((np.abs(4*X - 2) + a_expanded) / (a_expanded + 1), axis=1)

# Analytical variance-based estimates (D_i and D)
D_i = 1 / (3 * (1 + a)**2)  # Sensitivity estimates
D = np.sum(D_i)             # Sum over all parameters

# Sensitivity of each parameter based on D_i
ST_an = 1 / D * D_i.T        # Sensitivity analysis results

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

if __name__ == '__main__':
	# Run the HDMR_EXT toolbox
	S, Ss, Fx, Em, Xy = HDMR_EXT(X, y, options)

	# Numerical & Analytic sensitivity estimates matching
	S_values = S[1:d+1, 6]  # Extract the relevant values from S (column 7, 1-indexed to 6)
	ST_an_flat = ST_an.flatten()  # Flatten ST_an to match dimensions

	# Calculate the error between numerical and analytic estimates
	err = np.sum(np.abs(S_values - ST_an_flat)) / d
	print(f"Error between numerical and analytic sensitivity estimates: {err}")
