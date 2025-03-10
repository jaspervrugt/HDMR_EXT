# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   11   0000   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       11  00  00  #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE    11  00  00  #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE       11  00  00  #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   11   0000   #
#                                                                         #
# ####################################################################### #
#
# Example 2 of the following paper
#  Chastaing, G., F. Gamboa, and C. Prieur (2012), Generalized Hoeffding-
#      Sobol decomposition for dependent variables - application to 
#      sensitivity analysis, Electron. J. Statistics, 6, pp. 2420-2448, 
#      https://doi.org/10.1214/12-EJS749

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
from user_draw2 import user_draw2

d = 4               # Number of parameters
N = 1000            # Number of samples
mu1 = np.zeros(d)   # Sample means for first normal (mu1)
mu2 = np.zeros(d)   # Sample means for second normal (mu2)
w1 = 0.2            # Weight of the first normal
r12_1 = 0.4         # Covariance X1, X2 for 1st normal
r12_2 = 0.37        # Covariance X1, X2 for 2nd normal

# Covariance matrices for each normal distribution
C1 = np.array([[0.5, r12_1], [r12_1, 0.5]])
C2 = np.array([[0.7, r12_2], [r12_2, 0.3]])

# Draw N samples from the normal mixture using the user_draw2 function
X = user_draw2(mu1[:2], mu2[:2], C1, C2, w1, N)

# Compute the function y = f(X)
y = 5*X[:, 0] + 4*X[:, 1] + 3*X[:, 2] + 2*X[:, 3]

# Normalize the X-values
X_min = np.min(X, axis=0)
X_max = np.max(X, axis=0)
X = (X - X_min) / (X_max - X_min)

# Add random noise with variance scaled by sigma2
sigma2 = np.var(y) / 100
y += np.random.normal(0, np.sqrt(sigma2), N)

# Define HDMR_EXT options
options = {
    'graphics': 1,    # Graphics flag [screen output=1]
    'basis': 3,       # choice of orthonormal basis [orthonormal=1, legendre=2, hermite=3]
    'maxorder': 2,    # maximum order for HDMR_EXT expansion
    'm': 3,           # Polynomial degree
    'K': 50,          # Number of bootstrap trials
    'R': N/2,         # Number of samples for each bootstrap trial
    'alfa': 0.01,     # Significance level for F-test [component function significance]
    'method': 1,      # Selection method [forward=1, backward=2]
}

if __name__ == '__main__':
	# Run the HDMR_EXT toolbox
	S, Ss, Fx, Em, Xy = HDMR_EXT(X, y, options)

	# Print the results (you can modify this to fit your use case)
	print("Sensitivity indices (S):", S)
	print("Sensitivity indices (Ss):", Ss)
	print("Function values (Fx):", Fx)
	print("Errors (Em):", Em)
	print("Model outputs (Xy):", Xy)
