# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE        1111   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE           11 11   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       11  11   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              11   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE          11   #
#                                                                         #
# ####################################################################### #

# Example from the following paper
#  Li, G., and H. Rabitz (2012), General formulation of HDMR component 
#      functions with independent and correlated variables, J. 
#      Math. Chem., 50, pp. 99-130

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

d = 5                      # Number of dimensions
N = 5000                   # Number of samples
mu = 0.5 * np.ones(d)      # Sample mean, µ
C = np.eye(d)              # Covariance matrix, Σ (Identity matrix)
  
# Draw N samples from N(µ, Σ) (Multivariate Normal Distribution)
X = np.random.multivariate_normal(mu, C, N)

# y = f(x) function: Sum of all columns in X (row-wise sum)
y = np.sum(X[:, :d], axis=1)

# Variance of random error
sigma2 = np.var(y) / 100

# Add random error (Noise)
y = y + np.random.normal(0, np.sqrt(sigma2), N)

# HDMR_EXT options 
options = {
    'graphics': 1,
    'basis': 1,
    'maxorder': 2,
    'm': 3,
    'K': 1,
    'R': 1000,
    'alfa': 0.01,
    'method': 1,
    'tolopt': 1e-3
}

if __name__ == '__main__':
	# Run the HDMR toolbox
	S, Ss, Fx, Em, Xy = HDMR_EXT(X, y, options)
