# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      888888   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          88  88   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       888888   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE          88  88   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      888888   #
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
from corr2cov import corr2cov

d = 3    # Number of parameters
N = 5000  # Number of samples

# Mean vector (Sample means µ)
mu = 0.5 * np.ones(d)  # Creates a vector of 0.5s of length d

# Correlation matrix
R = np.array([[1, 0.5, 0.4], 
              [0.5, 1, 0.7], 
              [0.4, 0.7, 1]])

# Sample standard deviations (assuming 1 for each)
std_devs = np.ones(d)

# Covariance matrix
C = corr2cov(std_devs, R)

# Draw N samples from N(mu, C)
X = np.random.multivariate_normal(mu, C, N)

# Normalize X (min-max scaling)
X = (X - np.min(X, axis=0)) / (np.max(X, axis=0) - np.min(X, axis=0))

# Define the function y = f(X)
y = 2 * X[:, 0] + X[:, 1] + 3 * np.exp(X[:, 0] * X[:, 1] * X[:, 2])

# Add random error
sigma2 = np.var(y) / 100  # Variance of random error
y += np.random.normal(0, np.sqrt(sigma2), N)

# HDMR_EXT options as a dictionary
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

	# Output for debugging or verification (displaying the first few rows of the result)
	print("Sensitivity results (S):\n", S)
