# ####################################################################### #
#                                                                         #
#  EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   222222 88888 #
#  EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       22 22  88 88 #
#  EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE      22   88888 #
#  EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        22    88 88 #
#  EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   222222 88888 #
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

# Set parameters
d = 2  # Number of parameters
N = 5000  # Number of samples to be used

# Generate random samples from normal distribution (mean=0, std=1)
X = np.random.normal(0, 1, (N, 1))

# Create second column as square of first column
X = np.hstack([X, X[:, 0:1]**2])

# Define the true function
y = X[:, 0] + X[:, 1]

# Variance of random error
sigma2 = np.var(y) / 100

# Add random error
y = y + np.random.normal(0, np.sqrt(sigma2), N)

# Define options for HDMR_EXT
options = {
    'graphics': 1,    # Graphics flag [screen output=1]
    'basis': 1,       # choice of orthonormal basis [orthonormal=1, legendre=2, hermite=3]
    'maxorder': 2,    # maximum order for HDMR_EXT expansion
    'm': 3,           # Polynomial degree
    'K': 1,           # Number of bootstrap trials
    'R': N//5,        # Number of samples for each bootstrap trial
    'alfa': 0.01,     # Significance level for F-test [component function significance]
    'method': 1,      # Selection method [forward=1, backward=2]
}


if __name__ == '__main__':
	# Run the HDMR_EXT toolbox
	S, Ss, Fx, Em, Xy = HDMR_EXT(X, y, options)

	# Output the results (optional)
	print(S, Ss, Fx, Em, Xy)
