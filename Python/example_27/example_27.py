# ####################################################################### #
#                                                                         #
#  EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   222222 77777 #
#  EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       22 22     77 #
#  EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE      22     77  #
#  EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        22     77   #
#  EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   222222 77    #
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

# Define the number of samples
N = int(1e5)

# Sample from U(0,1) N times (uniform distribution)
X = np.random.rand(N, 3)

# Define the function f(x)
def f(x):
    return x[:, 0] + 2 * x[:, 1] + 3 * x[:, 2] + x[:, 0] * x[:, 1]

# Compute the function output y
y = f(X)

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

	# Output the results (optional)
	print(S, Ss, Fx, Em, Xy)
