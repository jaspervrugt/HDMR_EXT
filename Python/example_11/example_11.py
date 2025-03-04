# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE    11  11     #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE        11  11     #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE     11  11     #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        11  11     #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE    11  11     #
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
from test_function import test_function

# Modified Sobol test function

#    Reference:
#    Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto,
#    M., and Tarantola, S. (2010). Variance-based sensitivity
#    analysis of model output. Design and estimator for the
#    total sensitivity index. Computer Physics Communications,
#    181, 259-270.

#    Configuration Used:
#    Campolongo, F., Saltelli, A., and Cariboni, J. (2011).
#    From screening to quantitative sensitivity analysis. A
#    unified approach. Computer Physics Communications,
#    182, 978-988.

d = 20  	# Number of parameters
N = 5000  	# Number of samples

# Sample X-values from U(0, 1)
X = np.random.rand(N, d)

# Initialize y values
y = np.empty(N)

# Loop through each sample to evaluate the test function
for i in range(N):
    y[i] = test_function(X[i])

# HDMR_EXT options
options = {
    'graphics': 1,    # Graphics flag [screen output=1]
    'basis': 3,       # choice of orthonormal basis [orthonormal=1, legendre=2, hermite=3]
    'maxorder': 2,    # maximum order for HDMR_EXT expansion
    'm': 3,           # Polynomial degree
    'K': 10,          # Number of bootstrap trials
    'R': N/5,         # Number of samples for each bootstrap trial
    'alfa': 0.01,     # Significance level for F-test [component function significance]
    'method': 1,      # Selection method [forward=1, backward=2]
}

if __name__ == '__main__':
	# Run the HDMR_EXT toolbox
	S, Ss, Fx, Em, Xy = HDMR_EXT(X, y, options)

	# Print the results
	print(S, Ss, Fx, Em, Xy)

