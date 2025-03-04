# ####################################################################### #
#                                                                         #
#  HHH   HHH DDDDDDDD  MMM    MMM RRRRRRR    EEEEEEEE XXX   XXX TTTTTTTTT #
#  HHH   HHH DDDDDDDDD MMM    MMM RRRRRRRR   EEEEEEEE XXX   XXX TTTTTTTTT #
#  HHH   HHH DDD   DDD MMM    MMM RRR   RRR  EEE       XXX XXX  TT TTT TT #
#  HHH   HHH DDD   DDD MMMM  MMMM RRR   RRR  EEE       XXX XXX  T  TTT  T #
#  HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRRR   EEEEEE     XXXXX      TTT    #
#  HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRR    EEEEEE     XXXXX      TTT    #
#  HHH   HHH DDD   DDD MMM    MMM RRRRRRR    EEE       XXX XXX     TTT    #
#  HHH   HHH DDD   DDD MMM    MMM RRR  RRR   EEE       XXX XXX     TTT    #
#  HHH   HHH DDDDDDDDD MMM    MMM RRR   RRR  EEEEEEEE XXX   XXX    TTT    #
#  HHH   HHH DDDDDDDD  MMM    MMM RRR   RRR  EEEEEEEE XXX   XXX    TTT    #
#                                                                         #
# ####################################################################### #
#                                                                         #
# Extended Bases High-Dimensional Model Representation (HDMR) uses three  #
# different type of orthonormal polynomials (Classical, Legendre, and     #
# Hermite) for variance-based global sensitivity analysis (GSA) with      #
# correlated and uncorrelated inputs. This function uses N x d matrix of  #
# N different vector d-vectors of model inputs (factors/parameters) and a #
# N x 1 vector of corresponding model outputs and returns to the user     #
# each factor's first, second, and third order sensitivity coefficient    #
# (separated in total, structural and correlative contributions), an      #
# estimate of their 95# confidence intervals (from bootstrap method)      #
# and the coefficients of the extended bases functions that govern output,#
# Y emulator (determined by an F-test of the error residuals of  the HDMR #
# model with/without a given first, second and/or third order components).#
# These coefficients define an emulator that can be used to predict the   #
# output, Y, of the original model for any d-vector of model inputs. For  #
# uncorrelated model inputs (columns of X are independent), the HDMR      #
# sensitivity indices reduce to a single index (= structural contribution)#
# consistent with their values derived from commonly used variance-based  #
# GSA methods.                                                            #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  MAIN REFERENCE                                                         #
#   Li, G. H. Rabitz, P.E. Yelvington, O.O. Oluwole, F. Bacon,            #
#       C.E. Kolb, and J. Schoendorf (2010), Global sensitivity analysis  #
#       for systems with independent and/or correlated inputs, Journal of #
#       Physical Chemistry A, Vol. 114 (19), pp. 6022 - 6032, 2010        #
#   Gao, Y., A. Sahin, and J.A. Vrugt (2023), Probabilistic Sensitivity   #
#       Analysis With Dependent Variables: Covariance-Based Decomposition #
#       of Hydrologic Models, Water Resources Research, 59 (4),           #
#       e2022WR0328346, https://doi.org/10.1029/2022WR032834              #
#                                                                         #
#  PYTHON CODE                                                            #
#  © Written by Jasper A. Vrugt using GPT-4 OpenAI's language model       #
#    University of California Irvine                                      #
#  Version 2.0    Dec 2024                                                #
#                                                                         #
# ####################################################################### #

import numpy as np
import sys
import os
import itertools
from scipy.stats import f
from scipy.special import erfcinv
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from pandas.plotting import table
from matplotlib.backends.backend_pdf import PdfPages
from screeninfo import get_monitors

def HDMR_EXT_setup(X, y, user_options):
    """
    Setup function for HDMR_EXT, checks and sets up input arguments.

    Parameters:
    - X (numpy.ndarray): Nxd matrix of N vectors of d parameters
    - y (numpy.ndarray): Nx1 vector with one model output for each row of X
    - user_options (dict): Dictionary containing the options for HDMR variables

    Returns:
    - N (int): Number of rows in X
    - d (int): Number of columns in X
    - graphics (int): Option for graphics (0 or 1)
    - basis (int): Basis function option (1, 2, or 3)
    - maxorder (int): Maximum order of the interaction (1, 2, or 3)
    - m (int): Parameter for the number of samples (1 to 12)
    - K (int): Parameter K (1 to 500)
    - R (int): Number of rows for R (between 100 and N)
    - method (int): Method for forward/backward selection (1 or 2)
    - alfa (float): Significance level (0.01 to 0.1 for 99 to 90% intervals)
    """

    # Determine the size of matrix X
    N, d = X.shape

    if N < 300:
        raise ValueError('HDMR_EXT ERROR: Number of samples, N, of N x d matrix X is insufficient')

    if y.shape[0] != N:
        raise ValueError('HDMR_EXT ERROR: Dimension mismatch. The number of rows, N of Y should match number of rows of X')

#    if y.shape[1] != 1:
#        raise ValueError('HDMR_EXT ERROR: Y should be a N x 1 vector with one simulated output for each N parameter vectors')

    # Define default options
    def_options = {
        'graphics': 1,
        'basis': 1,
        'maxorder': 2,
        'm': 3,
        'K': 1,
        'R': N // 2,
        'alfa': 0.01,
        'method': 1,
    }

    # Open file for settings.txt
    with open('HDMR_EXT_settings.txt', 'w') as fid:
        fid.write('|-------------------------------------------------------------------------|\n')
        fid.write('|                                                                         |\n')
        fid.write('| HHH   HHH DDDDDDDD  MMM    MMM RRRRRRR     EEEEEEEE XXX   XXX TTTTTTTTT |\n')
        fid.write('| HHH   HHH DDDDDDDDD MMM    MMM RRRRRRRR    EEEEEEEE XXX   XXX TTTTTTTTT |\n')
        fid.write('| HHH   HHH DDD   DDD MMM    MMM RRR   RRR   EEE       XXX XXX  TT TTT TT |\n')
        fid.write('| HHH   HHH DDD   DDD MMMM  MMMM RRR   RRR   EEE       XXX XXX  T  TTT  T |\n')
        fid.write('| HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRRR    EEEEEE     XXXXX      TTT    |\n')
        fid.write('| HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRR     EEEEEE     XXXXX      TTT    |\n')
        fid.write('| HHH   HHH DDD   DDD MMM    MMM RRRRRRR     EEE       XXX XXX     TTT    |\n')
        fid.write('| HHH   HHH DDD   DDD MMM    MMM RRR  RRR    EEE       XXX XXX     TTT    |\n')
        fid.write('| HHH   HHH DDDDDDDDD MMM    MMM RRR   RRR   EEEEEEEE XXX   XXX    TTT    |\n')
        fid.write('| HHH   HHH DDDDDDDD  MMM    MMM RRR   RRR   EEEEEEEE XXX   XXX    TTT    |\n')
        fid.write('|                                                                         |\n')
        fid.write('|  SSSSSSS EEEEEEEE TTTTTTTT TTTTTTTT IIIIIIII NN    NN  GGGGGG   SSSSSSS |\n')
        fid.write('| SSSSSSS  EEEEEEEE TTTTTTTT TTTTTTTT  IIIIII  NNN   NN GGGGGG   SSSSSSS  |\n')
        fid.write('| SS       EE       TT TT TT TT TT TT    II    NNNN  NN GG       SS       |\n')
        fid.write('| SSSSSS   EEEEE    T  TT  T T  TT  T    II    NN NN NN GG  GGGG SSSSSS   |\n')
        fid.write('| SSSSSSS  EEEEE       TT       TT       II    NN NN NN GG   GGG SSSSSSS  |\n')
        fid.write('|       SS EE          TT       TT       II    NN  NNNN GG    GG       SS |\n')
        fid.write('|  SSSSSSS EEEEEEEE    TT       TT     IIIIII  NN   NNN GGGGGGGG  SSSSSSS |\n')
        fid.write('| SSSSSSS  EEEEEEEE    TT       TT    IIIIIIII NN    NN  GGGGGG  SSSSSSS  |\n')
        fid.write('|                                                                         |\n')
        fid.write('|-------------------------------------------------------------------------|\n')

        # Unpack structure user_options and apply defaults
        graphics = int(user_options.get('graphics', def_options['graphics']))
        if graphics not in [0, 1]:
            raise ValueError('HDMR_EXT ERROR: Field "graphics" of options should take on the value of 0 or 1')

        maxorder = int(user_options.get('maxorder', def_options['maxorder']))
        if maxorder not in [1, 2, 3]:
            raise ValueError('HDMR_EXT ERROR: Field "maxorder" of options should be an integer with values of 1, 2, or 3')

        basis = int(user_options.get('basis', def_options['basis']))
        if basis not in [1, 2, 3]:
            raise ValueError('HDMR_EXT ERROR: Field "basis" of options should be an integer with values of 1, 2, or 3')

        m = int(user_options.get('m', def_options['m']))
        if not (1 <= m <= 12):
            raise ValueError('HDMR_EXT ERROR: Field "m" of options should be an integer between 1 and 12')

        K = int(user_options.get('K', def_options['K']))
        if not (1 <= K <= 500):
            raise ValueError('HDMR_EXT ERROR: Field "K" of options should be an integer between 1 and 500')

        R = int(user_options.get('R', def_options['R']))
        if not (100 <= R <= N):
            raise ValueError(f'HDMR_EXT ERROR: Field "R" of options should be an integer between 100 and N, where N = {N}')

        method = user_options.get('method', def_options['method'])
        if method not in [1, 2]:
            raise ValueError('HDMR_EXT ERROR: Field "method" of options should be 1 (forward selection) or 2 (backward elimination)')

        alfa = user_options.get('alfa', def_options['alfa'])
        if isinstance(alfa, str):
            raise ValueError('HDMR_EXT ERROR: Field "alfa" should not be a string but a numerical value')
        if not (0 <= alfa <= 0.1):
            raise ValueError('HDMR_EXT ERROR: Field "alfa" should be between 0.0 and 0.1')

        # Save options to file
        fid.write(f'\n          |===================================|\n')
        fid.write(f'          |  field of options      value      |\n')
        fid.write(f'          |-----------------------------------|\n')
        fid.write(f'          |    graphics   \t {graphics:8d}     |\n')
        fid.write(f'          |    basis      \t {basis:8d}     |\n')
        fid.write(f'          |    maxorder   \t {maxorder:8d}     |\n')
        fid.write(f'          |    m          \t {m:8d}     |\n')
        fid.write(f'          |    K          \t {K:8d}     |\n')
        fid.write(f'          |    R          \t {R:8d}     |\n')
        fid.write(f'          |    alfa       \t {alfa:8f}     |\n')
        fid.write(f'          |    method     \t {method:8f}     |\n')
        fid.write(f'          |===================================|\n')

    return N, d, graphics, basis, maxorder, m, K, R, method, alfa


def HDMR_EXT_initialize(X, y, N, d, K, R, m, maxorder, basis):
    """
    Translated Python function for HDMR_EXT_initialize
    
    Parameters:
        X (numpy.ndarray): Nxd matrix of N samples of d parameters.
        y (numpy.ndarray): Nx1 vector of model outputs.
        N (int): Number of input samples.
        d (int): Number of input parameters (dimensions).
        K (int): Number of bootstrap samples.
        R (int): Number of realizations (models).
        m (int): Degree for the polynomial basis functions.
        maxorder (int): Maximum order for HDMR expansion (1, 2, or 3).
        basis (int): Type of basis functions to use (1 = orthonormal, 2 = Legendre, 3 = Hermite).
    
    Returns:
        Xy (dict): Dictionary with normalized inputs and other metadata.
        Em (dict): Dictionary containing HDMR model parameters and coefficients.
        S (dict): Dictionary for storing sensitivity analysis results.
        O (numpy.ndarray): Matrix related to B-splines.
        n_ns (numpy.ndarray): Number of non-singular values.
        ct (numpy.ndarray): Counting terms for higher-order interactions.
    """
    # Set random seed for reproducibility
    np.random.seed(1 + round(100 * np.random.rand()))
    
    # Structure Xy: Define content
    id = np.argsort(np.random.rand(N, K), axis=0)
    Xy = {
        'X_n': np.nan * np.ones((N, d)),
        'minX': np.min(X, axis=0),
        'maxX': np.max(X, axis=0),
        'y': y,
        'R': R,
        'id': id[:R, :K]
        }

    # Basis functions: Select parameters for the chosen basis
    if basis == 1:      ## Orthonormal polynomial
        a, b = 0, 1
    elif basis == 2:    ## Legendre polynomial
        a, b = -1, 1
    elif basis == 3:    ## Hermite polynomial
        a, b = -2, 2
    
    # Normalize the input values
    Xy['X_n'] = a + (b - a) * (X[:N, :d] - Xy['minX']) / (Xy['maxX'] - Xy['minX'])
    
    # Calculate polynomial coefficients (C) based on the basis type
    if basis == 1:      ## Orthonormal polynomials
        C = orthopoly(Xy['X_n'], N, d, m)
    elif basis == 2:    ## Legendre polynomials
        C = legendre(m, d)
    elif basis == 3:    ## Hermite polynomials
        C = hermite(m, d)

    # Combinations for parameter interactions (1st, 2nd, and 3rd order)
    n1 = d  # 1st order terms
    n2 = n3 = beta = gamma = ct = 0
    n_ns = np.nan * np.ones(K)
    
    if maxorder > 1:
        beta = np.array(list(itertools.combinations(range(0, d),2)))  # 2nd order
        n2 = beta.shape[0]

    if maxorder == 3:
        gamma = np.array(list(itertools.combinations(range(0, d),3)))  # 2nd order
        n3 = gamma.shape[0]
    
    n = n1 + n2 + n3
    m1 = m                          # Terms for 1st order
    m2 = 2 * m + m**2               # Terms for 2nd order
    m3 = 3 * m + 3 * m**2 + m**3    # Terms for 3rd order
    
    # Number of terms per order
    k1 = d * m
    k2 = m2 * (d * (d - 1) // 2)
    k3 = m3 * (d * (d - 1) * (d - 2) // 6)
    k = k1 + k2 + k3

    # Initialize structure Em
    if maxorder == 1:
        Em = {
            'nterms': np.nan * np.ones(K),
            'p0': np.nan * np.ones(K),
            'RMSE': np.nan * np.ones(K),
            'm': m,
            'Y_e': np.nan * np.ones((R, K)),
            'f0': np.nan * np.ones(K),
            'm1': m1,
            'n1': d,
            'm2': m2,
            'n2': n2,
            'm3': m3,
            'n3': n3,
            'n': n,
            'maxorder': maxorder,
            'select': np.nan * np.ones((n, K)),
            'k1': k1,
            'k2': k2,
            'k3': k3,
            'k': k1,
            'beta': beta,
            'gamma': gamma,
            'C': np.nan * np.ones((k1, K)),
            'EB': np.zeros((N, k1)),
            'RT': np.zeros(K)
        }
    elif maxorder == 2:
        Em = {
            'nterms': np.nan * np.ones(K),
            'p0': np.nan * np.ones(K),
            'RMSE': np.nan * np.ones(K),
            'm': m,
            'Y_e': np.nan * np.ones((R, K)),
            'f0': np.nan * np.ones(K),
            'm1': m1,
            'n1': n1,
            'm2': m2,
            'n2': n2,
            'm3': m3,
            'n3': n3,
            'n': n,
            'maxorder': maxorder,
            'select': np.nan * np.ones((n, K)),
            'k1': k1,
            'k2': k2,
            'k3': k3,
            'k': k1 + k2,
            'beta': beta,
            'gamma': gamma,
            'C': np.nan * np.ones((k1 + k2, K)),
            'EB': np.zeros((N, k1 + k2)),
            'RT': np.zeros(K)
        }
    elif maxorder == 3:
        Em = {
            'nterms': np.nan * np.ones(K),
            'p0': np.nan * np.ones(K),
            'RMSE': np.nan * np.ones(K),
            'm': m,
            'Y_e': np.nan * np.ones((R, K)),
            'f0': np.nan * np.ones(K),
            'm1': m1,
            'n1': n1,
            'm2': m2,
            'n2': n2,
            'm3': m3,
            'n3': n3,
            'n': n,
            'maxorder': maxorder,
            'select': np.nan * np.ones((n, K)),
            'k1': k1,
            'k2': k2,
            'k3': k3,
            'k': k1 + k2 + k3,
            'beta': beta,
            'gamma': gamma,
            'C': np.nan * np.ones((k, K)),
            'EB': np.zeros((N, k)),
            'RT': np.zeros(K)
        }
  
    # Compute B-splines for all N samples of X_n
    Em['EB'], O = HDMR_EXT_construct_EB(Em['EB'], Xy['X_n'][:N, :d], N, d, m, C, maxorder)

    # Initialize structure S for sensitivity analysis
    S = {
        'S': np.nan * np.ones((Em['n'], K)),
        'Sa': np.nan * np.ones((Em['n'], K)),
        'Sb': np.nan * np.ones((Em['n'], K)),
        'ST': np.nan * np.ones((d, 1)),
        'V_em': np.nan * np.ones((Em['n'], K)),
        'V_y': np.nan * np.ones(K)
    }

    # Print header information
    print('  -----------------------------------------------------------------------            ')
    print('  HHH   HHH DDDDDDDD  MMM    MMM RRRRRRR     EEEEEEEE XXX   XXX TTTTTTTTT            ')
    print('  HHH   HHH DDDDDDDDD MMM    MMM RRRRRRRR    EEEEEEEE XXX   XXX TTTTTTTTT            ')
    print('  HHH   HHH DDD   DDD MMM    MMM RRR   RRR   EEE       XXX XXX  TT TTT TT            ')
    print('  HHH   HHH DDD   DDD MMMM  MMMM RRR   RRR   EEE       XXX XXX  T  TTT  T            ')
    print('  HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRRR    EEEEEE     XXXXX      TTT        /^ ^\  ')
    print('  HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRR     EEEEEE     XXXXX      TTT       / 0 0 \ ')
    print('  HHH   HHH DDD   DDD MMM    MMM RRRRRRR     EEE       XXX XXX     TTT       V\ Y /V ')
    print('  HHH   HHH DDD   DDD MMM    MMM RRR  RRR    EEE       XXX XXX     TTT        / - \  ')
    print('  HHH   HHH DDDDDDDDD MMM    MMM RRR   RRR   EEEEEEEE XXX   XXX    TTT       /     | ')
    print('  HHH   HHH DDDDDDDD  MMM    MMM RRR   RRR   EEEEEEEE XXX   XXX    TTT       V__) || ')
    print('  -----------------------------------------------------------------------            ')
    print('  © Jasper A. Vrugt, University of California Irvine & GPT-4 OpenAI''s language model')
    print('    ________________________________________________________________________')
    print('    Version 2.0, Dec. 2024, Beta-release: MATLAB implementation is benchmark')
    print('\n')

    return Xy, Em, S, n_ns, ct

def orthopoly(X, N, d, m):
    """
    Computes the coefficients of an orthonormal polynomial for a given matrix X.

    Parameters:
    X (ndarray): Input matrix of size (N, d)
    N (int): Number of rows in X
    d (int): Number of columns in X
    m (int): Degree of the polynomial

    Returns:
    C (ndarray): Coefficients of the orthonormal polynomial, shape (m, m+1, d)
    """
    # Initialize moment matrix M
    k = 0
    M = np.zeros((m+1, m+1, d))
    # Compute the moment matrix for each dimension of X
    for i in range(d):
        for j in range(m+1):
            for z in range(m+1):
                M[j, z, i] = np.sum(X[:, i]**k,axis=0) / N
                k += 1
#            k = j
            k = j + 1
        k = 0

    # Initialize coefficient matrix C
    C = np.zeros((m, m+1, d))
    # Calculate the coefficients of the orthonormal polynomial for each dimension of X
    for i in range(d):                 
        for j in range(m):             
            for k in range(j+2):      
                z = list(range(j+2))  
                z.remove(k)
                det_ij = np.linalg.det(M[:j+1, :j+1, i]) * np.linalg.det(M[:j+2, :j+2, i])
                C[j, j+1-k, i] = (-1)**(j+k+3) * np.linalg.det(M[:j+1, z, i]) / np.sqrt(det_ij)

    return C

def hermite(m, d):
    """
    Returns the value of an m-th degree Hermite polynomial at x for a given dimension d.

    Parameters:
    m (int): Degree of the Hermite polynomial
    d (int): Number of dimensions (or instances) to return the polynomial values for

    Returns:
    C (ndarray): A (m, m+1, d) array containing the Hermite polynomials for each dimension
    """
    # Initialize the coefficient matrix C
    C = np.zeros((m, m+1, d))

    # Hermite polynomial coefficients based on the degree m
    if m == 1:
        p = np.array([2, 0])
    elif m == 2:
        p = np.array([[2, 0, 0], [4, 0, -2]])
    elif m == 3:
        p = np.array([[2, 0, 0, 0], [4, 0, -2, 0], [8, 0, -12, 0]])
    elif m == 4:
        p = np.array([[2, 0, 0, 0, 0], [4, 0, -2, 0, 0], [8, 0, -12, 0, 0],
                      [16, 0, -48, 0, 12]])
    elif m == 5:
        p = np.array([[2, 0, 0, 0, 0, 0], [4, 0, -2, 0, 0, 0], [8, 0, -12, 0, 0, 0],
                      [16, 0, -48, 0, 12, 0], [32, 0, -160, 0, 120, 0]])
    elif m == 6:
        p = np.array([[2, 0, 0, 0, 0, 0, 0], [4, 0, -2, 0, 0, 0, 0], [8, 0, -12, 0, 0, 0, 0],
                      [16, 0, -48, 0, 12, 0, 0], [32, 0, -160, 0, 120, 0, 0], [64, 0, -480, 0, 720, 0, -120]])
    elif m == 7:
        p = np.array([[2, 0, 0, 0, 0, 0, 0, 0], [4, 0, -2, 0, 0, 0, 0, 0], [8, 0, -12, 0, 0, 0, 0, 0],
                      [16, 0, -48, 0, 12, 0, 0, 0], [32, 0, -160, 0, 120, 0, 0, 0], [64, 0, -480, 0, 720, 0, -120, 0],
                      [128, 0, -1344, 0, 3360, 0, -1680, 0]])
    elif m == 8:
        p = np.array([[2, 0, 0, 0, 0, 0, 0, 0, 0], [4, 0, -2, 0, 0, 0, 0, 0, 0], [8, 0, -12, 0, 0, 0, 0, 0, 0],
                      [16, 0, -48, 0, 12, 0, 0, 0, 0], [32, 0, -160, 0, 120, 0, 0, 0, 0], [64, 0, -480, 0, 720, 0, -120, 0, 0],
                      [128, 0, -1344, 0, 3360, 0, -1680, 0, 0], [256, 0, -3584, 0, 13440, 0, -13440, 0, 1680]])
    elif m == 9:
        p = np.array([[2, 0, 0, 0, 0, 0, 0, 0, 0, 0], [4, 0, -2, 0, 0, 0, 0, 0, 0, 0], [8, 0, -12, 0, 0, 0, 0, 0, 0, 0],
                      [16, 0, -48, 0, 12, 0, 0, 0, 0, 0], [32, 0, -160, 0, 120, 0, 0, 0, 0, 0], [64, 0, -480, 0, 720, 0, -120, 0, 0, 0],
                      [128, 0, -1344, 0, 3360, 0, -1680, 0, 0, 0], [256, 0, -3584, 0, 13440, 0, -13440, 0, 1680, 0], [512, 0, -9216, 0, 48384, 0, -80640, 0, 30240, 0]])
    elif m == 10:
        p = np.array([[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [4, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0], [8, 0, -12, 0, 0, 0, 0, 0, 0, 0, 0],
                      [16, 0, -48, 0, 12, 0, 0, 0, 0, 0, 0], [32, 0, -160, 0, 120, 0, 0, 0, 0, 0, 0], [64, 0, -480, 0, 720, 0, -120, 0, 0, 0, 0],
                      [128, 0, -1344, 0, 3360, 0, -1680, 0, 0, 0, 0], [256, 0, -3584, 0, 13440, 0, -13440, 0, 1680, 0, 0], [512, 0, -9216, 0, 48384, 0, -80640, 0, 30240, 0, 0],
                      [1024, 0, -23040, 0, 161280, 0, -403200, 0, 302400, 0, -30240]])
    elif m == 11:
        p = np.array([[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [4, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0], [8, 0, -12, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      [16, 0, -48, 0, 12, 0, 0, 0, 0, 0, 0, 0], [32, 0, -160, 0, 120, 0, 0, 0, 0, 0, 0, 0], [64, 0, -480, 0, 720, 0, -120, 0, 0, 0, 0, 0],
                      [128, 0, -1344, 0, 3360, 0, -1680, 0, 0, 0, 0, 0], [256, 0, -3584, 0, 13440, 0, -13440, 0, 1680, 0, 0, 0], [512, 0, -9216, 0, 48384, 0, -80640, 0, 30240, 0, 0, 0],
                      [1024, 0, -23040, 0, 161280, 0, -403200, 0, 302400, 0, -30240, 0], [2048, 0, -56320, 0, 506880, 0, -1774080, 0, 2217600, 0, -665280, 0]])

    elif m == 12:
        p = np.array([[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [4, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      [8, 0, -12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [16, 0, -48, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0],
                      [32, 0, -160, 0, 120, 0, 0, 0, 0, 0, 0, 0, 0], [64, 0, -480, 0, 720, 0, -120, 0, 0, 0, 0, 0, 0],
                      [128, 0, -1344, 0, 3360, 0, -1680, 0, 0, 0, 0, 0, 0], [256, 0, -3584, 0, 13440, 0, -13440, 0, 1680, 0, 0, 0, 0],
                      [512, 0, -9216, 0, 48384, 0, -80640, 0, 30240, 0, 0, 0, 0], [1024, 0, -23040, 0, 161280, 0, -403200, 0, 302400, 0, -30240, 0, 0],
                      [2048, 0, -56320, 0, 506880, 0, -1774080, 0, 2217600, 0, -665280, 0, 0],
                      [4096, 0, -135168, 0, 1520640, 0, -7096320, 0, 13305600, 0, -7983360, 0, 665280]])

    # Assign Hermite coefficients for all dimensions
    for i in range(d):
        C[:, :m+1, i] = p
    
    return C


def legendre(m, d):
    """
    Constructs Legendre Polynomials much faster than legendreP.

    Args:
    - m (int): Degree of the Legendre polynomial.
    - d (int): Number of slices to create in the third dimension of the result.

    Returns:
    - C (numpy.ndarray): A 3D array containing the Legendre polynomials.
    """

    # Define the method for construction of the polynomials
    method = 'backward'  # 'forward' or 'backward'
    C = np.zeros((m, m+1, d))

    # Define the polynomial coefficients based on 'method' and degree 'm'
    if method == 'forward':
        if m == 1:
            p = np.array([[1, 0]])  # +3
        elif m == 2:
            p = np.array([[1, 0, 0], [1.5, 0, -0.5]])  # +4
        elif m == 3:
            p = np.array([[1, 0, 0, 0], [1.5, 0, -0.5, 0], [2.5, 0, -1.5, 0]])  # +5
        elif m == 4:
            p = np.array([[1, 0, 0, 0, 0], [1.5, 0, -0.5, 0, 0], [2.5, 0, -1.5, 0, 0],
                          [4.375, 0, -3.75, 0, 0.375]])  # +6
        elif m == 5:
            p = np.array([[1, 0, 0, 0, 0, 0], [1.5, 0, -0.5, 0, 0, 0], [2.5, 0, -1.5, 0, 0, 0],
                          [4.375, 0, -3.75, 0, 0.375, 0], [7.875, 0, -8.75, 0, 1.875, 0]])  # +7
        elif m == 6:
            p = np.array([[1, 0, 0, 0, 0, 0, 0], [1.5, 0, -0.5, 0, 0, 0, 0], [2.5, 0, -1.5, 0, 0, 0, 0],
                          [4.375, 0, -3.75, 0, 0.375, 0, 0], [7.875, 0, -8.75, 0, 1.875, 0, 0],
                          [14.4375, 0, -19.6875, 0, 6.5625, 0, -0.3125]])  # +8
        elif m == 7:
            p = np.array([[1, 0, 0, 0, 0, 0, 0, 0], [1.5, 0, -0.5, 0, 0, 0, 0, 0], [2.5, 0, -1.5, 0, 0, 0, 0, 0],
                          [4.375, 0, -3.75, 0, 0.375, 0, 0, 0], [7.875, 0, -8.75, 0, 1.875, 0, 0, 0],
                          [14.4375, 0, -19.6875, 0, 6.5625, 0, -0.3125, 0], [26.8125, 0, -43.3125, 0, 19.6875, 0, -2.1875, 0]])  # +9
        elif m == 8:
            p = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0], [1.5, 0, -0.5, 0, 0, 0, 0, 0, 0], [2.5, 0, -1.5, 0, 0, 0, 0, 0, 0],
                          [4.375, 0, -3.75, 0, 0.375, 0, 0, 0, 0], [7.875, 0, -8.75, 0, 1.875, 0, 0, 0, 0],
                          [14.4375, 0, -19.6875, 0, 6.5625, 0, -0.3125, 0, 0], [26.8125, 0, -43.3125, 0, 19.6875, 0, -2.1875, 0, 0],
                          [50.2734375, 0, -93.84375, 0, 54.140625, 0, -9.84375, 0, 0.2734375]])  # +10
        elif m == 9:
            p = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1.5, 0, -0.5, 0, 0, 0, 0, 0, 0, 0],
                          [2.5, 0, -1.5, 0, 0, 0, 0, 0, 0, 0], [4.375, 0, -3.75, 0, 0.375, 0, 0, 0, 0, 0],
                          [7.875, 0, -8.75, 0, 1.875, 0, 0, 0, 0, 0], [14.4375, 0, -19.6875, 0, 6.5625, 0, -0.3125, 0, 0, 0],
                          [26.8125, 0, -43.3125, 0, 19.6875, 0, -2.1875, 0, 0, 0], [50.2734375, 0, -93.84375, 0, 54.140625, 0, -9.84375, 0, 0.2734375, 0],
                          [94.9609375, 0, -201.09375, 0, 140.765625, 0, -36.09375, 0, 2.4609375, 0]])  # +11 
        elif m == 10:
            p = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1.5, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, 0],
                          [2.5, 0, -1.5, 0, 0, 0, 0, 0, 0, 0, 0], [4.375, 0, -3.75, 0, 0.375, 0, 0, 0, 0, 0, 0],
                          [7.875, 0, -8.75, 0, 1.875, 0, 0, 0, 0, 0, 0], [14.4375, 0, -19.6875, 0, 6.5625, 0, -0.3125, 0, 0, 0, 0],
                          [26.8125, 0, -43.3125, 0, 19.6875, 0, -2.1875, 0, 0, 0, 0], [50.2734375, 0, -93.84375, 0, 54.140625, 0, -9.84375, 0, 0.2734375, 0],
                          [94.9609375, 0, -201.09375, 0, 140.765625, 0, -36.09375, 0, 2.4609375, 0],
                          [180.4257813, 0, -427.3242188, 0, 351.9140625, 0, -117.3046875, 0, 13.53515625, 0]])  # +12
        elif m == 11:
            p = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1.5, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [2.5, 0, -1.5, 0, 0, 0, 0, 0, 0, 0, 0, 0], [4.375, 0, -3.75, 0, 0.375, 0, 0, 0, 0, 0, 0, 0],
                          [7.875, 0, -8.75, 0, 1.875, 0, 0, 0, 0, 0, 0, 0], [14.4375, 0, -19.6875, 0, 6.5625, 0, -0.3125, 0, 0, 0, 0, 0],
                          [26.8125, 0, -43.3125, 0, 19.6875, 0, -2.1875, 0, 0, 0, 0, 0],
                          [50.2734375, 0, -93.84375, 0, 54.140625, 0, -9.84375, 0, 0.2734375, 0],
                          [94.9609375, 0, -201.09375, 0, 140.765625, 0, -36.09375, 0, 2.4609375, 0],
                          [180.4257813, 0, -427.3242188, 0, 351.9140625, 0, -117.3046875, 0, 13.53515625, 0],
                          [344.4492188, 0, -902.1289063, 0, 854.6484375, 0, -351.9140625, 0, 58.65234375, 0]])   # +13
        elif m == 12:
            p = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [1.5, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [2.5, 0, -1.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [4.375, 0, -3.75, 0, 0.375, 0, 0, 0, 0, 0, 0, 0, 0],
                          [7.875, 0, -8.75, 0, 1.875, 0, 0, 0, 0, 0, 0, 0, 0],
                          [14.4375, 0, -19.6875, 0, 6.5625, 0, -0.3125, 0, 0, 0, 0, 0, 0],
                          [26.8125, 0, -43.3125, 0, 19.6875, 0, -2.1875, 0, 0, 0, 0, 0, 0],
                          [50.2734375, 0, -93.84375, 0, 54.140625, 0, -9.84375, 0, 0.2734375, 0],
                          [94.9609375, 0, -201.09375, 0, 140.765625, 0, -36.09375, 0, 2.4609375, 0],
                          [180.4257813, 0, -427.3242188, 0, 351.9140625, 0, -117.3046875, 0, 13.53515625, 0],
                          [344.4492188, 0, -902.1289063, 0, 854.6484375, 0, -351.9140625, 0, 58.65234375, 0],
                          [653.9140625, 0, -1672.5390625, 0, 1771.953125, 0, -840.34375, 0, 163.96484375, 0]])  # +14
    
    elif method == 'backward':
        # Start from case 12, then remove rows and columns for higher m values
        p = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1.5, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      [2.5, 0, -1.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [4.375, 0, -3.75, 0, 0.375, 0, 0, 0, 0, 0, 0, 0, 0],
                      [7.875, 0, -8.75, 0, 1.875, 0, 0, 0, 0, 0, 0, 0, 0], [14.4375, 0, -19.6875, 0, 6.5625, 0, -0.3125, 0, 0, 0, 0, 0, 0],
                      [26.8125, 0, -43.3125, 0, 19.6875, 0, -2.1875, 0, 0, 0, 0, 0, 0], [50.2734375, 0, -93.84375, 0, 54.140625, 0, -9.84375, 0, 0.2734375, 0, 0, 0, 0],
                      [94.9609375, 0, -201.09375, 0, 140.765625, 0, -36.09375, 0, 2.4609375, 0, 0, 0, 0],
                      [180.4257813, 0, -427.3242188, 0, 351.9140625, 0, -117.3046875, 0, 13.53515625, 0, -0.24609375, 0, 0],
                      [344.4492188, 0, -902.1289063, 0, 854.6484375, 0, -351.9140625, 0, 58.65234375, 0, 2.70703125, 0, 0],
                      [660.1943359, 0, -1894.470703, 0, 2029.790039, 0, -997.0898438, 0, 219.9462891, 0, -17.59570313, 0, 0.2255859375]])
        p = p[:m,:m+1]  # Adjust matrix size for the required m
        
    # Duplicate each slice in the third dimension
    for i in range(d):
        C[:m, :m+1, i] = p

    return C


def HDMR_EXT_construct_EB(EB, X, N, d, m, C, maxorder):
    """
    This function constructs the extended bases matrix, EB.

    Arguments:
    EB -- Extended bases matrix (output)
    X -- Input data (matrix of size N x d)
    N -- Number of samples
    d -- Number of dimensions
    m -- Number of basis functions
    C -- Coefficients for the orthogonal polynomials (3D array)
    maxorder -- Maximum order of interactions to consider

    Returns:
    EB -- Extended bases matrix
    O -- Indices of interactions
    """
    
    k = 0
    O = []
    
    # First-order terms (main effects)
    for i in range(d):
        for j in range(m):
            k += 1
#            EB[:, k-1] = orthopolyeval(C[j, :, i], X[:N, i], j)
            EB[:, k-1] = orthopolyeval(C[j, :, i], X[:N, i], j+1)
            O.append([j])

    if maxorder > 1:
        # Second-order terms (pairwise interactions)
        for i in range(d-1):
            for j in range(i+1, d):
                for z in range(m):
                    k += 1
#                    EB[:, k-1] = orthopolyeval(C[z, :, i], X[:N, i], z)
                    EB[:, k-1] = orthopolyeval(C[z, :, i], X[:N, i], z+1)
                    O.append([z])
                for z in range(m):
                    k += 1
#                    EB[:, k-1] = orthopolyeval(C[z, :, j], X[:N, j], z)
                    EB[:, k-1] = orthopolyeval(C[z, :, j], X[:N, j], z+1)
                    O.append([z])
                for z in range(m):
                    for r in range(m):
                        k += 1
#                        EB[:, k-1] = orthopolyeval(C[z, :, i], X[:N, i], z) * orthopolyeval(C[r, :, j], X[:N, j], r)
                        EB[:, k-1] = orthopolyeval(C[z, :, i], X[:N, i], z+1) * orthopolyeval(C[r, :, j], X[:N, j], r+1)
                        O.append([z, r])

    if maxorder == 3:
        # Third-order terms (three-way interactions)
        for i in range(d-2):
            for j in range(i+1, d-1):
                for z in range(j+1, d):
                    for r in range(m):
                        k += 1
#                        EB[:, k-1] = orthopolyeval(C[r, :, i], X[:N, i], r)
                        EB[:, k-1] = orthopolyeval(C[r, :, i], X[:N, i], r+1)
                        O.append([r])
                    for r in range(m):
                        k += 1
#                        EB[:, k-1] = orthopolyeval(C[r, :, j], X[:N, j], r)
                        EB[:, k-1] = orthopolyeval(C[r, :, j], X[:N, j], r+1)
                        O.append([r])
                    for r in range(m):
                        k += 1
#                        EB[:, k-1] = orthopolyeval(C[r, :, z], X[:N, z], r)
                        EB[:, k-1] = orthopolyeval(C[r, :, z], X[:N, z], r+1)
                        O.append([r])
                    for t in range(m):
                        for r in range(m):
                            k += 1
#                            EB[:, k-1] = orthopolyeval(C[t, :, i], X[:N, i], t) * orthopolyeval(C[r, :, j], X[:N, j], r)
                            EB[:, k-1] = orthopolyeval(C[t, :, i], X[:N, i], t+1) * orthopolyeval(C[r, :, j], X[:N, j], r+1)
                            O.append([t, r])
                    for t in range(m):
                        for r in range(m):
                            k += 1
#                            EB[:, k-1] = orthopolyeval(C[t, :, i], X[:N, i], t) * orthopolyeval(C[r, :, z], X[:N, z], r)
                            EB[:, k-1] = orthopolyeval(C[t, :, i], X[:N, i], t+1) * orthopolyeval(C[r, :, z], X[:N, z], r+1)
                            O.append([t, r])
                    for t in range(m):
                        for r in range(m):
                            k += 1
#                            EB[:, k-1] = orthopolyeval(C[t, :, j], X[:N, j], t) * orthopolyeval(C[r, :, z], X[:N, z], r)
                            EB[:, k-1] = orthopolyeval(C[t, :, j], X[:N, j], t+1) * orthopolyeval(C[r, :, z], X[:N, z], r+1)
                            O.append([t, r])
                    for i1 in range(m):
                        for i2 in range(m):
                            for i3 in range(m):
                                k += 1
#                                EB[:, k-1] = (orthopolyeval(C[i1, :, i], X[:N, i], i1) *
#                                               orthopolyeval(C[i2, :, j], X[:N, j], i2) *
#                                               orthopolyeval(C[i3, :, z], X[:N, z], i3))
                                EB[:, k-1] = (orthopolyeval(C[i1, :, i], X[:N, i], i1+1) *
                                               orthopolyeval(C[i2, :, j], X[:N, j], i2+1) *
                                               orthopolyeval(C[i3, :, z], X[:N, z], i3+1))
                                O.append([i1, i2, i3])

    return EB, O


def orthopolyeval(c, x, m):
    """
    This function evaluates the orthonormal polynomial for a given coefficient matrix.
    
    Arguments:
    c -- Coefficients of the polynomial (2D array with dimensions (1, m+1))
    x -- Input values (array or vector)
    m -- Degree of the polynomial
    
    Returns:
    phi -- The evaluated orthonormal polynomial
    """

    # Initialize the polynomial value with the first coefficient
#    phi = c[0, 0]  
    phi = c[0]

    # Evaluate the polynomial using the recurrence relation
    for i in range(1, m + 1):
     #   phi = x * phi + c[0, i]
         phi = x * phi + c[i]

    return phi

def HDMR_EXT_construct_B(EB, R, d, m, m2, m3, k1, k12, k, maxorder):
    # Preallocate memory for C2, C3, C_22, and C_33
    C_2 = np.full((2 * m + 1, m2), np.nan)
    C_22 = np.full((m2, m2, d * (d - 1) // 2), np.nan)
#    C_3 = np.full((3 * m + 3 * m**2, m3), np.nan)
    C_3 = np.full((3 * m + 3 * m**2 + 1, m3), np.nan)
    C_33 = np.full((m3, m3, d * (d - 1) * (d - 2) // 6), np.nan)
    
    # Compute C_0 and C
#    C_0 = np.sum(EB) / R
    C_0 = np.sum(EB,axis=0) / R
    C = np.dot(EB.T, EB) / R

    # Compute C_22
    l = 0

    for i in range(d - 1):
        for j in range(i + 1, d):
            l += 1
#            C_2[0, 0:m2] = C_0[k1 + (l - 1) * m2 + 1: k1 + l * m2]
            C_2[0, 0:m2] = C_0[k1 + (l - 1) * m2: k1 + l * m2]
            for z in range(2 * m):
#                C_2[z + 1, 0:m2] = C[k1 + (l - 1) * m2 + z, k1 + (l - 1) * m2 + 1: k1 + l * m2]
                C_2[z + 1, 0:m2] = C[k1 + (l - 1) * m2 + z, k1 + (l - 1) * m2: k1 + l * m2]
            C_22[:, :, l - 1] = np.dot(C_2.T, C_2)

    # Compute C_33 if maxorder == 3
    if maxorder == 3:
        l = 0
        for i in range(d - 2):
            for j in range(i + 1, d - 1):
                for z in range(j + 1, d):
                    l += 1
                    C_3[0, 0:m3] = C_0[k12 + (l - 1) * m3: k12 + l * m3]
                    for ii in range(3 * m + 3 * m**2):
#                        C_3[ii + 1, 0:m3] = C[k12 + (l - 1) * m3 + ii, k12 + (l - 1) * m3 + 1: k12 + l * m3]
                        C_3[ii + 1, 0:m3] = C[k12 + (l - 1) * m3 + ii, k12 + (l - 1) * m3: k12 + l * m3]
                    C_33[:, :, l - 1] = np.dot(C_3.T, C_3)
    else:
        C_33 = np.zeros_like(C_33)

    # Construct B matrix
    B = np.zeros((k, k))
    
    # Fill in C_22
    for l in range(d * (d - 1) // 2):
#        B[k1 + (l - 1) * m2: k1 + l * m2, k1 + (l - 1) * m2: k1 + l * m2] = C_22[:, :, l]
        B[k1 + (l) * m2: k1 + (l+1) * m2, k1 + (l) * m2: k1 + (l + 1) * m2] = C_22[:, :, l]

    # Fill in C_33 if maxorder == 3
    if maxorder == 3:
        for l in range(d * (d - 1) * (d - 2) // 6):
#            B[k12 + (l - 1) * m3: k12 + l * m3, k12 + (l - 1) * m3: k12 + l * m3] = C_33[:, :, l]
            B[k12 + (l) * m3: k12 + (l+1) * m3, k12 + (l) * m3: k12 + (l + 1) * m3] = C_33[:, :, l]

    return B, C

def HDMR_EXT_dmorph(EB, B, C, Y_res, R, k1, k):

    # Remove all 1st order component functions
    C = C[k1:, :]
    # Right hand side of the algebraic equation
    D = np.dot(EB.T, Y_res) / R
#    D = D[k1:, :]
    D = D[k1:] 

    # Generalized inverse of coefficient matrix
    C_inv = np.linalg.pinv(C)
    # Least squares solution
    Cf_0 = np.dot(C_inv, D)
    # Cf_0 is correct

    Ia = np.eye(k)
    Pr = Ia - np.dot(C_inv, C)
    PB = np.dot(Pr, B)

    # Singular value decomposition
#    U, S, Vt = np.linalg.svd(PB)  # Vt is the transpose of V
#    s = np.diag(S)  # Extract vector with singular values [not needed in python]
    U, s, Vt = np.linalg.svd(PB)  # Vt is the transpose of V
    
    # Normalized difference of adjacent singular values
    diff_s = -np.diff(s) / s[1:k]  # diff(s) divided by s[1:k]
    tol = 0  # Initialize tolerance
    n_ns = int(1)  # Index non-singular values not satisfying tol
    
    # Determine the number of non-singular values based on the tolerance
    for i in range(k - 1):
        if diff_s[i] > tol:
            tol = diff_s[i]
            n_ns = i + 1  # Adjust for 0-based index

    V = Vt.T # transpose back the Vt matrix [JAV]
    # Remove columns of singular values from U and V
    U = U[:, n_ns:]
    #Vt = Vt[n_ns:, :]
    V = V[:, n_ns:]

    # Compute UV and Q matrices
    UV = np.dot(U.T, V)
    #Q = np.dot(Vt.T, np.linalg.pinv(UV))
    #Q = np.dot(V, np.linalg.pinv(UV))
    Q = V @ np.linalg.pinv(UV) @ U.T
    # Compute the D-Morph Regression result
    Cf = np.dot(Q, Cf_0)
    
    return Cf, n_ns


def HDMR_EXT_comp_func(EB, w, R, d, m, k1, k2, k3, k, n1, n2, n3, maxorder):
    # Ensure that w is a row vector (in MATLAB: w(:)' turns w into a row vector)
    w = np.array(w).flatten()  # Convert to a 1D array
    w = w.reshape(1, -1)  # Ensure it's a row vector
    
    # Pre-allocate memory
    Y_t = np.zeros((R, k))
#   Y_em = np.zeros((R, d))
    if maxorder == 1:
        Y_em = np.zeros((R, n1))
    elif maxorder == 2:
        Y_em = np.zeros((R, n1 + n2))
    else:
        Y_em = np.zeros((R, n1 + n2 + n3))
    
    # Compute temporary matrix Y_t
    for i in range(R):
        Y_t[i, :] = EB[i, :] * w
    
    # First order component functions
    for i in range(n1):
#        Y_em[:R, i] = np.sum(Y_t[:R, i*m - m + 1:i*m], axis=1)
        Y_em[:R, i] = np.sum(Y_t[:R, (i+1)*m - m:(i+1)*m], axis=1)
    
#for i = 1:n1
#    Y_em(1:R,i) = sum(Y_t(1:R,i*m-m+1:i*m),2);
#end
#MATLAB: i = 1, 1:m, i = 2, m+1: 2*m
#PYTHON: i = 0 --> 0:m
#        i = 1 --> m:2*m

    # Second order component functions
    if maxorder > 1:
        m = k2 // n2  # Update m based on k2 and n2
        for i in range(n2):
#            Y_em[:R, i + n1] = np.sum(Y_t[:, i*m - m + k1 + 1:i*m + k1], axis=1)
            Y_em[:R, i + n1] = np.sum(Y_t[:, (i+1)*m - m + k1:(i+1)*m + k1], axis=1)
    
    # Third order component functions
    if maxorder == 3:
        m = k3 // n3  # Update m based on k3 and n3
        for i in range(n3):
#            Y_em[:R, i + n1 + n2] = np.sum(Y_t[:, i*m - m + k1 + k2 + 1:i*m + k1 + k2], axis=1)
            Y_em[:R, i + n1 + n2] = np.sum(Y_t[:, (i+1)*m - m + k1 + k2:(i+1)*m + k1 + k2], axis=1)

    return Y_em


def HDMR_EXT_F_test(y, f0, Y_bf, R, alfa, m1, m2, m3, n1, n2, n3, n, method):
    # Initialize ind with zeros (all terms insignificant)
    ind = np.zeros(n)
    
    # Determine the significant components of the HDMR model via the F-test
    if method == 1:  # Forward selection of terms (start with Y = f0)
        Y_res0 = y - f0
        SSR0 = np.sum(Y_res0**2)
        p0 = 0
        
        for i in range(n):
            # Model with the ith term included
            Y_res1 = Y_res0 - Y_bf[:, i]
            
            # Number of parameters of the proposed model (order dependent)
            if i < n1:
                p1 = p0 + m1  # 1st order
            elif i >= n1 and i < n1 + n2:
                p1 = p0 + m2  # 2nd order
            else:
                p1 = p0 + m3  # 3rd order
            
            # Calculate SSR of Y1
            SSR1 = np.sum(Y_res1**2)
            
            # Now calculate the F_stat (F_stat > 0 -> SSR1 < SSR0)
            F_stat = ((SSR0 - SSR1) / (p1 - p0)) / (SSR1 / (R - p1))
            
            # Now calculate critical F value at confidence level 1-alfa
            F_crit = f.ppf(1 - alfa, p1 - p0, R - p1)
            
            # Now determine whether to accept the ith component into the model
            if F_stat > F_crit:
                # ith term is significant and should be included in model
                ind[i] = 1
                Y_res0 = Y_res1
                SSR0 = SSR1
                p0 = p1

    elif method == 2:  # Backward elimination of terms (start with Y = f0 + sum(all_terms))
        # NOTE: ONLY POSSIBLE IF R - p1 > 0 OTHERWISE MUST DO FORWARD SELECTION
        Y_res0 = y - f0 - np.sum(Y_bf, axis=1)
        SSR0 = np.sum(Y_res0**2)
        
        # Determine number of parameters of full model
        p0 = n1 * m1 + n2 * m2 + n3 * m3
        
        for i in range(n):
            # Previous model with ith term excluded
            Y_res1 = Y_res0 + Y_bf[:, i]
            
            # Number of parameters of the proposed model (order dependent)
            if i < n1:
                p1 = p0 - m1  # 1st order
            elif i >= n1 and i < n1 + n2:
                p1 = p0 - m2  # 2nd order
            else:
                p1 = p0 - m3  # 3rd order
            
            # Calculate SSR of Y1
            SSR1 = np.sum(Y_res1**2)
            
            # Now calculate the F_stat (F_stat > 0 if SSR1 > SSR0 as p1 < p0)
            F_stat = ((SSR0 - SSR1) / (p1 - p0)) / (SSR1 / (R - p1))
            
            # Now calculate critical F value at confidence level 1 - alfa
            F_crit = f.ppf(1 - alfa, p0 - p1, R - p1)
            
            # Now determine whether ith component is significant
            if F_stat > F_crit:
                # ith term is significant and should retain in model
                ind[i] = 1
            else:
                # ith term is insignificant and will be removed from model
                p0 = p1
                SSR0 = SSR1
                Y_res0 = Y_res1

    # Now return the number of terms of the final HDMR model
    nterms = np.sum(ind)
    rmse = np.sqrt(np.sum((y - (f0 + np.sum(Y_bf[:, ind == 1], axis=1)))**2) / R)

    return ind, nterms, p0, rmse


def ANCOVA(y, Y_em, V_y, R, n):
    """
    ANALYSIS OF COVARIANCE: Returns sensitivity indices for each emulator term.

    Parameters:
        y (numpy array): R x 1 vector of the model output for each row of X.
        Y_em (numpy array): R x n matrix of HDMR model terms from backfitting.
        V_y (float): Scalar with the variance of the model output.
        R (int): Number of samples of X.
        n (int): Number of terms in the HDMR model.

    Returns:
        S (numpy array): n x 1 vector of total sensitivity for each emulator term.
        S_a (numpy array): n x 1 vector of structural sensitivity for each emulator term.
        S_b (numpy array): n x 1 vector of correlative sensitivity for each emulator term.
        V_em (numpy array): n x 1 vector of variance for each emulator term.
    """
    # Now compute the sum of all Y_em terms (i.e., sum over the rows for each column)
    Y0 = np.sum(Y_em[:R, :n], axis=1)
    
    # Initialize each variable
    S = np.full((n, 1), np.nan)
    S_a = np.full((n, 1), np.nan)
    S_b = np.full((n, 1), np.nan)
    V_em = np.full((n, 1), np.nan)

    # Analysis of covariance -> extension of analysis of variance
    for j in range(n):
        # Covariance matrix of jth term of Y_em and actual Y
        C = np.cov(Y_em[:R, j], y, rowvar=False)
        S[j, 0] = C[0, 1] / V_y  # Total sensitivity of jth term (Eq. 19 of Li et al)
        
        # Covariance matrix of jth term with emulator Y without jth term
        C = np.cov(Y_em[:R, j], Y0 - Y_em[:R, j], rowvar=False)
        S_a[j, 0] = C[0, 0] / V_y  # Structural contribution of jth term (Eq. 20 of Li et al)
        S_b[j, 0] = C[0, 1] / V_y  # Correlative contribution of jth term (Eq. 21 of Li et al)
        
        # Variance in Y of jth term (S_a * V_y + S_b * V_y)
        V_em[j, 0] = np.sum(C[0, :2])  # (S_a + S_b) * V_y = C(0,1) from Li et al Eq. 37

    S = S.reshape(-1)
    S_a = S_a.reshape(-1)
    S_b = S_b.reshape(-1)
    V_em = V_em.reshape(-1)

    return S, S_a, S_b, V_em


def HDMR_EXT_end(SA, nterms, p0, RMSE, select, K, C2, C3, n1, n2, n3, n, n_ns):
    """
    Prepares the return arguments of HDMR_EXT.
    
    Arguments:
    SA -- Sensitivity Analysis table (Matrix)
    nterms -- Number of terms in the model
    p0 -- Initial guess for the model parameters
    RMSE -- Root Mean Square Error
    select -- Selected terms for bootstrap sampling
    K -- Number of bootstrap trials
    C2, C3 -- Some specific constants (not used explicitly in the original code)
    n1, n2, n3 -- Parameters related to term breakdown
    n -- Number of terms in the HDMR model
    n_ns -- Not used explicitly in the original code
    
    Returns:
    SA -- Updated Sensitivity Analysis table
    SA_sig -- Table with only significant sensitivity terms
    Fx -- Other sensitivity results (e.g., term coefficients)
    """

    # Set the alpha (significance level) for bootstrap confidence intervals
    alfa = 0.05
    
    # ----------------------------------------------------------------------- #
    # Now compile table, SI, with all sensitivity results
    # ----------------------------------------------------------------------- #
    
    # Take sum of all selected terms of all K bootstrap trials
    f_ret = np.sum(select, axis=1)
    
    # Now construct two tables with results
    SA, Fx = construct_tables(SA, nterms, p0, RMSE, alfa, f_ret, K, C2, C3, n1, n2, n3, n, n_ns)
    
    # Compute some variables used to write tables
    nr = len(SA)
    id_sig = np.concatenate(([1], np.where(f_ret > 0)[0] + 2, [nr])) - 1
    SA_sig = [SA[i][:13] for i in id_sig]    
    
    # ----------------------------------------------------------------------- #
    # Print a Table with emulator results
    # ----------------------------------------------------------------------- #
    
    with open('HDMR_EXT_results.txt', 'w') as fid:
        # Print header (equivalent to the MATLAB block of ASCII art)
        fid.write('|-------------------------------------------------------------------------|\n')
        fid.write('|                                                                         |\n')
        fid.write('| HHH   HHH DDDDDDDD  MMM    MMM RRRRRRR     EEEEEEEE XXX   XXX TTTTTTTTT |\n')
        fid.write('| HHH   HHH DDDDDDDDD MMM    MMM RRRRRRRR    EEEEEEEE XXX   XXX TTTTTTTTT |\n')
        fid.write('| HHH   HHH DDD   DDD MMM    MMM RRR   RRR   EEE       XXX XXX  TT TTT TT |\n')
        fid.write('| HHH   HHH DDD   DDD MMMM  MMMM RRR   RRR   EEE       XXX XXX  T  TTT  T |\n')
        fid.write('| HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRRR    EEEEEE     XXXXX      TTT    |\n')
        fid.write('| HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRR     EEEEEE     XXXXX      TTT    |\n')
        fid.write('| HHH   HHH DDD   DDD MMM    MMM RRRRRRR     EEE       XXX XXX     TTT    |\n')
        fid.write('| HHH   HHH DDD   DDD MMM    MMM RRR  RRR    EEE       XXX XXX     TTT    |\n')
        fid.write('| HHH   HHH DDDDDDDDD MMM    MMM RRR   RRR   EEEEEEEE XXX   XXX    TTT    |\n')
        fid.write('| HHH   HHH DDDDDDDD  MMM    MMM RRR   RRR   EEEEEEEE XXX   XXX    TTT    |\n')
        fid.write('|                                                                         |\n')
        fid.write('|-------------------------------------------------------------------------|\n')
        fid.write('\n')

        # Print table header
        if K == 1:
            fid.write('Table 1: Properties of HDMR_EXT emulator, y = f(x), for training data set (no bootstrapping)\n')
        else:
            fid.write(f'Table 1: Properties of HDMR_EXT emulator, y = f(x), for randomized training data set of {K} bootstrap trials\n')
        fid.write('===============================================\n')
        fid.write('Emulator    # terms   # coefs.   RMSE     ns   \n')
        fid.write('-----------------------------------------------\n')

        # Format for printing table
        fmt_1 = '   {:<7}    {:>3}       {:>3}     {:>6.3f}   {:>3} \n'
        for k in range(1, K + 1):
            # Assuming Fx is a list of lists, print based on the structure
            fid.write(fmt_1.format(int(Fx[k][0]), Fx[k][1], Fx[k][2], float(Fx[k][3]), Fx[k][4]))
        fid.write('===============================================\n')
        fid.write('\n')

        # Print sensitivity analysis results
        for pr_tab in range(1, 3):
            fid.write('\n\n')
            if pr_tab == 1:
                id_table = np.arange(nr)
                SAtable = SA
            else:
                id_table = id_sig
                nr = len(id_sig)
                SAtable = SA_sig

            # Print header for each table
            if K == 1:
                if pr_tab == 1:
                    fid.write('Table 2: HDMR_EXT results for all model components using all X and Y data (no bootstrapping)\n')
                else:
                    fid.write('Table 3: HDMR_EXT results for significant model components only using all X and Y data (no bootstrapping)\n')
            elif K > 1:
                if pr_tab == 1:
                    fid.write(f'Table 2: HDMR_EXT results for all model components using {K} bootstrap trials\n')
                else:
                    fid.write(f'Table 3: HDMR_EXT results for significant model components only using {K} bootstrap trials\n')
            
            fid.write('============================================================================================================ \n')
            fid.write('                                         ( = JOINT DETERMINATION )                                    \n')
            fid.write('                  ------------------------------------------------------------------------------------------ \n')
            fid.write('Term       Order        S^a             S^b              S               ST          V(i)/V(Y)     #select   \n')
            fid.write('------------------------------------------------------------------------------------------------------------ \n')

            # Define formats for printing the table
            if K == 1:
                fmt_1 = '{:<11}  {:>1}     {:>6.3f} (-----)  {:>6.3f} (-----)  {:>6.3f} (-----)  {:>6.3f} (-----)  {:>6.3f} (-----)    {:<2}\n'
                fmt_2 = '{:<11}  {:>1}     {:>6.3f} (-----)  {:>6.3f} (-----)  {:>6.3f} (-----)         (-----)  {:>6.3f} (-----)    {:<2}\n'
                
                for r in range(1, nr - 1):
                    if id_table[r - 1] >= n1: # columns 3,5,7,11,13 in matlab = 2,4,6,10,12 in python
                        fid.write(fmt_2.format(str(SAtable[r][0]), int(SAtable[r][1]), safe_float(SAtable[r][2]), safe_float(SAtable[r][4]), safe_float(SAtable[r][6]), safe_float(SAtable[r][10]), safe_int(SAtable[r][12])))
                    else: # columns 3,5,7,9,11,13 in matlab = 2,4,6,8,10,12 in python 
                        fid.write(fmt_1.format(str(SAtable[r][0]), int(SAtable[r][1]), safe_float(SAtable[r][2]), safe_float(SAtable[r][4]), safe_float(SAtable[r][6]), safe_float(SAtable[r][8]), safe_float(SAtable[r][10]), safe_int(SAtable[r][12])))
            else:
                fmt_1 = '{:<11}  {:>1}     {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})    {:<2}\n'
                fmt_2 = '{:<11}  {:>1}     {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})                  {:>6.3f} (\261{:>.2f})    {:<2}\n'

                for r in range(1, nr - 1):
                    if id_table[r - 1] >= n1: # columns 3:8,11:13 in matlab = 2:7,10,11,12 in python
                        fid.write(fmt_2.format(str(SAtable[r][0]), int(SAtable[r][1]), safe_float(SAtable[r][2]), abs(safe_float(SAtable[r][3])), safe_float(SAtable[r][4]), abs(safe_float(SAtable[r][5])), safe_float(SAtable[r][6]), abs(safe_float(SAtable[r][7])), safe_float(SAtable[r][10]), abs(safe_float(SAtable[r][11])), safe_int(SAtable[r][12])))
                    else: # columns 3:13 in matlab = 2:12 in python
                        fid.write(fmt_1.format(str(SAtable[r][0]), int(SAtable[r][1]), safe_float(SAtable[r][2]), abs(safe_float(SAtable[r][3])), safe_float(SAtable[r][4]), abs(safe_float(SAtable[r][5])), safe_float(SAtable[r][6]), abs(safe_float(SAtable[r][7])), safe_float(SAtable[r][8]), abs(safe_float(SAtable[r][9])), safe_float(SAtable[r][10]), abs(safe_float(SAtable[r][11])), safe_int(SAtable[r][12])))

            fid.write('------------------------------------------------------------------------------------------------------------ \n')
            
            if K == 1: # columns 3,5,7,11 in matlab = 2,4,6,10 in python
                fmt = '{:<11}  {:>1}     {:>6.3f} (-----)  {:>6.3f} (-----)  {:>6.3f} (-----)         (-----)  {:>6.3f} (-----)\n'
#               data = np.array([SAtable[nr-1][i] for i in range(2, 10, 2)])
                fid.write(fmt.format(str(SAtable[nr-1][0]), '', safe_float(SAtable[nr-1][2]), safe_float(SAtable[nr-1][4]), safe_float(SAtable[nr-1][6]), safe_float(SAtable[nr-1][10])))
            else: # columns 3:8,11:12 in matlab = 2:7,10,11 in python
                fmt = '{:<11}  {:>1}     {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})                  {:>6.3f} (\261{:>.2f}) \n'
                fid.write(fmt.format(str(SAtable[nr - 1][0]), '', safe_float(SAtable[nr-1][2]), abs(safe_float(SAtable[nr-1][3])), safe_float(SAtable[nr-1][4]), abs(safe_float(SAtable[nr-1][5])), safe_float(SAtable[nr-1][6]), safe_float(SAtable[nr-1][7]), safe_float(SAtable[nr-1][10]), safe_float(SAtable[nr-1][11])))

            fid.write('============================================================================================================ \n')
            fid.write(' S^a: Structural sensitivity index of individual terms\n')
            fid.write(' S^b: Correlative sensitivity index of individual terms\n')
            fid.write(' S: Total sensitivity index of individual terms\n')
            fid.write(' ST: Total sensitivity index\n')
            if K == 1:
                fid.write(' (--): Cannot compute confidence intervals of listed statistics with K = 1\n')
            else:
                fid.write(f' (\261): {int(100*(1-alfa))}% confidence intervals derived from bootstrapping\n')
                fid.write(' V(i)/V_Y: Relative contribution of each term to model output variance ( = var(Y))\n')
            if K == 1:
                fid.write(' #select: 0 (if term is insignificant) and 1 (significant term)\n')
            else:
                fid.write(' #select: Number of bootstrap trials that identifies respective term as significant\n')

            fid.write('\n')

    if sys.platform == "win32" or sys.platform == "darwin":
        os.system("start HDMR_EXT_settings.txt")
        os.system("start HDMR_EXT_results.txt")

    return SA, SA_sig, Fx


def safe_int(value):
    return int(value) if value else 0  # or return a default value if empty


def safe_float(value):
    return float(value) if value else 0.0  # or return a default value if empty


def construct_tables(SA, nterms, p0, RMSE, alfa, f_ret, K, C2, C3, n1, n2, n3, n, n_ns):
    # This function returns a Table with a) sensitivity, and b) emulator results
    # Same as contruct_tables of HDMR package except n_ns is extra input argument
 
    # ----------------------------------------------------------------------- #
    #                 FIRST POSTPROCESS SENSITIVITY ESTIMATES                 #
    # ----------------------------------------------------------------------- #

    # Compute average sensitivity values
    S_m = np.mean(SA['S'], axis=1)
    Sa_m = np.mean(SA['Sa'], axis=1)
    Sb_m = np.mean(SA['Sb'], axis=1)
    
    # Now calculate sum of each statistic
    S_sum = np.sum(SA['S'], axis=0)
    Sa_sum = np.sum(SA['Sa'], axis=0)
    Sb_sum = np.sum(SA['Sb'], axis=0)
    
    # Now calculate associated std's
    Z = lambda p: -np.sqrt(2) * erfcinv(p * 2)
    m = Z(1 - alfa/2)  # Multiplier, alfa is significance level   

    # Compute output statistics for Y (variance)
    V_em_div_V_y = SA['V_em'] / SA['V_y']
    V_em_rat = np.mean(V_em_div_V_y, axis=1)
    V_em_div_V_Y_sum = np.sum(V_em_div_V_y, axis=0)
    
    # Compute standard deviation of bootstrap results
    if K > 1:
        S_range = m * np.std(SA['S'], axis=1, ddof=0)
        Sa_range = m * np.std(SA['Sa'], axis=1, ddof=0)
        Sb_range = m * np.std(SA['Sb'], axis=1, ddof=0)
        V_em_range = m * np.std(V_em_div_V_y, axis=1, ddof=0)
        S_sum_range = m * np.std(S_sum)
        Sa_sum_range = m * np.std(Sa_sum)
        Sb_sum_range = m * np.std(Sb_sum)
        V_em_div_V_Y_sum_range = m * np.std(V_em_div_V_Y_sum)
    else:
        S_range, Sa_range, Sb_range, V_em_range = [np.nan] * n, [np.nan] * n, [np.nan] * n, [np.nan] * n

    #ST = np.full(n1, np.nan)
    #ST_range = np.full(n1, np.nan)
    ST, ST_range = [np.nan] * n1, [np.nan] * n1

    # Now compute the total sensitivity of each parameter/coefficient
    for r in range(n1):
        ij = n1 + np.where(np.sum(C2 == r, axis=1) == 1)[0]
        if not isinstance(C3, (int, float)):
            ijk = n1 + n2 + np.where(np.sum(C3 == r, axis=1) == 1)[0]
            TS = np.sum(SA['S'][np.r_[r, ij, ijk], :K], axis=0)
        else:
            TS = np.sum(SA['S'][np.r_[r, ij], :K], axis=0)

        # use all bootstrap trials to determine total sensitivity + range!
        ST[r] = np.mean(TS)
        if K > 1:
            ST_range[r] = m * np.std(TS)

    # ----------------------------------------------------------------------- #
    #               NOW CONSTRUCT TABULATED TABLE WITH RESULTS                #
    # ----------------------------------------------------------------------- #

    # how many rows of this table
    nr = n1 + n2 + n3 + 1
    
    # initialize row_names of Table and f_i ( = order)
    row_names = [''] * nr
    f_ord = np.full(nr, np.nan)
    
    # now create row_names + f_i
    for r in range(n1):
        f_ord[r] = 1
        row_names[r] = f'x{r+1}'

    for i in range(n2):
        r = i + n1
        f_ord[r] = 2
        row_names[r] = f'x{C2[i, 0]+1}/x{C2[i, 1]+1}'

    for i in range(n3):
        r = i + n1 + n2
        f_ord[r] = 3
        row_names[r] = f'x{C3[i, 0]+1}/x{C3[i, 1]+1}/x{C3[i, 2]+1}'
   
    # add as last row name the sum of the previous rows
    row_names[nr-1] = 'sum'
    
    # now create column names
    col_names = ['term', 'order', 'S_a', 'std.', 'S_b', 'std.', 'S', 'std.', 'S_T', 'std.', 'V(i)/V(Y)', 'std.', '#select']
    nc = len(col_names)

    # Reinitialize SA to become a list of sensitivity analysis results
    SA_table = [[''] * nc for _ in range(nr+1)]
    SA_table[0][:] = col_names
    
    # first column stores row_names
    for i in range(nr):
        SA_table[i+1][0] = row_names[i]

    # Fill columns 2 - 8 of table -> [ f_ord, S_a, std., S_b, std., S, std. ]
    j = [0, 1, 2, 3, 4, 5]
    for r in range(1, nr):
        SA_table[r][1] = int(f_ord[r-1])
        data = np.array([Sa_m[r-1], Sa_range[r-1], Sb_m[r-1], Sb_range[r-1], S_m[r-1], S_range[r-1]])
        for idx in range(len(j)):
            SA_table[r][2 + idx] = f'{data[j[idx]]:4.3f}'

    # Fill columns 9 - 10 of table -> [ ST, std. ]
    j = [0, 1]
    for r in range(1, n1+1):
        data = np.array([ST[r-1], ST_range[r-1]])
        for idx in range(len(j)):
            SA_table[r][8 + idx] = f'{data[j[idx]]:4.3f}'

    # Fill columns 11 - 12 of table -> [ V(i)/V_Y, std.]
    j = [6, 7]
    for r in range(1, nr):
        data = np.array([V_em_rat[r-1], V_em_range[r-1]])
        for idx in range(len(j)):
            SA_table[r][10 + idx] = f'{data[j[idx]-6]:4.3f}'

    # Fill column 13 of table -> [ FX.select ]
    for r in range(1, nr):
        SA_table[r][12] = int(f_ret[r-1])

    # Fill the last row of table (sum)
    if K > 1:
        in_T = [
            np.sum(Sa_m), Sa_sum_range, np.sum(Sb_m), Sb_sum_range, np.sum(S_m), S_sum_range,
            np.nan, np.nan, np.sum(V_em_rat), V_em_div_V_Y_sum_range, np.nan
        ]
    elif K == 1:
        in_T = [
            np.sum(Sa_m), np.nan, np.sum(Sb_m), np.nan, np.sum(S_m), np.nan, np.nan, np.nan,
            np.sum(V_em_rat), np.nan, np.nan
        ]
    
    j = 2 + np.where(~np.isnan(in_T))[0]
    for idx in range(len(j)):
        SA_table[nr][j[idx]] = f'{in_T[j[idx]-2]:4.3f}'

    # ----------------------------------------------------------------------- #
    #               NOW CONSTRUCT TABULATED RESULTS EMULATORS                 #
    # ----------------------------------------------------------------------- #

    row_names = [str(r+1) for r in range(K)]
    col_names = ['emulator', '# terms', '# coefs.', 'RMSE', '# nonsing.']
    nc = len(col_names)

    Fx = [[''] * nc for _ in range(K+1)]
    Fx[0][:] = col_names
    j = [0, 1]
    for k in range(1, K+1):
        Fx[k][0] = int(k)
        data = np.array([nterms[k-1], p0[k-1]])
        for idx in range(len(j)):
           # Fx[k][idx+1] = f'{data[j[idx]]:4.3f}'
            Fx[k][idx+1] = int(data[j[idx]])

        Fx[k][nc-1] = int(n_ns[k-1])
        Fx[k][nc-2] = f'{RMSE[k-1]:4.3f}'

    return SA_table, Fx

def HDMR_EXT_plot(SA_sig, Fx, y, Y_e, select, p0, n_ns, id, R, K):
    # Print wait statement to the screen
    print('HDMR_EXT PLOTTING: PLEASE WAIT ...')
    
    # Define name of program
    n_program = 'HDMR_EXT'
    
    # Define name of figures file
    file_name = f'{n_program}_figures.pdf'
    
    # Define names of tables
    table1_name = f'Table_1_{n_program}.pdf'
    table2_name = f'Table_2_{n_program}.pdf'
    
    # Determine the entries without pdf
    id_pdf = range(1, table1_name.find('pdf') - 2)
    
    # Determine screen size (using matplotlib to get screen dimensions)
    monitor = get_monitors()[0]
    screen_width = monitor.width
    screen_height = monitor.height
    x_mult = screen_width / 1920
    y_mult = screen_height / 1080
    t_mult = min(x_mult, y_mult)

    # Define fontsize for figures
    fontsize_xylabel = 18 * t_mult
    fontsize_axis = 16 * t_mult
    fontsize_legend = 14 * t_mult
    fontsize_text = 14 * t_mult
    fontsize_table = 16 * t_mult
    fontsize_titlepage = 30 * t_mult
    
    # ----------------------------------------------------------------------- %
    # Now plot empty figure for PDF file
    # ----------------------------------------------------------------------- %

    with PdfPages(file_name) as pdf:
    
        ### Plot Empty Figure for PDF
        plt.figure(figsize=(12, 6))
        plt.plot([], [], 'ro')  # Empty plot
        plt.axis([0, 1, 0, 1])
        plt.gca().set_facecolor('w')
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])        
        plt.gca().set_xticks([])
        plt.gca().set_yticks([])
        plt.text(0.1 * x_mult, 0.6 * y_mult, r'${\rm Visual \; results \; of \; HDMR_{EXT} \; toolbox}$', fontsize=fontsize_titlepage) #, ha='center', va='center')
        plt.text(0.11 * x_mult, 0.5 * y_mult, r'${\rm Tables \; are \; printed \; to \; this \; PDF \; file}$', fontsize=fontsize_titlepage) #, ha='center', va='center') #, fontweight='bold')
        ax = plt.gca()  # Get current axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        pdf.savefig()    
        plt.show()

        #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        # Now plot results of emulator of K trials
        #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        Y_plot = np.linspace(min(y), max(y))
    
        # Number of rows and columns for subplots
        row = 2
        col = 4
    
        for k in range(K):
            if (k % (row * col) == 0):
                plt.figure(figsize=(15, 10))
                c = 1
                ii = 1
            else:
                ii = ii+1

            if c <= col:
                r = 1
                plot_c = c
            else:
                r = 2
                plot_c = c - col

            id_R = id[0:R, k]
            ax1 = plt.subplot(row, col, ii)
            ax1.plot(y[id_R], Y_e[0:R, k], 'rs', markerfacecolor='r', markeredgecolor='r')
            # Now plot the 1:1 Line in gray
            ax1.plot(Y_plot, Y_plot, '-', color=[0.5, 0.5, 0.5], linewidth=2)
            ax1.set_xlabel(r"$y$", fontsize=fontsize_xylabel)
            if k % col == 0:
                ax1.set_ylabel(r"$y = f(\mathbf{x})$", fontsize=fontsize_xylabel)
            ax1.legend(['Training', '1:1 Line'], loc='lower right', fontsize=fontsize_legend, frameon=False)
            ax1.tick_params(axis='both', labelsize=fontsize_axis)

            RMSE_train = np.sqrt(1 / R * np.sum((y[id_R] - Y_e[0:R, k]) ** 2))
            RMSE_1 = f'{RMSE_train:.4f}'
            ax1.set_title(f'Bootstrap trial: {k + 1}', fontweight='bold', fontsize=fontsize_text)

            ax1.text(0.05, 0.92, f'$\# terms.: {int(np.sum(select[:, k]))}$', transform=ax1.transAxes, fontsize=fontsize_text)
            ax1.text(0.05, 0.85, f'$\# coef.: {int(p0[k])}$', transform=ax1.transAxes, fontsize=fontsize_text)
            ax1.text(0.05, 0.78, f'$\# nonsgl.: {int(n_ns[k])}$', fontsize=fontsize_text, transform=ax1.transAxes)
            ax1.text(0.05, 0.70, r"$\text{RMSE}_\text{T}:$"f'{RMSE_train:.4f}', transform=ax1.transAxes, fontsize=fontsize_text)
            c += 1
            if ii == row*col or k == K-1:
                pdf.savefig()
                plt.close()

        plt.close()

        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        # Now plot Fx to table in figure
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        col_names = Fx[0]
        n_col = len(col_names)

        # Format column names
        col_names = [f'<html><center /><font size=5>{name}</font></html>' for name in col_names]

    #   Insert from HDMR.py
        Fx = np.array(Fx)
        Fx_df = pd.DataFrame(Fx[0:K + 1, 0:5]) #, columns=[Fx[0, i] for i in range(5)])
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.axis('tight')
        ax.axis('off')
        tbl = table(ax, Fx_df, loc='center', cellLoc='center', colWidths=[0.1]*len(Fx_df.columns))
        # Adjust font size and table appearance
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(6)
        tbl.scale(1, 1)  # Adjust table scaling if needed
        ax.set_title(f'Table 1: Emulator performance training data set', fontsize=fontsize_xylabel, pad=20)

        # Save the table to a PDF file
        pdf.savefig(fig)  # Save the current figure
        plt.close(fig)  # Close the plot to release memory
    
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        # Now plot SA_sig to table in figure
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        SA_sig = np.array(SA_sig)
        SA_sig_df = pd.DataFrame(SA_sig, columns=[SA_sig[0, i] for i in range(0, SA_sig.shape[1])])

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.axis('tight')
        ax.axis('off')
        tbl = table(ax, SA_sig_df, loc='center', cellLoc='center', colWidths=[0.08]*len(SA_sig_df.columns))
        # Adjust font size and table appearance
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(6)
        tbl.scale(1, 1)  # Adjust

        #ax.table(cellText=Fx_df[1:,:].values, colLabels=Fx_df.columns, loc='center', cellLoc='center', colLoc='center')
        ax.set_title(f'Table 2: Variance-based decomposition and sensitivity coefficients', fontsize=fontsize_xylabel, pad=20)
    
        pdf.savefig(fig)  # Save the current figure
        plt.close(fig)  # Close the plot to release memory

    # Open the final PDFs
    os.startfile(file_name)
#    os.startfile(table1_name)
#    os.startfile(table2_name)
 
 #   os.system(f'open {table1_name}')
 #   os.system(f'open {table2_name}')

    # Merge figures and tables into a single PDF file
 #   output_pdf = f'{file_name}'
 #   input_pdfs = [file_name, table1_name, table2_name]
    
    # Append PDFs
 #   with PdfPages(output_pdf) as pdf:
 #       for pdf_file in input_pdfs:
 #           with open(pdf_file, 'rb') as file:
 #               pdf.savefig(file)
    
    # Delete the temporary tables
 #   for table_file in [table1_name, table2_name]:
 #       os.remove(table_file)
    
    # Open the final PDF
#    os.system(f'open {file_name}')

#    print(' DONE')

# Latin Hypercube Sampling function
def LH_sampling(mn, mx, N):
    """
    Latin Hypercube Sampling.
    
    Args:
        mn: Lower bound vector
        mx: Upper bound vector
        N: Number of samples to generate

    Returns:
        N x d matrix of Latin Hypercube samples
    """
    d = len(mn)  # Number of parameters
    rng = np.array(mx) - np.array(mn)  # 1 x d vector with parameter ranges
    y =  np.random.rand(N, d)  # N x d matrix with uniform random labels
    # really important change below so that X stays in bound! as list is from 0 - N-1 rather than 1 to N
    id_matrix = 1 + np.argsort(np.random.rand(N, d), axis=0)  # Random sort (1:N without replacement)
    M = (id_matrix - y) / N  # Multiplier matrix (y introduces randomness)
    R = np.add(np.multiply(M, rng), mn)  # N x d matrix of stratified LH samples
    
    return R
