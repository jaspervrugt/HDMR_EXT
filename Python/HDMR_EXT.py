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
#  SYNOPSIS                                                               #
#   [S,Ss,Fx,Em,Xy,RT] = HDMR_EXT(X,y);                                   #
#   [S,Ss,Fx,Em,Xy,RT] = HDMR_EXT(X,y,options);                           #
#  where                                                                  #
#   X           [input] Nxd matrix: N vectors of d parameters             #
#   y           [input] Nx1 vector: single model output each row X        #
#   options     [input] (optional) structure: HDMR variables              #
#    .graphics   integer [0,1]: graphical output?             (def: 1)    #
#    .basis      integer [1-3]: type of the orthonormal basis (def: 1)    #
#                         1: orthonormal polynomial                       #
#                         2: legendre polynomial                          #
#                         3: hermite polynomial                           #
#    .maxorder   integer [2-3]: max order emulator            (def: 3)    #
#    .m          integer [1-12]: polynomial degree            (def: 3)    #
#    .K          integer [1-500] # bootstrap iterations       (def: 100)  #
#    .R          integer [100-N/2] # bootstrap samples        (def: N/2)  #
#    .alfa       real [0.0-.1]: significance level F-test     (def: 0.01) #
#    .method     integer [1,2]: 1=forw. sel.; 2=backw. elim.  (def: 1)    #
#   def_options = struct('graphics',1,'basis','1','maxorder',3,...        #
#                     'm',3,'K',1,'R',N/2,'alfa',0.01,'method',1);        #
#   S           [outpt] Cell array: structural, correlative & total       #
#                       sensitivity each component function -> 1st, 2nd   #
#                       and 3rd order effects                             #
#                       (= Table 2 screen and in "HDMR_EXT_results.txt")  #
#   Ss          [outpt] Cell array: As "S" but lists only significant     #
#                       component functions determined via model          #
#                       selection using a F-test                          #
#                       (= Table 3 screen and in "HDMR_EXT_results.txt")  #
#   Fx          [outpt] Cell array: Tabulates emulator properties and     #
#                       performance on training data set for each of the  #
#                       K bootstrap trials                                #
#                       (= Table 1 screen and in "HDMR_EXT_results.txt" ) #
#   Em          [outpt] Structure array: Fields input/output K emulators  #
#    .EB         Nxk matrix: extended bases evaluated at X                #
#    .C          kx1 matrix: all order coefficients (via D-Morph)         #
#    .m          scalar: degree of orthonormal polynomials                #
#    .Y_e        RxK matrix: emulator predictions of K bootstraps         #
#    .RMSE       1xK vector: RMSE of emulator residuals of K bootstraps   #
#    .m1         1x1 scalar: # 1st order terms ( = m )                    #
#    .m2         1x1 scalar: # 2nd order terms ( = 2*m+m^2 )              #
#    .m3         1x1 scalar: # 3rd order terms ( = 3*m+3*m^2+m^3 )        #
#    .f0         1xK vector: mean y of each of K bootstrap trials         #
#    .n          scalar: ( = n1+n2+n3 ) total # of component functions    #
#    .n1         scalar: ( = d) with total number of 1st order terms      #
#    .n2         scalar: total number of 2nd order component functions    #
#    .n3         scalar: total number of 3rd order component functions    #
#    .k          scalar: total number of terms                            #
#    .k1         scalar: total number of 1st order terms                  #
#    .k2         scalar: total number of 2nd order terms                  #
#    .k3         scalar: total number of 3rd order terms                  #
#    .maxorder   scalar: maximum order of HDMR emulator                   #
#    .nterms     1xK vector: # significant terms of each of K emulators   #
#    .p0         1xK vector: # parameters of each of the K emulators      #
#    .RT         1xK vector: CPU time (in sec) to construct each emulator #
#    .select     nxK matrix: significant/insignifant terms each K trials  #
#                if Em.select(i,1) = 0 -> term i insignificant in trial 1 #
#                if Em.select(i,4) = 1 -> term i significant in trial 4   #
#   Xy          Structure: X/y samples and bootstrap information          #
#    .R          scalar: number of random samples each bootstrap trial    #
#    .X_n        Nxd matrix: normalized parameter vectors                 #
#                X_n(i,1:d) = (X(i,1:d)-X_min)./(X_max-X_min); i = 1,..,N #
#    .minX       1xd vector: min values each X column ( = input )         #
#    .maxX       1xd vector: max values each X column ( = input )         #
#    .y          Nx1 vector: y values supplied by user                    #
#    .id         RxK matrix: index R samples of X each bootstrap trial    #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  BUILT-IN CASE STUDIES                                                  #
#   Example 1   Multivariate normal benchmark study                       #
#   Example 2   Multivariate normal benchmark study with Σ ≠ identity     #
#   Example 3   Multivariate normal benchmark study with Σ ≠ identity     #
#   Example 4   Ishigami function                                         #
#   Example 5   Ishigami function with Σ ≠ identity                       #
#   Example 6   Function with Σ ≠ identity                                #
#   Example 7   Sobol g function                                          #
#   Example 8   Multivariate normal with Σ ≠ identity                     #
#   Example 9   Example 1 of Chastaing et al. (2012)                      #
#   Example 10  Example 2 of Chastaing et al. (2012)                      #
#   Example 21  Soil temperature modeling                                 #
#   Example 22  Rainfall runoff modeling: hmodel                          #
#   Example 23  Rainfall runoff modeling: SAC-SMA model                   #
#   Example 24  Dynamic foodweb: One-predator-one-prey model              #
#   Example 25  Dynamic foodweb: Two-predators-two-preys model            #
#   Example 26  Dynamic foodweb: Two-predators-two-preys model real data  #
#   Example 27  Simple multivariate function                              #
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
import time
from typing import Dict, Any
from scipy.stats import f
from HDMR_EXT_functions import *

def HDMR_EXT(X: np.ndarray, y: np.ndarray, options: Dict[str, Any] = None) -> tuple:
    """
    Extended Bases High-Dimensional Model Representation (HDMR)

    Args:
    - X: Nxd matrix of input parameters.
    - y: Nx1 vector of output values.
    - options: A dictionary containing optional parameters.
        - graphics (int): 0 or 1 (default: 1)
        - basis (int): Polynomial type (1: Classical, 2: Legendre, 3: Hermite, default: 1)
        - maxorder (int): Maximum order emulator (default: 3)
        - m (int): Polynomial degree (default: 3)
        - K (int): Number of bootstrap iterations (default: 100)
        - R (int): Number of bootstrap samples (default: N/2)
        - alfa (float): Confidence interval F-test (default: 0.99)
        - method (int): 1=Forward selection, 2=Backward elimination (default: 1)

    Returns:
    - S: Sensitivity indices (structural, correlative, and total).
    - Ss: Significant component functions.
    - Fx: Emulator properties and performance.
    - Em: Emulator details.
    - Xy: Information on X and y samples and bootstrap.
    """

    if options is None:
        options = {}

    # Set up default options if not provided
    options.setdefault('graphics', 1)
    options.setdefault('basis', 1)
    options.setdefault('maxorder', 3)
    options.setdefault('m', 2)
    options.setdefault('K', 100)
    options.setdefault('R', X.shape[0] // 2)
    options.setdefault('alfa', 0.01)
    options.setdefault('method', 1)

    # Ensure y is a column vector
    y_in = y.flatten()

    if len(y) < 2:
        raise ValueError("HDMR_EXT ERROR: Insufficient number of input arguments. Use 'help install_HDMR_EXT'")
    else:
        # Close all hidden (Not implemented in Python)
        pass

    y_in = np.array(y).flatten()  # Ensure y is a column vector and ct is zero
    N, d, graphics, basis, maxorder, m, K, R, method, alfa = HDMR_EXT_setup(X, y_in, options)

    # Initialize variables
    Xy, Em, S, n_ns, ct = HDMR_EXT_initialize(X, y_in, N, d, K, R, m, maxorder, basis)

    # Main computation loop (Bootstrap for confidence intervals)
    for k in range(K):
        start_time = time.time()
        y_bootstrap = y_in[Xy['id'][:R, k]]  # Extract the "right" Y values
        S['V_y'][k] = np.var(y_bootstrap)  # Variance of the model output
        Em['f0'][k] = np.mean(y_bootstrap)  # Mean of the output
        Y_res = y_bootstrap - Em['f0'][k]  # Compute residual

        # Construct cost matrix, B
        B, C = HDMR_EXT_construct_B(Em['EB'][Xy['id'][:R, k]], R, d, m, Em['m2'], Em['m3'], Em['k1'], Em['k1']+Em['k2'], Em['k'], maxorder)

        # Perform D-Morph regression
        Em['C'], n_ns[k] = HDMR_EXT_dmorph(Em['EB'][Xy['id'][:R, k]], B, C, Y_res, R, Em['k1'], Em['k'])

        # Compute component functions
        Y_em = HDMR_EXT_comp_func(Em['EB'][Xy['id'][:R, k]], Em['C'], R, d, m, Em['k1'], Em['k2'], Em['k3'], Em['k'], Em['n1'], Em['n2'], Em['n3'], maxorder)

        # Perform F-test for selecting significant terms
        Em['select'][:Em['n'], k], Em['nterms'][k], Em['p0'][k], Em['RMSE'][k] = HDMR_EXT_F_test(y_bootstrap, Em['f0'][k], Y_em, R, alfa, Em['m1'], Em['m2'], Em['m3'], Em['n1'], Em['n2'], Em['n3'], Em['n'], method)
        
        # Store emulated output
        Em['Y_e'][:, k] = Em['f0'][k] + np.sum(Y_em, axis=1)

        # Compute sensitivity indices (ANCOVA)
        S['S'][:Em['n'], k], S['Sa'][:Em['n'], k], S['Sb'][:Em['n'], k], S['V_em'][:Em['n'], k] = ANCOVA(y_bootstrap, Y_em, S['V_y'][k], R, Em['n'])

        # Record CPU time for this bootstrap iteration
        Em['RT'][k] = time.time() - start_time
        
        # Print progress
        print(f"HDMR_EXT CALCULATING: FINISHED TRIAL {k+1} OF {K} AND {100 * (k+1) / K:.2f}% DONE", end='\r')

    # Postprocessing: Generate output tables and graphical output
    S, Ss, Fx = HDMR_EXT_end(S, Em['nterms'], Em['p0'], Em['RMSE'], Em['select'], K, Em['beta'], Em['gamma'], Em['n1'], Em['n2'], Em['n3'], Em['n'], n_ns)

    if graphics == 1:
        HDMR_EXT_plot(Ss, Fx, y_in, Em['Y_e'], Em['select'], Em['p0'], n_ns, Xy['id'], R, K)

    return S, Ss, Fx, Em, Xy