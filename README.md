# High Dimensional Model Representation with extended bases (HDMR_EXT): (Co)variance-based Sensitivity Analysis with Uncorrelated and Correlated (Dependent) Input Variables: MATLAB and Python Toolboxes

## Description

High-Dimensional Model Representation with extended bases (HDMR_EXT) uses three different type of orthonormal polynomials (Classical, Legendre, and Hermite) for variance-based global sensitivity analysis (GSA) with correlated and uncorrelated inputs. This function uses N x d matrix of N different vector d-vectors of model inputs (factors/parameters) and a N x 1 vector of corresponding model outputs and returns to the user each factor's first, second, and third order sensitivity coefficient (separated in total, structural and correlative contributions), an estimate of their 95% confidence intervals (from bootstrap method) and the coefficients of the extended bases functions that govern output, Y emulator (determined by an F-test of the error residuals of the HDMR model with/without a given first, second and/or third order components). These coefficients define an emulator that can be used to predict the output, Y, of the original model for any d-vector of model inputs. For uncorrelated model inputs (columns of X are independent), the HDMR sensitivity indices reduce to a single index (= structural contribution) consistent with their values derived from commonly used variance-based GSA methods.

## Getting Started

### Installing: MATLAB

* Download and unzip the zip file 'MATLAB_code_HDMR_EXT_V2.0.zip' in a directory 'HDMR_EXT'
* Add the toolbox to your MATLAB search path by running the script 'install_HDMR_EXT.m' available in the root directory
* You are ready to run the examples.

### Executing program

* After intalling, you can simply direct to each example folder and execute the local 'example_X.m' file.
* Please make sure you read carefully the instructions (i.e., green comments) in 'install_HDMR_EXT.m'  

### Installing: Python

* Download and unzip the zip file 'Python_code_HDMR_EXT_V2.0.zip' to a directory called 'HDMR_EXT'.

### Executing program

* Go to Command Prompt and directory of example_X in the root of 'HDMR_EXT'
* Now you can execute this example by typing 'python example_X.py'
* Instructions can be found in the file 'HDMR_EXT.py' 
  
## Authors

* Vrugt, Jasper A. (jasper@uci.edu) 

## Version History

* 1.0
    * Initial Release
* 2.0
    * Python implementation
    * New built-in case studies

## Acknowledgments
The MATLAB toolbox is based on the published works of Dr. Genyuan Li from Princeton University. We greatly appreciate his diligent responses to our questions. 

