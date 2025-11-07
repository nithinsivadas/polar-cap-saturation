# Regression to the Mean of Extreme Geomagnetic Storms

This repository contains code that is associated with the paper titled "Regression to the Mean of Extreme Geomagnetic Storms"

# View the Code

Three of the main codes used to calculate all the results in the paper are available as .html files, so that the logic and results can be followed by anyone without having to use MATLAB or Python.

## Code used for the different results presented in the paper

1. [Code 1: Monte Carlo Simulation of the Solar Wind Error Model](Code1_Monte_Carlo_Simulation_of_Solar_Wind_Uncertainty.html). Here we solve the solar wind driver error model using a Monte Carlo simulation, compare and validate the model results with data, and also correct the data using 'regression calibration' to reveal that geomagnetic response to solar wind driving is linear.
2. [Code 2: Analytical derivation of non-linear bias from uncertainty in time](Code2_Analytical_derivation_of_time_uncertainty.html). Here we use the analytical formula for regression bias stemming from uncertainty in the time of the measurement of a stochastic process, derived in the paper. It is in the form of an integral, and hence in the code, we complete the numerical integration of the formula to show the non-linear regression bias it creates. Furthermore, we validate this analytical solution independently, by creating an error model of time uncertainty, and using a Monte Carlo simulation to solve it. The analytical solution and the Monte Carlo solution match very well, proving that the analytical derivation is correct. 
3. [Code 3: Calculating the magnitude uncertainty](Code3_Calculating_magnitude_uncertainty.html). Here we use data from simultaneous magnetosheath satellite measurements and spacecraft at L1, to calculate the variance of the L1 solar wind driver estimates for a given value of the same driver close to magnetopause subsolar point. The resulting variance is the magnitude uncertainty used in Code 1, as part of the error model. This magnitude uncertainty varies with the strength of the shocked solar wind driver, and hence is an heteroskedastic error, unlike what is normally assumed.
4. [Code 4: Regression to the mean](Code_4_Demonstrating_regression_bias.html). Here we demonstrate that regression to the mean is fundamental property of the relation between measurement and the true value corresponding to it. We also show how linear regression bias, and non-linear regression bias is created due to this property. Part of this code was first published by [Sivadas & Sibeck, 2022](https://doi.org/10.3389/fspas.2022.924976).

# Run the Code

If you wish to run the MATLAB and Python Code,

The Data files used by the code are stored in this [zenodo repository](https://doi.org/10.5281/zenodo.17546719), and can be downloaded from there. 

For running [Code 1: Monte Carlo Simulation of the Solar Wind Error Model](Code1_Monte_Carlo_Simulation_of_Solar_Wind_Uncertainty.html)
1. Download **Data.zip** from [Zenodo-record](https://doi.org/10.5281/zenodo.17546719) and unzip into this gitrepo on your local machine, and point MATLAB Code 1 variable `DataDir` to a corresponding folder in the gitrepo: **Data**. Example `DataDir = 'C:\Users\username\Documents\GitHub\polar-cap-saturation\Data\';`

For running [Code 3: Calculating the magnitude uncertainty](Code3_Calculating_magnitude_uncertainty.html)

2. Download *Datasets.zip* from [Zenodo-record](https://doi.org/10.5281/zenodo.17546719) and unzip into your local machine within the local copy of this git-repo. The python notebook `Code3_Calculating_magnitude_uncertainty.ipynb` uses this folder for its input data. So the code should be able to access this folder as `./Datasets/`. So ensure the data files are inside the your gitrepo like follows: `~\polar-cap-saturation\Datasets\`. 

[Code 2](Code2_Analytical_derivation_of_time_uncertainty.html) and [Code 4](Code_4_Demonstrating_regression_bias.html) do not require any data files to run. 

If there are any issues running the code, you can contact me at nithinsivadas90@gmail.com



