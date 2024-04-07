# FMCC
Code used in "Flexible modeling of a control function approach subject to different types of dependent censoring"

This repository contains all R code used in simulations and data applications of the above mentioned paper (https://arxiv.org/abs/2403.11860).

The scripts __Functions_ad.R__, __Functions_cmprsk.R__ and __Goodness-of-fit-test_functions.R__ contain the various low level and high level functions used to implement the methodology and to run the simulations/data applications. __Functions_ad.R__ does this for the main model in the paper, while __Functions_cmprsk.R__ focusses on the extension of the model to the competing risks framework. __Goodness-of-fit-test_functions.R__ contains all of the code related to the goodness-of-fit test.

These functions are then called in the scripts performing the simulations (__Code_main_ad.R__ and __Code_main_cmprsk.R__) and the data applications (__Data_application_ad.R__ and __Data_application_cmprsk.R__). The folder __Goodness of Fit__ contains some more scripts, which are of lesser importance.

The cleaned JTPA dataset was obtained from the paper by _Crommen et al. (2024)_ and can be found, alongside some further documentation, on their GitHub page (https://github.com/GillesCrommen/DCC).
