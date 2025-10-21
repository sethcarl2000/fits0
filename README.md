# Fiting Exercises

## LSQFit.py (Python) 
Made a general function ```fit_fcn_to_data.py``` which takes an array of functions, each of which is its own term in the fit (i.e., f_1(x), f_2(x), f_3(x), ...). and also takes arrays of x-values, y-values and y-errors. 
To see the results of this fit, see ```fits-100-pts-py.png``` 

## LSQFit.cpp (C++) 
Same scheme as above, see ```fit_fcn_to_data()``` fcn. It uses the ```<functional>``` header to take a vector of functions as input (each element of the vector represents a term in the fit), and returns a 'FitResult_t' struct, which gives you the coefficients and chi2 of the fit. 

This is demonstrated in the plots ```fits-*-pts.png``` files, where * is the number of fit-points per trial, 4-128. each point value has 10^5 trails. 

## plots_and_comments.pdf
Here are the plots and comments. 
