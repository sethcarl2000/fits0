import numpy as np 
from typing import Callable
from numpy.linalg import inv


def fit_fcn_to_data(fcn : Callable[[float],float], data):
    '''given a function, and a numpy 2d array, fit the fcn.'''
