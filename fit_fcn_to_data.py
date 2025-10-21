import numpy as np 
from numpy.linalg import inv

# print(f'fcn: {fit_terms[2](np.exp(2.))}')
def fit_fcn_to_data(fcn_terms, points_x, points_y, points_y_err):
    '''given an array of function terms, and some numpy arrays, fit the fcn.'''

    n_dof = len(fcn_terms)

    Xi = np.empty([n_dof, len(points_x)])

    # print(f"shape(Xi): {np.shape(Xi)}\nshape(Xi[0]): {np.shape(Xi[0])}")

    for i in range(0, n_dof): 
        Xi[i] = fcn_terms[i](points_x)

    #now, create the A-matrix
    A = np.zeros([n_dof, n_dof])
    B = np.zeros(n_dof)

    points_sigma = points_y_err ** 2

    for i in range(0, n_dof):

        B[i] = np.sum(points_y * Xi[i] / points_sigma )

        for j in range(0, n_dof):
            A[i,j] = np.sum( Xi[i] * Xi[j] / points_sigma )

    coeffs = inv(A) @ B

    # print(f'shape coeffs: {np.shape(coeffs)}')
    # compute chi2
    
    y_model = np.zeros(len(points_x))
    for i in range(0, len(y_model)):
        for j in range(0, n_dof):
            y_model[i] += Xi[j,i] * coeffs[j]

    chi2 = np.sum( (( points_y - y_model ) ** 2) / points_sigma )

    return coeffs, chi2 

