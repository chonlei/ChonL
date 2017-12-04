#!/usr/bin/env python2
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

tol_err = np.loadtxt("tol_R_Chaste.txt", delimiter=',')
tol = np.arange(6,6+10)

plt.plot(tol, np.log(tol_err[:,0]), 'o-', label='CellML/Chaste')
plt.plot(tol, np.log(tol_err[:,1]), 'o--', label='C/R (CiPA)')
plt.xlabel('|log(tolerance)|')
plt.ylabel('log(||state variable difference||)')
plt.legend()
plt.savefig('solver_tol_test.png', bbox_inch='tight')
