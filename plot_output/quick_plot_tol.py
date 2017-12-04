#!/usr/bin/env python2
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

tol_test = np.loadtxt("solver_tol_test.txt")
tol = tol_test[:,0]
log_tol = np.abs(np.log(tol))
tol_err = tol_test[:,1:]
del(tol_test)

plt.plot(log_tol, np.log(tol_err[:,0]), 'o-', label='CellML/Chaste')
plt.plot(log_tol, np.log(tol_err[:,1]), 'o--', label='C/R (CiPA)')
plt.xlabel('|log(tolerance)|')
plt.ylabel('log(||state variable difference||)')
plt.legend()
plt.savefig('solver_tol_test.png', bbox_inch='tight')
