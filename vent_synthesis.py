# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 15:27:28 2020

@author: UPB
"""

# =============================================================================
# Modules
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import vent_analysis as vt
from scipy.optimize import fsolve

# =============================================================================
# Model for synthesis
# =============================================================================
def rrr_res(x, r1x, r1y, r2, q21, q22, q41, q42):
    r3 = x[0]
    r4 = x[1]
    q31 = x[2]
    q32 = x[3]
    
    res = np.zeros(4)
    res[0] = r2 * np.cos(q21) + r3 * np.cos(q31) - r1x - r4 * np.cos(q41)
    res[1] = r2 * np.sin(q21) + r3 * np.sin(q31) - r1y - r4 * np.sin(q41)
    res[2] = r2 * np.cos(q22) + r3 * np.cos(q32) - r1x - r4 * np.cos(q42)
    res[3] = r2 * np.sin(q22) + r3 * np.sin(q32) - r1y - r4 * np.sin(q42)
    
    return res

# =============================================================================
# Synthesis
# =============================================================================
plt.close('all')

# =============================================================================
# Mechanism A - synthesis
# =============================================================================
# Rocker conditions
delta4 = 25. * np.pi / 180.
q41 = -120. * np.pi / 180.
q42 = q41 + delta4
# Input conditions
q21A = -135. * np.pi / 180.
q22A = -90. * np.pi / 180.
# Chasis
r1Ax = 167.5
r1Ay = 19.
# Input link
r2 = 25.

# Initialization for solver
r30 = 165.
r40 = 40.
q310 = 0.75
q320 = 1.41
x0 = np.array([r30, r40, q310, q320])

# Solver
# args = (r1x, r1y, r2, q21, q22, q41, q42)
sol = fsolve(rrr_res, x0, args = (r1Ax, r1Ay, r2, q21A, q22A, q41, q42))
r3 = sol[0]
r4 = sol[1]
print(r3, r4)

mecA = vt.RRR(r2, r3, r4, r1Ax, r1Ay, inv = 1.)

# =============================================================================
# Mechanism B - synthesis
# =============================================================================
# Rocker conditions
delta4 = 25. * np.pi / 180.
q41 = -214. * np.pi / 180.
q42 = q41 + delta4
# Input conditions
q21B = 135. * np.pi / 180.
q22B = 180. * np.pi / 180.
# Chasis
r1Bx = 11.
r1By = 147.4
# Input link
r2 = 25.

# Initialization for solver
r30 = 154.
r40 = 40.
q310 = 0.75
q320 = 1.41
x0 = np.array([r30, r40, q310, q320])

# Solver
# args = (r1x, r1y, r2, q21, q22, q41, q42)
sol = fsolve(rrr_res, x0, args = (r1Bx, r1By, r2, q21B, q22B, q41, q42))
r3 = sol[0]
r4 = sol[1]
print(r3, r4)

mecB = vt.RRR(r2, r3, r4, r1Bx, r1By, inv = -1.)

# =============================================================================
# Drawing of extreme positions
# =============================================================================
plt.figure()
mecA.update_posi(q21A)
mecA.draw_posi()
mecA.update_posi(q22A)
mecA.draw_posi()

mecB.update_posi(q21B)
mecB.draw_posi()
mecB.update_posi(q22B)
mecB.draw_posi()

plt.axis('equal')

