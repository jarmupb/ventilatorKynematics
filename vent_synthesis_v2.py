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
# Model for simultaneous mechanism synthesis
# =============================================================================
def rrr_res(x, r1Ax, r1Ay, r2, q21A, q22A, q21B, q22B, q4mB, delta4):
    r3A = x[0]
    r3B = x[1]
    r4 = x[2]
    q31A = x[3]
    q32A = x[4]
    q31B = x[5]
    q32B = x[6]
    q4mA = x[7]
    
    res = np.zeros(8)
    
    # Mechanism A
    res[0] = r2 * np.cos(q21A) + r3A * np.cos(q31A) \
        - r1Ax - r4 * np.cos(q4mA - delta4 / 2.)
    res[1] = r2 * np.sin(q21A) + r3A * np.sin(q31A) \
        - r1Ay - r4 * np.sin(q4mA - delta4 / 2.)
    res[2] = r2 * np.cos(q22A) + r3A * np.cos(q32A) \
        - r1Ax - r4 * np.cos(q4mA + delta4 / 2.)
    res[3] = r2 * np.sin(q22A) + r3A * np.sin(q32A) \
        - r1Ay - r4 * np.sin(q4mA + delta4 / 2.)
    
    # Mechanism B
    res[4] = r2 * np.cos(q21B) + r3B * np.cos(q31B) \
        - r1Bx - r4 * np.cos(q4mB - delta4 / 2.)
    res[5] = r2 * np.sin(q21B) + r3B * np.sin(q31B) \
        - r1By - r4 * np.sin(q4mB - delta4 / 2.)
    res[6] = r2 * np.cos(q22B) + r3B * np.cos(q32B) \
        - r1Bx - r4 * np.cos(q4mB + delta4 / 2.)
    res[7] = r2 * np.sin(q22B) + r3B * np.sin(q32B) \
        - r1By - r4 * np.sin(q4mB + delta4 / 2.)
    
    return res

# =============================================================================
# Synthesis opt 1
# =============================================================================
plt.close('all')
# =============================================================================
# General parameters
# =============================================================================
# Input link
r2 = 25.
delta2 = 45. * np.pi / 180.
# Rocker
delta4 = 27. * np.pi / 180.
# =============================================================================
# Mechanism A - parameters
# =============================================================================
# Input conditions
q2mA = -90. * np.pi / 180.
q21A = q2mA - delta2 / 2.
q22A = q2mA + delta2 / 2.
# Rocker conditions
q4mA = -75. * np.pi / 180.
q41A = q4mA - delta4 / 2.
q42A = q4mA + delta4 / 2.
# Chasis
r1Ax = 167.5 + 40.
r1Ay = 19. + 40.
print('Angle_A:', np.arctan(r1Ay / r1Ax ) * 180. / np.pi)

# =============================================================================
# Mechanism B - parameters
# =============================================================================
# Input conditions
q2mB = 180. * np.pi / 180.
q21B = q2mB - delta2 / 2.
q22B = q2mB + delta2 / 2.
# Rocker conditions
q4mB = 165. * np.pi / 180.
q41B = q4mB - delta4 / 2.
q42B = q4mB + delta4 / 2.
# Chasis
r1Bx = 11. + 40.
r1By = 147.4 + 40.
print('Angle_B:', np.arctan(r1Bx / r1By) * 180. / np.pi)

# Initialization for solver
r3A0 = 165.
r3B0 = 154.
r40 = 40.
q31A0 = 0.
q32A0 = 0.
q31B0 = 1.5
q32B0 = 1.5
x0 = np.array([r3A0, r3B0, r40, q31A0, q32A0, q31B0, q32B0, q4mA])

# Solver
# rrr_res(x, r1Ax, r1Ay, r2, q21A, q22A, q21B, q22B, q4Bm, delta4)
sol = fsolve(rrr_res, x0,
             args = (r1Ax, r1Ay, r2, q21A, q22A, q21B, q22B, q4mB, delta4))
r3A = sol[0]
r3B = sol[1]
r4 = sol[2]
q4mA = sol[7]
print(r3A, r3B, r4, q4mA * 180. / np.pi)

mecA = vt.RRR(r2, r3A, r4, r1Ax, r1Ay, inv = 1.)
mecB = vt.RRR(r2, r3B, r4, r1Bx, r1By, inv = -1.)

# Input angle config
delta_ang = 1.
q = np.arange(22.5, 67.5 + delta_ang, delta_ang) * np.pi / 180. 
    
# Input link position
q2A = q - 3. * np.pi / 4.
q2B = q + 3. * np.pi / 4.

# Compute a list of positions
mecA.posi_list(q2A)
mecB.posi_list(q2B)

# Compute a list of vel coefficients    
mecA.coeff_list(q2A)
mecB.coeff_list(q2B)

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
plt.grid()
plt.xlabel('$x$ (mm) - horz')
plt.ylabel('$y$ (mm) - vert')

# Rocker angle
plt.figure()
plt.plot(q * 180. / np.pi, mecA.q4_list * 180. / np.pi)
plt.plot(q * 180. / np.pi, mecB.q4_list * 180. / np.pi)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$\\theta_4$ (deg)')
plt.grid()
plt.legend(['A', 'B'])

# Rocker angle delta
plt.figure()
plt.plot(q * 180. / np.pi, (mecA.q4_list - mecA.q4_list[0]) * 180. / np.pi)
plt.plot(q * 180. / np.pi, (mecB.q4_list - mecB.q4_list[0]) * 180. / np.pi)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$\Delta \\theta_4$ (deg)')
plt.grid()
plt.legend(['A', 'B'])

# Velocity coeffs
plt.figure()
plt.plot(q * 180. / np.pi, mecA.k4_list)
plt.plot(q * 180. / np.pi, mecB.k4_list)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$K_4 = 1 / VM$')
plt.legend(['A', 'B'])

# VM
plt.figure()
plt.plot(q * 180. / np.pi, 1. / mecA.k4_list, '.')
plt.plot(q * 180. / np.pi, 1. / mecB.k4_list, '+')
plt.plot(q * 180. / np.pi, 0.5 * (1. / mecA.k4_list + 1. / mecB.k4_list))
plt.legend(['A', 'B'])
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$VM$')
