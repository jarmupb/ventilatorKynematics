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
# Synthesis opt 1
# =============================================================================
plt.close('all')

# =============================================================================
# Mechanism A - synthesis
# =============================================================================
# Rocker conditions
delta4 = 25. * np.pi / 180.
q4mA = -90. * np.pi / 180.
q41A = q4mA - delta4 / 2.
q42A = q4mA + delta4 / 2.
# Input conditions
delta2 = 45. * np.pi / 180.
q2mA = -90. * np.pi / 180.
q21A = q2mA - delta2 / 2.
q22A = q2mA + delta2 / 2.
# Chasis
r1Ax = 167.5 + 40.
r1Ay = 19. + 40.
# Input link
r2 = 25.

# Initialization for solver
r30 = 165.
r40 = 40.
q310 = 0.
q320 = 0.
x0 = np.array([r30, r40, q310, q320])

# Solver
# args = (r1x, r1y, r2, q21, q22, q41, q42)
sol = fsolve(rrr_res, x0, args = (r1Ax, r1Ay, r2, q21A, q22A, q41A, q42A))
r3A = sol[0]
r4A = sol[1]
print(r3A, r4A)

mecA = vt.RRR(r2, r3A, r4A, r1Ax, r1Ay, inv = 1.)

# =============================================================================
# Mechanism B - synthesis
# =============================================================================
# Rocker conditions
delta4 = 25. * np.pi / 180.
q4mB = 183.5 * np.pi / 180.
q41B = q4mB - delta4 / 2.
q42B = q4mB + delta4 / 2.
# Input conditions
delta2 = 45. * np.pi / 180.
q2mB = 180. * np.pi / 180.
q21B = q2mB - delta2 / 2.
q22B = q2mB + delta2 / 2.

# delta4 = 25. * np.pi / 180.
# q41B = 167.5 * np.pi / 180.
# q42B = 167.5 * np.pi / 180. + delta4
# # Input conditions
# q21B = 157.5 * np.pi / 180.
# q22B = 202.5 * np.pi / 180.
# Chasis
r1Bx = 11. + 40.
r1By = 147.4 + 40.
# Input link
r2 = 25.

# Initialization for solver
r30 = 154.
r40 = 45.
q310 = 1.5
q320 = 1.5
x0 = np.array([r30, r40, q310, q320])

# Solver
# args = (r1x, r1y, r2, q21, q22, q41, q42)
sol = fsolve(rrr_res, x0, args = (r1Bx, r1By, r2, q21B, q22B, q41B, q42B))
r3B = sol[0]
r4B = sol[1]
print(r3B, r4B)
print(sol)

mecB = vt.RRR(r2, r3B, r4B, r1Bx, r1By, inv = -1.)
# mecB.coeff_list()

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
plt.legend(['A', 'B'])
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$VM$')


# =============================================================================
# Synthesis - opt 2
# =============================================================================
# =============================================================================
# Mechanism A - synthesis
# =============================================================================
# Rocker conditions
delta4 = 25. * np.pi / 180.
q4mA = -107. * np.pi / 180.
q41A = q4mA - delta4 / 2.
q42A = q4mA + delta4 / 2.
# Input conditions
delta2 = 45. * np.pi / 180.
q2mA = -67.5 * np.pi / 180.
q21A = q2mA - delta2 / 2.
q22A = q2mA + delta2 / 2.
# Chasis
r1Ax = 167.5 + 40.
r1Ay = 19. + 40.
print('Angle:', np.arctan(r1Ay / r1Ax) * 180. / np.pi)
# Input link
r2 = 25.

# Initialization for solver
r30 = 165.
r40 = 40.
q310 = 0.
q320 = 0.
x0 = np.array([r30, r40, q310, q320])

# Solver
# args = (r1x, r1y, r2, q21, q22, q41, q42)
sol = fsolve(rrr_res, x0, args = (r1Ax, r1Ay, r2, q21A, q22A, q41A, q42A))
r3A = sol[0]
r4A = sol[1]
print(r3A, r4A)

mecA = vt.RRR(r2, r3A, r4A, r1Ax, r1Ay, inv = 1.)

# =============================================================================
# Mechanism B - synthesis
# =============================================================================
# Rocker conditions
delta4 = 25. * np.pi / 180.
q4mB = 159. * np.pi / 180.
q41B = q4mB - delta4 / 2.
q42B = q4mB + delta4 / 2.
# Input conditions
delta2 = 45. * np.pi / 180.
q2mB = 202.5 * np.pi / 180.
q21B = q2mB - delta2 / 2.
q22B = q2mB + delta2 / 2.

# delta4 = 25. * np.pi / 180.
# q41B = 167.5 * np.pi / 180.
# q42B = 167.5 * np.pi / 180. + delta4
# # Input conditions
# q21B = 157.5 * np.pi / 180.
# q22B = 202.5 * np.pi / 180.
# Chasis
r1Bx = 11. + 40.
r1By = 147.4 + 40
print('Angle:', np.arctan(r1Bx / r1By) * 180. / np.pi)
# Input link
r2 = 25.

# Initialization for solver
r30 = 154.
r40 = 45.
q310 = 1.5
q320 = 1.5
x0 = np.array([r30, r40, q310, q320])

# Solver
# args = (r1x, r1y, r2, q21, q22, q41, q42)
sol = fsolve(rrr_res, x0, args = (r1Bx, r1By, r2, q21B, q22B, q41B, q42B))
r3B = sol[0]
r4B = sol[1]
print(r3B, r4B)
print(sol)

mecB = vt.RRR(r2, r3B, r4B, r1Bx, r1By, inv = -1.)
# mecB.coeff_list()


# Input angle config
delta_ang = 1.
q = np.arange(45., 90. + delta_ang, delta_ang) * np.pi / 180. 
    
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
plt.legend(['A', 'B'])
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$VM$')