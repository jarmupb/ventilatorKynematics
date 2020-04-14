# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 18:37:14 2020

@author: UPB
"""
# =============================================================================
# Modules
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import vent_analysis as vt

# =============================================================================
# Object creation
# =============================================================================
r1Ax = 207.5
r1Ay = 59.
r1Bx = 51.
r1By = 187.4
mecA0 = vt.RRR(r2 = 25., r3 = 224., r4 = 40., r1x = r1Ax, r1y = r1Ay, inv = 1.)
mecB0 = vt.RRR(r2 = 25., r3 = 200., r4 = 40., r1x = r1Bx, r1y = r1By, inv = -1.)

# Cam parameters
r4c = 150. # cam length (drawing only)
# angle between rocker and cam
delta4A = 176. * np.pi / 180.
delta4B = 117. * np.pi / 180.

# =============================================================================
# Position analysis - compute a set of angles
# =============================================================================
# Input angle config
delta_ang = 1.
q = np.arange(22.5, 67.5 + delta_ang, delta_ang) * np.pi / 180.
    
# Input link position
q2A = q - 3. * np.pi / 4.
q2B = q + 3. * np.pi / 4.

# Compute a list of positions
mecA0.posi_list(q2A)
mecB0.posi_list(q2B)

# Compute a list of vel coefficients    
mecA0.coeff_list(q2A)
mecB0.coeff_list(q2B)

# =============================================================================
# Plot
# =============================================================================
plt.close('all')

# Rocker angle
plt.figure()
plt.plot(q * 180. / np.pi, mecA0.q4_list * 180. / np.pi)
plt.plot(q * 180. / np.pi, mecB0.q4_list * 180. / np.pi)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$\\theta_4$ (deg)')
plt.grid()
plt.legend(['A', 'B'])

# Rocker angle delta
plt.figure()
plt.plot((q - q[0]) * 180. / np.pi, (mecA0.q4_list - mecA0.q4_list[0]) * 180. / np.pi)
plt.plot((q - q[0]) * 180. / np.pi, (mecB0.q4_list - mecB0.q4_list[0]) * 180. / np.pi)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$\Delta \\theta_4$ (deg)')
plt.grid()
plt.legend(['A', 'B'])

# Velocity coeffs
plt.figure()
plt.plot(q * 180. / np.pi, mecA0.k4_list)
plt.plot(q * 180. / np.pi, mecB0.k4_list)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$K_4 = 1 / VM$')
plt.legend(['A', 'B'])
plt.grid()

# VM
plt.figure()
plt.plot((q - q[0]) * 180. / np.pi, 1. / mecA0.k4_list)
plt.plot((q - q[0]) * 180. / np.pi, 1. / mecB0.k4_list)
plt.legend(['A', 'B'])
plt.xlabel('$\Delta \\theta_2$ (deg)')
plt.ylabel('$VM$')
plt.grid()

# Positions drawing
plt.figure()
for k in range(0, q.size, 10):
    mecA0.update_posi(q2A[k])
    mecA0.draw_posi()
    mecB0.update_posi(q2B[k])
    mecB0.draw_posi()
    
    R4A_c = r4c * vt.vunitary(mecA0.q4_list[k] + delta4A)
    R4B_c = r4c * vt.vunitary(mecB0.q4_list[k] + delta4B)
    
    # Cam follower
    plt.plot([mecA0.r1x, mecA0.r1x + R4A_c[0]],
             [mecA0.r1y, mecA0.r1y + R4A_c[1]], 'C2.--')
    plt.plot([mecB0.r1x, mecB0.r1x + R4B_c[0]],
             [mecB0.r1y, mecB0.r1y + R4B_c[1]], 'C2.--')

plt.text(0,0, '(0,0)')
plt.text(r1Ax, r1Ay, '(' + str(r1Ax) + ', ' + str(r1Ay) + ')')
plt.text(r1Bx, r1By, '(' + str(r1Bx) + ', ' + str(r1By) + ')')
plt.text(mecA0.R2[0] + 0.5 * mecA0.R3[0],
         mecA0.R2[1] + 0.5 * mecA0.R3[1], str(mecA0.r3))
plt.text(mecB0.R2[0] + 0.5 * mecB0.R3[0],
         mecB0.R2[1] + 0.5 * mecB0.R3[1], str(mecB0.r3))
plt.minorticks_on()
plt.grid(b = True, which = 'both', axis='both')
plt.axis('equal')
plt.xlabel('$x$ (mm) - hor')
plt.ylabel('$y$ (mm) - vert')


# =============================================================================
# Position analysis - compute a set of angles 60 deg version
# =============================================================================
# Input angle config
delta_ang = 1.
q = np.arange(22.5, 67.5 + 15. + delta_ang, delta_ang) * np.pi / 180.
    
# Input link position
q2A = q - 3. * np.pi / 4.
q2B = q + 3. * np.pi / 4.

# Compute a list of positions
mecA0.posi_list(q2A)
mecB0.posi_list(q2B)

# Compute a list of vel coefficients    
mecA0.coeff_list(q2A)
mecB0.coeff_list(q2B)

# =============================================================================
# Plot
# =============================================================================
plt.close('all')

# Rocker angle
plt.figure()
plt.plot(q * 180. / np.pi, mecA0.q4_list * 180. / np.pi)
plt.plot(q * 180. / np.pi, mecB0.q4_list * 180. / np.pi)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$\\theta_4$ (deg)')
plt.grid()
plt.legend(['A', 'B'])

# Rocker angle delta
plt.figure()
plt.plot((q - q[0]) * 180. / np.pi, (mecA0.q4_list - mecA0.q4_list[0]) * 180. / np.pi)
plt.plot((q - q[0]) * 180. / np.pi, (mecB0.q4_list - mecB0.q4_list[0]) * 180. / np.pi)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$\Delta \\theta_4$ (deg)')
plt.grid()
plt.legend(['A', 'B'])

# Velocity coeffs
plt.figure()
plt.plot(q * 180. / np.pi, mecA0.k4_list)
plt.plot(q * 180. / np.pi, mecB0.k4_list)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$K_4 = 1 / VM$')
plt.legend(['A', 'B'])
plt.grid()

# VM
plt.figure()
plt.plot((q - q[0]) * 180. / np.pi, 1. / mecA0.k4_list)
plt.plot((q - q[0]) * 180. / np.pi, 1. / mecB0.k4_list)
plt.legend(['A', 'B'])
plt.xlabel('$\Delta \\theta_2$ (deg)')
plt.ylabel('$VM$')
plt.grid()

# Positions drawing
plt.figure()
for k in range(0, q.size, 10):
    mecA0.update_posi(q2A[k])
    mecA0.draw_posi()
    mecB0.update_posi(q2B[k])
    mecB0.draw_posi()
    
    R4A_c = r4c * vt.vunitary(mecA0.q4_list[k] + delta4A)
    R4B_c = r4c * vt.vunitary(mecB0.q4_list[k] + delta4B)
    
    # Cam follower
    plt.plot([mecA0.r1x, mecA0.r1x + R4A_c[0]],
             [mecA0.r1y, mecA0.r1y + R4A_c[1]], 'C2.--')
    plt.plot([mecB0.r1x, mecB0.r1x + R4B_c[0]],
             [mecB0.r1y, mecB0.r1y + R4B_c[1]], 'C2.--')

plt.text(0,0, '(0,0)')
plt.text(r1Ax, r1Ay, '(' + str(r1Ax) + ', ' + str(r1Ay) + ')')
plt.text(r1Bx, r1By, '(' + str(r1Bx) + ', ' + str(r1By) + ')')
plt.text(mecA0.R2[0] + 0.5 * mecA0.R3[0],
         mecA0.R2[1] + 0.5 * mecA0.R3[1], str(mecA0.r3))
plt.text(mecB0.R2[0] + 0.5 * mecB0.R3[0],
         mecB0.R2[1] + 0.5 * mecB0.R3[1], str(mecB0.r3))
plt.minorticks_on()
plt.grid(b = True, which = 'both', axis='both')
plt.axis('equal')
plt.xlabel('$x$ (mm) - hor')
plt.ylabel('$y$ (mm) - vert')