# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 12:14:05 2020

@author: UPB
"""

# =============================================================================
# Modules
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Helper functions
# =============================================================================
# Triangle Cosine theorem
def tcos(a, b, c):
    return np.arccos((b ** 2 + c ** 2 - a ** 2) / 2. / b / c)

# Angle of a vector for a dim 2 array of two vectors
def vangle(R):
    return np.arctan2(R[1,:], R[0,:])

# Magnitude of a vector for a dim 2 array of two vectors
def vmag(R):
    return np.sqrt(R[1,:] ** 2 + R[0,:] ** 2)

def vmagangle(R):
    return vmag(R), vangle(R)

# Director cosines of a vector for a dim 2 array of two vectors
def vunitary(angle):
    return np.concatenate((np.cos(angle), np.sin(angle)), axis = 0)

def rrr(q2, r1x, r1y, r2, r3, r4, inv = 1.):
    # Input vector
    e2 = vunitary(q2)
    R2 = r2 * e2
    
    # Rocker position
    H = - np.array([[r1x], [r1y]]) + R2
    h, dh = vmagangle(H)    
    delta = tcos(r3, h, r4)    
    q4 = (dh + inv * delta).reshape(1, N)    
    e4 = vunitary(q4)    
    R4 = r4 * e4
    
    # Coupler position
    R3 = R4 - H
    q3 = vangle(R3)
    
    return q3, q4, R2, R3, R4

def rrr_vec(q2, r1x, r1y, r2, r3, r4):
    pass
# =============================================================================
# Parameters
# =============================================================================
# Input
r2 = 25.
# Coupler
r3A = 164.
r3B = 154.
# Rocker
r4 = 40.
# Cam
r4c = 150.
delta4A = 150. * np.pi / 180.
delta4B = 130. * np.pi / 180.
# Rocker axis position
# A
r1Ax = 167.5
r1Ay = 19.
# B
r1Bx = 11.
r1By = 147.4

# =============================================================================
# Position analysis
# =============================================================================
# Input angle config
N = 51
q = np.linspace(0., 45., N).reshape([1, N]) * np.pi / 180. 

# Input link position
q2A = q - 3. * np.pi / 4.
q2B = q + 3. * np.pi / 4.

q3A, q4A, R2A, R3A, R4A = rrr(q2A, r1Ax, r1Ay, r2, r3A, r4)
q3B, q4B, R2B, R3B, R4B = rrr(q2B, r1Bx, r1By, r2, r3B, r4, inv = -1.)

# Cam
R4A_c = r4c * vunitary(q4A - delta4A)
R4B_c = r4c * vunitary(q4B + delta4B)

print('delta q4A:', q4A[0,0] * 180. / np.pi, q4A[0,-1] * 180. / np.pi)
print('delta q4B:', q4B[0,0] * 180. / np.pi, q4B[0,-1] * 180. / np.pi)
print('delta q4A:', (q4A[0,-1] - q4A[0,0]) * 180. / np.pi)
print('delta q4B:',(q4B[0,-1] - q4B[0,0]) * 180. / np.pi)
print('q4A_m:', 0.5 * (q4A[0,-1] + q4A[0,0]) * 180. / np.pi)
print('q4B_m:', 0.5 * (q4B[0,-1] + q4B[0,0]) * 180. / np.pi)


# =============================================================================
# Velocity analysis
# =============================================================================
detA = r3A * np.sin(q3A) * r4 * np.cos(q4A) - r3A * np.cos(q3A) * r4 * np.sin(q4A)
detB = r3B * np.sin(q3B) * r4 * np.cos(q4B) - r3B * np.cos(q3B) * r4 * np.sin(q4B)

k3A = (-r2 * np.sin(q2A) *  r4 * np.cos(q4A) + r2 * np.cos(q2A) * r4 * np.sin(q4A)) / detA
k3B = (-r2 * np.sin(q2B) *  r4 * np.cos(q4B) + r2 * np.cos(q2B) * r4 * np.sin(q4B)) / detA

k4A = (r3A * np.sin(q3A) * r2 * np.cos(q2A) - r3A * np.cos(q3A) * r2 * np.sin(q2A)) / detA
k4B = (r3B * np.sin(q3B) * r2 * np.cos(q2B) - r3B * np.cos(q3B) * r2 * np.sin(q2B)) / detB

# =============================================================================
# Plot
# =============================================================================
plt.close('all')

plt.figure()
plt.plot(q[0,:] * 180. / np.pi, q4A[0,:] * 180. / np.pi)
plt.plot(q[0,:] * 180. / np.pi, q4B[0,:] * 180. / np.pi)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$\\theta_4$ (deg)')
plt.grid()
plt.legend(['A', 'B'])

plt.figure()
plt.plot(q[0,:] * 180. / np.pi, (q4A[0,:] - q4A[0,0]) * 180. / np.pi)
plt.plot(q[0,:] * 180. / np.pi, (q4B[0,:] - q4B[0,0]) * 180. / np.pi)
plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$\Delta \\theta_4$ (deg)')
plt.grid()
plt.legend(['A', 'B'])

# plt.figure()
# plt.plot(q[0,:] * 180. / np.pi, q3A[0,:] * 180. / np.pi)
# plt.plot(q[0,:] * 180. / np.pi, q3B[0,:] * 180. / np.pi - 90)

plt.figure()
plt.plot(q[0,:] * 180. / np.pi, k4A[0,:])
plt.plot(q[0,:] * 180. / np.pi, k4B[0,:])

plt.plot(q[0,:] * 180. / np.pi, np.gradient(q4A[0,:], q2A[0,:]), ':')
plt.plot(q[0,:] * 180. / np.pi, np.gradient(q4B[0,:], q2B[0,:]), ':')

plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$K_4 = 1 / VM$')
plt.legend(['A', 'B'])

plt.figure()
plt.plot(q[0,:] * 180. / np.pi, 1 / k4A[0,:], '.')
plt.plot(q[0,:] * 180. / np.pi, 1 / k4B[0,:], '+')
plt.legend(['A', 'B'])

# plt.plot(q[0,:] * 180. / np.pi, 1 / np.gradient(q4A[0,:], q2A[0,:]), ':')
# plt.plot(q[0,:] * 180. / np.pi, 1 / np.gradient(q4B[0,:], q2B[0,:]), ':')

plt.xlabel('$\\theta_2$ (deg)')
plt.ylabel('$VM$')
# plt.ylim([1., 3.])


plt.figure()
for k in range(0, N, 10):
    # 2
    plt.plot([0, R2A[0, k]], [0, R2A[1, k]], 'C0.--')
    plt.plot([0, R2B[0, k]], [0, R2B[1, k]], 'C0.-')
    
    # 3
    plt.plot([R2A[0, k], R2A[0, k] + R3A[0, k]],
             [R2A[1, k], R2A[1, k] + R3A[1, k]], 'C1.--')
    plt.plot([R2B[0, k], R2B[0, k] + R3B[0, k]],
             [R2B[1, k], R2B[1, k] + R3B[1, k]], 'C1.-')
    
    # 4
    plt.plot([r1Ax, R2A[0, k] + R3A[0, k]],
             [r1Ay, R2A[1, k] + R3A[1, k]], 'C2.--')
    plt.plot([r1Bx, R2B[0, k] + R3B[0, k]],
             [r1By, R2B[1, k] + R3B[1, k]], 'C2.-')
    
    # Cam follower
    plt.plot([r1Ax, r1Ax + R4A_c[0, k]],
             [r1Ay, r1Ay + R4A_c[1, k]], 'C2.--')
    plt.plot([r1Bx, r1Bx + R4B_c[0, k]],
             [r1By, r1By + R4B_c[1, k]], 'C2.--')
    
plt.grid()
plt.axis('equal')
plt.xlabel('$x$ (mm) - hor')
plt.ylabel('$y$ (mm) - vert')