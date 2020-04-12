# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 12:14:05 2020

@author: UPB
"""
# =============================================================================
# Ventilator kinematic analusys (Module)
# =============================================================================
# This code defines a class for computing a kinematic analysis of a RRRR
# mechanism
# 
# Author: Juan A. Ramírez-Macías, PhD
# email: juan.ramirez@upb.edu.co
# Universidad Pontificia Bolivariana
# Medellín, COLOMBIA
# 2020
# =============================================================================


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
    return np.arctan2(R[1,0], R[0,0])

# Magnitude of a vector for a dim 2 array of two vectors
def vmag(R):
    return np.sqrt(R[1,0] ** 2 + R[0,0] ** 2)

def vmagangle(R):
    return vmag(R), vangle(R)

# Director cosines of a vector for a dim 2 array of two vectors
def vunitary(angle):
    return np.array([[np.cos(angle)], [np.sin(angle)]])

# RRR dyad solucion algorithm
def rrr(q2, r1x, r1y, r2, r3, r4, inv = 1.):
    # Input vector
    e2 = vunitary(q2)
    R2 = r2 * e2
    
    # Rocker position
    H = - np.array([[r1x], [r1y]]) + R2
    h, dh = vmagangle(H)    
    delta = tcos(r3, h, r4)    
    q4 = dh + inv * delta
    e4 = vunitary(q4)
    R4 = r4 * e4
    
    # Coupler position
    R3 = R4 - H
    q3 = vangle(R3)
    
    return q3, q4, R2, R3, R4

class RRR:
    def __init__(self, r2, r3, r4, r1x, r1y, inv = 1.):
        self.r1x = r1x
        self.r1y = r1y
        self.r2 = r2
        self.r3 = r3
        self.r4 = r4
        self.inv = inv # this is mechanism inversion
        
    def update_posi(self, q2):
        # Here the RRR dyad solution is repeated
        self.q2 = q2
        # Input vector
        self.R2 = self.r2 * vunitary(q2)
        
        # Rocker position
        H = - np.array([[self.r1x], [self.r1y]]) + self.R2
        
        h, dh = vmagangle(H)
        
        delta = tcos(self.r3, h, self.r4)    
        self.q4 = dh + self.inv * delta
        
        self.R4 = self.r4 * vunitary(self.q4)
        
        # Coupler position
        self.R3 = self.R4 - H
        self.q3 = vangle(self.R3)
        
    def posi_list(self, q2_list):
        N = q2_list.size
        self.R2_list = np.zeros([2, N])
        self.R3_list = np.zeros([2, N])
        self.R4_list = np.zeros([2, N])
        self.q2_list = q2_list
        self.q3_list = np.zeros_like(q2_list)
        self.q4_list = np.zeros_like(q2_list)
        
        for k in range(q2_list.size):
            self.update_posi(q2_list[k])
            self.R2_list[:, k] = self.R2[:,0]
            self.R3_list[:, k] = self.R3[:,0]
            self.R4_list[:, k] = self.R4[:,0]
            self.q3_list[k] = self.q3
            self.q4_list[k] = self.q4
        
    # =========================================================================
    # Velocity analysis (velocity coefficients)
    # =========================================================================
    def vel_coeffs(self):
        # Determinant of the Jacobian
        det = self.r3 * np.sin(self.q3) * self.r4 * np.cos(self.q4) \
            - self.r3 * np.cos(self.q3) * self.r4 * np.sin(self.q4)
        
        self.k3 = (-self.r2 * np.sin(self.q2) *  self.r4 * np.cos(self.q4)
                   + self.r2 * np.cos(self.q2) * self.r4 * np.sin(self.q4)) \
                    / det
        
        self.k4 = (self.r3 * np.sin(self.q3) * self.r2 * np.cos(self.q2)
                   - self.r3 * np.cos(self.q3) * self.r2 * np.sin(self.q2)) \
                    / det
                    
    def coeff_list(self, q2_list):
        N = q2_list.size
        self.k3_list = np.zeros(N)
        self.k4_list = np.zeros(N)
        
        for k in range(q2_list.size):
            self.update_posi(q2_list[k])
            self.vel_coeffs()
            self.k3_list[k] = self.k3
            self.k4_list[k] = self.k4
            
    def draw_posi(self):
        # 2
        self.h2 = plt.plot([0, self.R2[0]], [0, self.R2[1]], 'C0.-')
        
        # 3
        self.h3 = plt.plot([self.R2[0], self.R2[0] + self.R3[0]],
                           [self.R2[1], self.R2[1] + self.R3[1]], 'C1.-')
        
        # 4
        self.h4 = plt.plot([self.r1x, self.R2[0] + self.R3[0]],
                           [self.r1y, self.R2[1] + self.R3[1]], 'C2.-')    

# =============================================================================
# Main - test code
# =============================================================================
if __name__ == '__main__':
    # =============================================================================
    # Object creation
    # =============================================================================
    mecA0 = RRR(r2 = 25., r3 = 165., r4 = 40., r1x = 167.5, r1y = 19., inv = 1.)
    mecB0 = RRR(r2 = 25., r3 = 154., r4 = 40., r1x = 11., r1y = 147.4, inv = -1.)

    # Cam parameters
    r4c = 150. # cam length (drawing only)
    # angle between rocker and cam
    delta4A = 150. * np.pi / 180.
    delta4B = 130. * np.pi / 180.

    # =============================================================================
    # Position analysis - compute a set of angles
    # =============================================================================
    # Input angle config
    delta_ang = 1.
    q = np.arange(0., 45. + delta_ang, delta_ang) * np.pi / 180. 
        
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
    plt.plot(q * 180. / np.pi, (mecA0.q4_list - mecA0.q4_list[0]) * 180. / np.pi)
    plt.plot(q * 180. / np.pi, (mecB0.q4_list - mecB0.q4_list[0]) * 180. / np.pi)
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
    
    # VM
    plt.figure()
    plt.plot(q * 180. / np.pi, 1. / mecA0.k4_list, '.')
    plt.plot(q * 180. / np.pi, 1. / mecB0.k4_list, '+')
    plt.legend(['A', 'B'])
    plt.xlabel('$\\theta_2$ (deg)')
    plt.ylabel('$VM$')
    
    # Positions drawing
    plt.figure()
    for k in range(0, q.size, 10):
        mecA0.update_posi(q2A[k])
        mecA0.draw_posi()
        mecB0.update_posi(q2B[k])
        mecB0.draw_posi()
        
        R4A_c = r4c * vunitary(mecA0.q4_list[k] - delta4A)
        R4B_c = r4c * vunitary(mecB0.q4_list[k] + delta4B)
        
        # Cam follower
        plt.plot([mecA0.r1x, mecA0.r1x + R4A_c[0]],
                 [mecA0.r1y, mecA0.r1y + R4A_c[1]], 'C2.--')
        plt.plot([mecB0.r1x, mecB0.r1x + R4B_c[0]],
                 [mecB0.r1y, mecB0.r1y + R4B_c[1]], 'C2.--')
        
    plt.grid()
    plt.axis('equal')
    plt.xlabel('$x$ (mm) - hor')
    plt.ylabel('$y$ (mm) - vert')