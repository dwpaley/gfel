# -*- coding: utf-8 -*-
'''
*GSASIIlattice: Unit cells*
---------------------------

Perform lattice-related computations

Note that *G* is the reciprocal lattice tensor, and *g* is its inverse,
:math:`G = g^{-1}`, where 

  .. math::

   g = \\left( \\begin{matrix}
   a^2 & a b\\cos\gamma & a c\\cos\\beta \\\\
   a b\\cos\\gamma & b^2 & b c \cos\\alpha \\\\
   a c\\cos\\beta &  b c \\cos\\alpha & c^2
   \\end{matrix}\\right)

The "*A* tensor" terms are defined as
:math:`A = (\\begin{matrix} G_{11} & G_{22} & G_{33} & 2G_{12} & 2G_{13} & 2G_{23}\\end{matrix})` and *A* can be used in this fashion:
:math:`d^* = \sqrt {A_0 h^2 + A_1 k^2 + A_2 l^2 + A_3 hk + A_4 hl + A_5 kl}`, where
*d* is the d-spacing, and :math:`d^*` is the reciprocal lattice spacing, 
:math:`Q = 2 \\pi d^* = 2 \\pi / d`. 
Note that GSAS-II variables ``p::Ai`` (``i``=0,1,...5) and ``p`` is a phase number are 
used for the *Ai* values. See :func:`A2cell`, :func:`cell2A` for interconversion between A and 
unit cell parameters; :func:`cell2Gmat` :func:`Gmat2cell` for G and cell parameters. 

'''
from __future__ import division, print_function
import math
import sys
import random as ran
import numpy as np
import numpy.linalg as nl

from cctbx import sgtbx, uctbx, miller

# trig functions in degrees
cosd = lambda x: np.cos(x*np.pi/180.)
acosd = lambda x: 180.*np.arccos(x)/np.pi

def sec2HMS(sec):
    """Convert time in sec to H:M:S string
    
    :param sec: time in seconds
    :return: H:M:S string (to nearest 100th second)
    
    """
    H = int(sec//3600)
    M = int(sec//60-H*60)
    S = sec-3600*H-60*M
    return '%d:%2d:%.2f'%(H,M,S)
    
        
    
def fillgmat(cell):
    """Compute lattice metric tensor from unit cell constants

    :param cell: tuple with a,b,c,alpha, beta, gamma (degrees)
    :return: 3x3 numpy array

    """
    a,b,c,alp,bet,gam = cell
    g = np.array([
        [a*a,  a*b*cosd(gam),  a*c*cosd(bet)],
        [a*b*cosd(gam),  b*b,  b*c*cosd(alp)],
        [a*c*cosd(bet) ,b*c*cosd(alp),   c*c]])
    return g
           
def cell2Gmat(cell):
    """Compute real and reciprocal lattice metric tensor from unit cell constants

    :param cell: tuple with a,b,c,alpha, beta, gamma (degrees)
    :return: reciprocal (G) & real (g) metric tensors (list of two numpy 3x3 arrays)

    """
    g = fillgmat(cell)
    G = nl.inv(g)        
    return G,g

def A2Gmat(A,inverse=True):
    """Fill real & reciprocal metric tensor (G) from A.

    :param A: reciprocal metric tensor elements as [G11,G22,G33,2*G12,2*G13,2*G23]
    :param bool inverse: if True return both G and g; else just G
    :return: reciprocal (G) & real (g) metric tensors (list of two numpy 3x3 arrays)

    """
    G = np.array([
        [A[0],  A[3]/2.,  A[4]/2.], 
        [A[3]/2.,A[1],    A[5]/2.], 
        [A[4]/2.,A[5]/2.,    A[2]]])
    if inverse:
        g = nl.inv(G)
        return G,g
    else:
        return G

def Gmat2A(G):
    """Extract A from reciprocal metric tensor (G)

    :param G: reciprocal maetric tensor (3x3 numpy array
    :return: A = [G11,G22,G33,2*G12,2*G13,2*G23]

    """
    return [G[0][0],G[1][1],G[2][2],2.*G[0][1],2.*G[0][2],2.*G[1][2]]
    
def cell2A(cell):
    """Obtain A = [G11,G22,G33,2*G12,2*G13,2*G23] from lattice parameters

    :param cell: [a,b,c,alpha,beta,gamma] (degrees)
    :return: G reciprocal metric tensor as 3x3 numpy array

    """
    G,g = cell2Gmat(cell)
    return Gmat2A(G)

def A2cell(A):
    """Compute unit cell constants from A

    :param A: [G11,G22,G33,2*G12,2*G13,2*G23] G - reciprocal metric tensor
    :return: a,b,c,alpha, beta, gamma (degrees) - lattice parameters

    """
    G,g = A2Gmat(A)
    return Gmat2cell(g)

def Gmat2cell(g):
    """Compute real/reciprocal lattice parameters from real/reciprocal metric tensor (g/G)
    The math works the same either way.

    :param g (or G): real (or reciprocal) metric tensor 3x3 array
    :return: a,b,c,alpha, beta, gamma (degrees) (or a*,b*,c*,alpha*,beta*,gamma* degrees)

    """
    oldset = np.seterr('raise')
    a = np.sqrt(max(0,g[0][0]))
    b = np.sqrt(max(0,g[1][1]))
    c = np.sqrt(max(0,g[2][2]))
    alp = acosd(g[2][1]/(b*c))
    bet = acosd(g[2][0]/(a*c))
    gam = acosd(g[0][1]/(a*b))
    np.seterr(**oldset)
    return a,b,c,alp,bet,gam

            
def calc_rVsq(A):
    """Compute the square of the reciprocal lattice volume (1/V**2) from A'

    """
    G,g = A2Gmat(A)
    rVsq = nl.det(G)
    if rVsq < 0:
        return 1
    return rVsq
    
def calc_rV(A):
    """Compute the reciprocal lattice volume (V*) from A
    """
    return np.sqrt(calc_rVsq(A))
    
def calc_V(A):
    """Compute the real lattice volume (V) from A
    """
    return 1./calc_rV(A)
    

#reflection generation routines
#for these: H = [h,k,l]; A is as used in calc_rDsq; G - inv metric tensor, g - metric tensor; 
#           cell - a,b,c,alp,bet,gam in A & deg
                   
def calc_rDsq(H,A):
    'needs doc string'
    rdsq = H[0]*H[0]*A[0]+H[1]*H[1]*A[1]+H[2]*H[2]*A[2]+H[0]*H[1]*A[3]+H[0]*H[2]*A[4]+H[1]*H[2]*A[5]
    return rdsq

def make_sgtype(ibrav):
    symmorphic_sgs = ['F23', 'I23', 'P23', 'R3', 'P3', 'I4', 'P4', 'F222', 
        'I222', 'A222', 'B222', 'C222', 'P222', 'I2', 'C2', 'P2', 'P1']
    return sgtbx.space_group_type(symmorphic_sgs[ibrav])

def GenHBravais(dmin,Bravais,A, sg_type=None):
    """Generate the positionally unique powder diffraction reflections
     
    :param dmin: minimum d-spacing in A
    :param Bravais: lattice type (see GetBraviasNum). Bravais is one of:
    
            * 0 F cubic
            * 1 I cubic
            * 2 P cubic
            * 3 R hexagonal (trigonal not rhombohedral)
            * 4 P hexagonal
            * 5 I tetragonal
            * 6 P tetragonal
            * 7 F orthorhombic
            * 8 I orthorhombic
            * 9 A orthorhombic
            * 10 B orthorhombic
            * 11 C orthorhombic
            * 12 P orthorhombic
            * 13 I monoclinic
            * 14 C monoclinic
            * 15 P monoclinic
            * 16 P triclinic
            
    :param A: reciprocal metric tensor elements as [G11,G22,G33,2*G12,2*G13,2*G23]
    :param st_type: an sgtbx.space_group_type object. Constructing these is slow
      so it's good to precalculate if possible.
    :return: HKL unique d list of [h,k,l,d,-1] sorted with largest d first
            
    """
    g_inv = np.array([[A[0],   A[3]/2, A[4]/2],
                      [A[3]/2, A[1],   A[5]/2], 
                      [A[4]/2, A[5]/2, A[2]]])
    g = np.linalg.inv(g_inv)
    g_elems = (g[0][0], g[1][1], g[2][2], g[0][1], g[0][2], g[1][2])
    try:
        uc = uctbx.unit_cell(metrical_matrix=g_elems)
    except ValueError: # this function sometimes receives an A matrix that gives
                       # numbers <0 in the diagonal elems of g. Not sure why.
        return []
    if sg_type is None:
        sg_type = make_sgtype(Bravais)
    mig = miller.index_generator(uc, sg_type, 0, dmin)
    result = []
    for h,k,l in mig: 
      d = uc.d((h,k,l))
      result.append([h, k, l, d, -1])
    result.sort(key=lambda l: l[3], reverse=True)
    return result
    
    
    

    
    
        

        
    
    



    
   
    

    
    
    
    
    
    
        
        
    

if __name__ == '__main__':
    # run self-tests
    selftestquiet = False
    for test in selftestlist:
        test()
    print ("OK")
