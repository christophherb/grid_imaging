import numpy as np
import matplotlib.pyplot as plt

def calc_ks(pz, pr2, LStart, LEnd):
    """return k1 k2 and k3 such that r = sqrt(k1+z*k2+z*z*k3) from a point pz, pr**2 on the mirror surface
    and the position of both focal points

    Args:
        pz (float): z-coordinate of the point on the mirror
        pr2 (float): x**2+y**2, square of the distance of the point on the mirror to the optical axis
        LStart (float): first focal point
        LEnd (float): second focal point

    Returns:
        tuple: (k1, k2, k3) tuple such that r = sqrt(k1+z*k2+z*z*k3)
"""
    c = (LEnd-LStart)/2
    u = (pz+c-LEnd)
    a = (u*u+c*c+pr2+((u*u+c*c+pr2)**2-4*c*c*u*u)**0.5)**0.5/2**0.5
    k3 = c*c/(a*a)-1
    k2 = 2*k3*(c-LEnd)
    k1 = k3*(c-LEnd)*(c-LEnd)-c*c+a*a
    return k1, k2, k3

def return_r(z, k1, k2, k3):
    return (k1+ k2*z+k3*z**2)**0.5



#------------------------ Testzone -----------#
def return_all_b0s(nummirrors, z0, r0, LStart, LEnd, lStart, lEnd):
    r_at_z0 = []
    r_at_zend = []
    pz0 = z0
    pr2 = r0**2
    for mirror_ind in range(nummirrors):
        k1, k2, k3 = calc_ks(pz0, pr2, LStart=LStart, LEnd=LEnd)
        r_at_z0.append((k1+z0*k2+z0**2*k3)**0.5)
        r2lEnd = (k1 + lEnd*k2 + lEnd*lEnd*k3)
        rlEnd = r2lEnd**0.5
        r_at_zend.append(rlEnd)
        rlStart = rlEnd*(lStart-LStart)/(lEnd-LStart)
        pr2 = rlStart**2
        pz0 = lStart
    return r_at_z0, r_at_zend

def return_maxN(min_distance, z0, r0, LStart, LEnd, lStart, lEnd):
    """returns the maximum number of mirrors for a given geometry such that the distance between two mirrors
    is no lower than min_distance

    Args:
        min_distance (float): minimum distance between mirrors (m)
        z0 (float): z0
        r0 (float): distance of the outermost mirror
        LStart (float): distance of the first focal point
        LEnd (float): distance of the second focal point
        lStart (float): lStart
        lEnd (float): lStart

    Returns:
        int: number of the mirrors 
    """
    r_at_z0 = []
    r_at_zend = []
    pz0 = z0
    pr2 = r0**2
    while True:
        k1, k2, k3 = calc_ks(pz0, pr2, LStart=LStart, LEnd=LEnd)
        r_at_z0.append((k1+z0*k2+z0**2*k3)**0.5)
        r2lEnd = (k1 + lEnd*k2 + lEnd*lEnd*k3)
        rlEnd = r2lEnd**0.5
        try:
            if abs(rlEnd-r_at_zend[-1]) < min_distance:
                break
            else: pass
        except IndexError:
            pass
        r_at_zend.append(rlEnd)
        rlStart = rlEnd*(lStart-LStart)/(lEnd-LStart)
        pr2 = rlStart**2
        #print('no reachy')
        pz0 = lStart
    return len(r_at_zend)

def return_max_r0_mindistance(min_distance, max_r0, z0, r0, LStart, LEnd, lStart, lEnd):
    """returns the maximum number of mirrors for a given geometry such that the distance between two mirrors
    is no lower than min_distance and the outermost mirror sits no further than max_distance

    Args:
        min_distance (float): minimum distance between mirrors (m)
        max_r0 (float): maximum outer distance mirrors (m)
        z0 (float): z0
        r0 (float): distance of the outermost mirror
        LStart (float): distance of the first focal point
        LEnd (float): distance of the second focal point
        lStart (float): lStart
        lEnd (float): lStart

    Returns:
        int: number of the mirrors 
    """
    r_at_z0 = []
    r_at_zend = []
    pz0 = z0
    pr2 = r0**2
    while True:
        k1, k2, k3 = calc_ks(pz0, pr2, LStart=LStart, LEnd=LEnd)
        r_at_z0.append((k1+z0*k2+z0**2*k3)**0.5)
        r2lEnd = (k1 + lEnd*k2 + lEnd*lEnd*k3)
        rlEnd = r2lEnd**0.5
        try:
            if abs(rlEnd-r_at_zend[-1]) < min_distance:
                break
            else: pass
        except IndexError:
            pass
        r_at_zend.append(rlEnd)
        rlStart = rlEnd*(lStart-LStart)/(lEnd-LStart)
        pr2 = rlStart**2
        #print('no reachy')
        pz0 = lStart
    r_at_z0 = np.array(r_at_z0)[np.array(r_at_z0)<max_r0]
    
    return len(r_at_z0), r_at_z0[0]

# calculation of all the employed mirrors in the ess extraction mirror code
#test of the scaling factor
#print(return_maxN(0.0005, 0, 0.02076, -0.6, 0.6, -0.06, 0.06))
#mirror 1
#mirror_ind = 1
#print('this is mirror {}'.format(mirror_ind))
#mirror_ind+=1
#
#scaling_factor = 2.0
#focal_length = 0.6
#z0 = 0
#r0 = 0.02076*scaling_factor
#LStart = -focal_length*scaling_factor
#LEnd = focal_length*scaling_factor
#lStart = -0.06*scaling_factor
#lEnd = 0.06*scaling_factor
#
#print([(k, k/scaling_factor) for k in return_all_b0s(10, z0, r0, LStart, LEnd, lStart, lEnd)[0]])
#print([(k, k/scaling_factor) for k in return_all_b0s(10, z0, r0, LStart, LEnd, lStart, lEnd)[1]])
##mirror 2
#scaling_factor = 1
#focal_length = 0.6*scaling_factor
#incoming_length = 0.6*scaling_factor
#z0 = 0
#lStart = -0.06*scaling_factor
#lEnd = 0.06*scaling_factor
#print('this is mirror {}'.format(mirror_ind))
#mirror_ind+=1
#r0 = 0.02076*scaling_factor
#LStart = -focal_length*scaling_factor
#LEnd = focal_length*scaling_factor
#print([k/scaling_factor for k in return_all_b0s(10, z0, r0, LStart, LEnd, lStart, lEnd)[0]])
#print([k/scaling_factor for k in return_all_b0s(10, z0, r0, LStart, LEnd, lStart, lEnd)[1]])
#
