import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from CSubd2D import *
from BezierAvg import *
from BSplineAvg import *

#-----------------------------------------------------------------------------
def bopt_average_export(t0, t1, b_open, p0, p1, n0, n1):
    p, n = bopt_average_v1(t0, p0, p1, n0, n1)
    return p, n, p, 0,0,0

#-----------------------------------------------------------------------------
def bopt_average_v1(t0, p0, p1, n0, n1):
    ca_pt, ca_norm, ca_center, ca_radius, ca_beta0, ca_beta1 = \
      circle_avg(0.5, 0.5, True, p0, p1, n0, n1)

    theta = get_angle_between(n0, n1)
    theta_deg = theta * 180.0 / np.pi
    p0_p1_dist = get_dist( p0, p1 )
    if theta_deg < 5. or p0_p1_dist < 0.001 :
        rn = p1-p0
        res_norm = np.array([rn[1], -rn[0]])
        res_norm /= np.linalg.norm(res_norm)
        return ca_pt, res_norm

    der0_dir = np.array([n0[1], -n0[0]])
    der1_dir = np.array([n1[1], -n1[0]])

    initial_len = 4./3. * np.tan(theta/4.) * ca_radius
    
    def err_func(dir_length):
        a = p0
        b = p0 - der0_dir * dir_length[0]
        c = p1 + der1_dir * dir_length[0]#1]
        d = p1
        
        t_tange = np.linspace(0.0, 1.0, 10)
        err = 0.
        for curr_t in t_tange:
            curr_pt, _ = eval_cubic_bezier( a, b, c, d, curr_t )
            dist_to_ca_center = get_dist(curr_pt, ca_center)
            err += np.abs(dist_to_ca_center - ca_radius)

        return err

    x0 = np.array([initial_len])#, initial_len])
    res = minimize(err_func, x0, method='nelder-mead',
                   options={'disp': False})

    obtained_len0 = np.abs( res.x[0])
    obtained_len1 = np.abs( res.x[0])#1])
    #print 'Init length =', initial_len, 'obtained length =', obtained_len
    a = p0
    b = p0 - der0_dir * obtained_len0
    c = p1 + der1_dir * obtained_len1
    d = p1



    res_t = t0
    #if t0 < 0.:
    #    d, c, b, a = extrapolate_bezier(d, c, b, a, np.abs(t0))
    #    #a, b, c, d = extrapolate_bezier(d, c, b, a, np.abs(t0))
    #    res_t = 0.0
    #elif t0 > 1.:
    #    a, b, c, d = extrapolate_bezier(a, b, c, d, t0-1.)
    #    res_t = 1.

    res_pt, res_der = eval_cubic_bezier( a, b, c, d, 1 - res_t )
    res_norm = np.array([res_der[1], -res_der[0]])
    res_norm /= np.linalg.norm(res_norm)
    #if np.linalg.norm(res_norm - ca_norm) > \
    #   np.linalg.norm(res_norm + ca_norm):
    #    res_norm = -res_norm    

    DEBUG = 'IN_DEBUG' in globals()
    if DEBUG and IN_DEBUG:
        #appr_pt, appr_der = eval_cubic_bezier( a, b_orig, c_orig, d, t0 )
        #print 'Ver 2 Approxim   = (', appr_pt[0], ', ', appr_pt[1], ')'
        #print 'Ver 2 Point      = (', res_pt[0],', ',res_pt[1], ')'
        #print 'Ver 2 Normal     = (', res_norm[0],', ',res_norm[1], ')'
        #print 'Ver 2 Der Length = ', der_length
        verts = [
            (a[0], a[1]), # P0
            (b[0], b[1]), # P1
            (c[0], c[1]), # P2
            (d[0], d[1]), # P3
            ]
        
        codes = [Path.Path.MOVETO,
                 Path.Path.CURVE4,
                 Path.Path.CURVE4,
                 Path.Path.CURVE4,
                 ]
        
        crv = Path.Path(verts, codes)
        
        fig = plt.figure() 
        plt.xlim([(ca_center[0] - ca_radius)*1.2, (ca_center[0] + ca_radius)*1.2])
        plt.ylim([(ca_center[1] - ca_radius)*1.2, (ca_center[1] + ca_radius)*1.2])
        cr1 = plt.Circle( (ca_center[0], ca_center[1]), radius=ca_radius, fc='y', ec='none')
        nr000 = plt.Arrow(p0[0], p0[1], 
                          n0[0]*.1, n0[1]*.1, width=0.03, fc='b', ec='none' )
        nr100 = plt.Arrow(p1[0], p1[1], 
                          n1[0]*.1, n1[1]*.1, width=0.03, fc='g', ec='none' )
        nr_res = plt.Arrow(res_pt[0], res_pt[1], 
                          res_norm[0]*.1, res_norm[1]*.1, width=0.03, fc='r', ec='none' )
        plt.plot([a[0], b[0], c[0], d[0]], 
                 [a[1], b[1], c[1], d[1]], '-c' )
        plt.gca().add_patch(cr1)
        plt.gca().add_patch(nr000)
        plt.gca().add_patch(nr_res)
        plt.gca().add_patch(nr100)
        patch = patches.PathPatch(crv, facecolor='none', lw=1)
        plt.gca().add_patch(patch)
        plt.axis('equal')
        plt.show()
    return res_pt, res_norm

#-----------------------------------------------------------------------------
def one_pair_test():
    points  = [(  0.0, 0.0),
               (  0.0, 1.0)]
    normals = [(  1.0, 0.0),
               (  1.0, -1.0)]
    normals[0] /= np.linalg.norm(normals[0])
    normals[1] /= np.linalg.norm(normals[1])
    p0 = np.array([points[0][0], points[0][1]])
    p1 = np.array([points[1][0], points[1][1]])
    n0 = np.array([normals[0][0], normals[0][1]])
    n1 = np.array([normals[1][0], normals[1][1]])
    #bopt_average_v1(0.5, p0, p1, n0, n1)
    #bopt_average_v1(-0.125, p0, p1, n0, n1)
    bopt_average_v1(0.5, p0, p1, n0, n1)
    #bopt_average_v1(1.125, p0, p1, n0, n1)
#-----------------------------------------------------------------------------
if __name__ == "__main__":

    global IN_DEBUG

    IN_DEBUG = True
    IN_DEBUG = False
    #one_pair_test()
    build_curves()

#============================ END OF FILE ====================================
