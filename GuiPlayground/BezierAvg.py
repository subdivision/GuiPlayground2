import numpy as np
from scipy.misc import comb
from scipy import optimize
from svg.path import CubicBezier, QuadraticBezier

from matplotlib import pyplot as plt
import matplotlib.path as Path
import matplotlib.patches as patches
from CSubd2D import *

#//---------------------------------------------------------------------------
#// NOTES:       TOLERANCE is a maximum error ratio
#//                      if n_limit isn't a power of 2 it will be act like the next higher
#//                      power of two.
def Simpson (f, a, b, n_limit, TOLERANCE):
    n = 1;
    multiplier = (b - a)/6.0
    endsum = f(a) + f(b)
    interval = (b - a)/2.0
    asum = 0
    bsum = f(a + interval)
    est1 = multiplier * (endsum + 2 * asum + 4 * bsum)
    est0 = 2 * est1

    while (n < n_limit and (abs(est1) > 0 and abs((est1 - est0) / est1) > TOLERANCE)):
        n *= 2
        multiplier /= 2;
        interval /= 2;
        asum += bsum;
        bsum = 0;
        est0 = est1;
        interval_div_2n = interval / (2.0 * n);

        for i in range(1, 2 * n, 2):
            t = a + i * interval_div_2n
            bsum += f(t)
        
        est1 = multiplier*(endsum + 2*asum + 4*bsum);

    return est1;

##---------------------------------------------------------------------------
def BezierArcLength(p1, p2, p3, p4):
    k1 = -p1 + 3*(p2 - p3) + p4;
    k2 = 3*(p1 + p3) - 6*p2;
    k3 = 3*(p2 - p1);
    k4 = p1;

    q1 = 9.0  * ((k1[0])**2 + (k1[1])**2)
    q2 = 12.0 * (k1[0]*k2[0] + k1[1]*k2[1]);
    q3 = 3.0*(k1[0]*k3[0] + k1[1]*k3[1]) + 4.0*(k2[0]**2 + k2[1]**2)
    q4 = 4.0*(k2[0]*k3[0] + k2[1]*k3[1])
    q5 = k3[0]**2 + k3[1]**2
    try:
        balf = lambda t: (q5 + t*(q4 + t*(q3 + t*(q2 + t*q1))))**0.5
    except RuntimeError as e:
        print e

    result = Simpson(balf, 0., 1., 1024, 0.001)
    return result;

#------------------------------------------------------------------------------
def eval_quadr_bezier(a, b, c, t):
  t2 = t * t
  mt = 1 - t
  mt2 = mt * mt
  return a*mt2 + b*2*mt*t + c*t2, 2*(a-b)*mt + 2*(b-c)*t

#------------------------------------------------------------------------------
def split_cubic_bezier_get_left(a, p, q, g, t):
    ''' 
        http://math.stackexchange.com/questions/877725/\
        retrieve-the-initial-cubic-b%C3%A9zier-curve-\
        subdivided-in-two-b%C3%A9zier-curves
    '''
    u = t
    v = 1. - t
    b = u*p + v*a
    f = u*g + v*q
    r = u*q + v*p
    c = u*r + v*b
    e = u*f + v*r
    d = u*e + v*c
    new_a, new_b, new_c, new_d = np.copy(a), b, c, d
    return new_a, new_b, new_c, new_d

#------------------------------------------------------------------------------
def split_cubic_bezier_get_right(a, p, q, g, t):
    ''' 
        http://math.stackexchange.com/questions/877725/\
        retrieve-the-initial-cubic-b%C3%A9zier-curve-\
        subdivided-in-two-b%C3%A9zier-curves
    '''
    u = t
    v = 1. - t
    b = u*p + v*a
    f = u*g + v*q
    r = u*q + v*p
    c = u*r + v*b
    e = u*f + v*r
    d = u*e + v*c
    new_a, new_b, new_c, new_d = d, e, f, np.copy(g)
    return new_a, new_b, new_c, new_d

#------------------------------------------------------------------------------
def eval_cubic_bezier(a, b, c, d, t):
    t2 = t * t
    t3 = t2 * t
    mt = 1-t
    mt2 = mt * mt
    mt3 = mt2 * mt
    der1, _ = eval_quadr_bezier(3.*(b-a), 3.*(c-b), 3.*(d-c), t)
    return a*mt3 + b*3*mt2*t + c*3*mt*t2 + d*t3, der1

#----------------------------------------------------------------------------
def estimate_bezier_length(a, b, c, d):
    ad = np.linalg.norm(a-d)
    abcd = np.linalg.norm(a-b) + np.linalg.norm(b-c) + np.linalg.norm(c-d)
    return (ad + abcd) / 2.

#----------------------------------------------------------------------------
def estimate_bezier_length_v2(a, b, c, d):
    prev_pt = np.array([np._NoValue, np._NoValue])
    res = 0.0

    xs = []
    ys = []
    for t in np.linspace(0.0, 1.0, 100):
        curr_pt, _ = eval_cubic_bezier(a,b,c,d,t)
        if np._NoValue != prev_pt[0]:
            res += np.linalg.norm(prev_pt - curr_pt)
            xs.append(curr_pt[0])
            ys.append(curr_pt[1])
        prev_pt = curr_pt

    #plt.plot(xs, ys)
    #plt.show()

    return res

#----------------------------------------------------------------------------
def get_t_for_length(a, b, c, d, lng):
    t_min, t_max = 0., 1.
    res_t = -1.
    while t_min < t_max:
        curr_t = (t_max + t_min)/2.
        curr_a, curr_b, curr_c, curr_d = split_cubic_bezier_get_left(a,b,c,d, curr_t)
        curr_lng = estimate_bezier_length_v2(curr_a, curr_b, curr_c, curr_d)
        if np.abs(curr_lng-lng) < 0.001:
            res_t = curr_t 
            break
        elif curr_lng > lng:
            t_max = curr_t 
        else:
            t_min = curr_t
    return res_t

#----------------------------------------------------------------------------
def check_derivative_length_factor(a, b, c, d, der_factor, t):
    crv0_len = estimate_bezier_length_v2(a,b,c,d)
    a1,b1,c1,d1 = split_cubic_bezier_get_left(a, b, c, d, t)
    crv1_len = estimate_bezier_length_v2(a1,b1,c1,d1)
    res0_pt, res0_norm = eval_cubic_bezier(a, b, c, d, t/2.)
    res0_norm /= np.linalg.norm(res0_norm)
    res1_pt, res1_norm = eval_cubic_bezier(a1, b1, c1, d1, t)
    res1_norm /= np.linalg.norm(res1_norm)
    prim_vs_sec_pt_dist = np.linalg.norm(res0_pt - res1_pt)
    prim_vs_sec_norm_dist = np.linalg.norm(res0_norm - res1_norm)
    print 'Factor ', der_factor, ' pt-to-pt ', prim_vs_sec_pt_dist, \
          ' norm-to-norm ', prim_vs_sec_norm_dist,
    print ' crv0_len ', crv0_len, ' crv1_len ', crv1_len


#----------------------------------------------------------------------------
def bezier_average_export( t, p0, p1, n0, n1 ):
    a = p0
    g = p1
    d0 = np.array([n0[1], -n0[0]])    
    d1 = np.array([n1[1], -n1[0]])    
    p = a + d0
    q = g - d1

    new_a, new_b, new_c, new_d = split_cubic_bezier_get_left(a, p, q, g, t)
    #new_d, new_e, new_f, new_g = split_cubic_bezier_get_right(a, p, q, g, t)
    res_pt = new_d
    res_norm = new_c - new_d
    res_norm = np.array([res_norm[1], -res_norm[0]])
    res_norm /= np.linalg.norm(res_norm)
    print '---------------'
    print 't = ', t
    print 'n0 = ', n0
    print 'n1 = ', n1
    print 'nR = ', res_norm
    print 'n0 <> n1', get_angle_between(n0, n1)
    print 'n0 <> res_norm ', get_angle_between(n0, res_norm)
    print 'n1 <> res_norm ', get_angle_between(n1, res_norm)
    print 'der0 length ', np.linalg.norm(d0)
    print 'der1 length ', np.linalg.norm(d1)
    print 'res der0 length ', np.linalg.norm(new_b - new_a)
    print 'res der1 length ', np.linalg.norm(new_c - new_d)

    DEBUG = True
    #DEBUG = False
    if DEBUG:
        verts = [
            (a[0], a[1]), # P0
            (p[0], p[1]), # P1
            (q[0], q[1]), # P2
            (g[0], g[1]), # P3
            ]
        
        codes = [Path.Path.MOVETO,
                 Path.Path.CURVE4,
                 Path.Path.CURVE4,
                 Path.Path.CURVE4,
                 ]
        
        crv = Path.Path(verts, codes)
        
        fig = plt.figure()
        dr0 = plt.Arrow( p0[0], p0[1], d0[0]*.1, d0[1]*.1, width=0.03, fc='g', ec='none' )
        dr1 = plt.Arrow( p1[0], p1[1], d1[0]*.1, d1[1]*.1, width=0.03, fc='g', ec='none' )
        dn0 = plt.Arrow( p0[0], p0[1], n0[0]*.1, n0[1]*.1, width=0.03, fc='b', ec='none' )
        dn1 = plt.Arrow( p1[0], p1[1], n1[0]*.1, n1[1]*.1, width=0.03, fc='b', ec='none' )
        dn2 = plt.Arrow( res_pt[0], res_pt[1], res_norm[0]*.1, res_norm[1]*.1, 
                        width=0.05, fc='r', ec='none' )
        plt.gca().add_patch(dr0)
        plt.gca().add_patch(dr1)
        plt.gca().add_patch(dn0)
        plt.gca().add_patch(dn1)
        plt.gca().add_patch(dn2)
        patch = patches.PathPatch(crv, facecolor='none', lw=2)
        plt.gca().add_patch(patch)
        plt.axis('equal')
        plt.axis([-1.0, 2., -1.0, 2])
        plt.show()


    return res_pt, res_norm

#----------------------------------------------------------------------------
if __name__ == "__main__":
    #np.seterr(all='raise')
    p0 = np.array([ 0.0, 0.0])
    p1 = np.array([ 1.0, 0.0])

    n0 = np.array([-1.0, 1.0])
    n1 = np.array([ 1.0, 1.0])

    #n0 = np.array([ 1.0,-1.0])
    #n1 = np.array([ 1.0, 1.0])

    #n0 = np.array([-1.0, 1.0])
    #n1 = np.array([ 1.0, 1.0])

    n0 /= np.linalg.norm(n0)
    n1 /= np.linalg.norm(n1)
    p_res05,   n_res05   = bezier_average_export(0.5,  p0, p1,      n0, n1)
    p_res0505, n_res0505 = bezier_average_export(0.5,  p0, p_res05, n0, n_res05)
    p_res025,  n_res025  = bezier_average_export(0.25, p0, p1,      n0, n1)
    print 'P_05_05 = ', p_res0505, 'N_05_05 = ', n_res0505
    print 'P_025   = ', p_res025,  'N_025   = ', n_res025

#----------------------------------------------------------------------------
#def bezier_average_export_otimizer_version( t0, t1, b_open, p0, p1, n0, n1 ):
#    middle_pt = ( p1 + p0 )/2.0
#    p0_to_p1_dist = np.linalg.norm(p0-p1)        
#    a = p0
#    d = p1
#    der_0 = np.array([n0[1], -n0[0]])    
#    der_1 = np.array([n1[1], -n1[0]])    
#    #b = der_0 / 3. + a
#    #c = d - der_1 / 3.
#    VE = 0.0001
#    ad = np.linalg.norm(a-d)
#    #length_estimator = lambda x: (ad + \
#    #                             (((a[0]-x[0])**2 + (a[1]-x[1])**2)**0.5 + \
#    #                              ((x[2]-x[0])**2 + (x[3]-x[1])**2)**0.5 + \
#    #                              ((d[0]-x[2])**2 + (d[1]-x[3])**2)**0.5))/2.0
#    #length_estimator = lambda x: BezierArcLength(a, np.array([x[0], x[1]]), np.array([x[2],x[3]]), d )
#    #length_constraint_d0x = lambda x: VE - np.abs( x[0]/((x[0]**2+x[1]**2)**0.5) - der_0[0]) 
#    #length_constraint_d0y = lambda x: VE - np.abs( x[1]/((x[0]**2+x[1]**2)**0.5) - der_0[1])
#    #length_constraint_d1x = lambda x: VE - np.abs( x[2]/((x[2]**2+x[3]**2)**0.5) - der_1[0]) 
#    #length_constraint_d1y = lambda x: VE - np.abs( x[3]/((x[2]**2+x[3]**2)**0.5) - der_1[1])
                                         
#    #opt_constr = ({'type': 'eq', 'fun': length_constraint_d0x},
#    #              {'type': 'eq', 'fun': length_constraint_d0y},
#    #              {'type': 'eq', 'fun': length_constraint_d1x},
#    #              {'type': 'eq', 'fun': length_constraint_d1y})

#    opt_constr = ({'type': 'ineq', 'fun': lambda x: x[0] + VE },
#                  {'type': 'ineq', 'fun': lambda x: x[1] + VE })

#    length_estimator = lambda x: BezierArcLength(a, x[0]*der_0 + a, d - x[1]*der_1, d ) 
#    #length_estimator = lambda x: (ad + x[0] + np.linalg.norm( (a+x[0]*der_0) - (d-x[1]*der_1)) +  x[1])/2.0

#    r = optimize.minimize(length_estimator, x0=[0., 0.], method= 'SLSQP', constraints=opt_constr) 
#    #method= 'COBYLA'
#    l0 = r.x[0]
#    l1 = r.x[1]
#    b = a + l0 * der_0
#    c = d - l1 * der_1
#    #b = np.array([r.x[0], r.x[1]])
#    #c = np.array([r.x[2], r.x[3]])
#    crv_length = BezierArcLength(a,b,c,d)
#    #est_crv_length = estimate_bezier_length(a,b,c,d)

#    part_length = crv_length * t0
#    res_t = get_t_for_length(a, b, c, d, part_length)
#    res_pt, res_der = eval_cubic_bezier(a, b, c, d, res_t)
#    res_norm = np.array([-res_der[1], res_der[0]])
#    res_norm /= np.linalg.norm(res_norm)
    
#    DEBUG = True
#    #DEBUG = False
#    if DEBUG:
#        verts = [
#            (a[0], a[1]), # P0
#            (b[0], b[1]), # P1
#            (c[0], c[1]), # P2
#            (d[0], d[1]), # P3
#            ]
        
#        codes = [Path.Path.MOVETO,
#                 Path.Path.CURVE4,
#                 Path.Path.CURVE4,
#                 Path.Path.CURVE4,
#                 ]
        
#        crv = Path.Path(verts, codes)
        
#        fig = plt.figure()
#        #ax = fig.add_subplot(111)
#        nr0 = plt.Arrow( p0[0], p0[1], n0[0]*.1, n0[1]*.1, width=0.03, fc='b', ec='none' )
#        nr1 = plt.Arrow( p1[0], p1[1], n1[0]*.1, n1[1]*.1, width=0.03, fc='b', ec='none' )
#        dr0 = plt.Arrow( p0[0], p0[1], der_0[0]*.1, der_0[1]*.1, width=0.03, fc='g', ec='none' )
#        dr1 = plt.Arrow( p1[0], p1[1], der_1[0]*.1, der_1[1]*.1, width=0.03, fc='g', ec='none' )
#        nr2 = plt.Arrow( res_pt[0], res_pt[1], res_norm[0]*.1, res_norm[1]*.1, 
#                        width=0.05, fc='r', ec='none' )
#        plt.gca().add_patch(nr0)
#        plt.gca().add_patch(nr1)
#        plt.gca().add_patch(nr2)
#        plt.gca().add_patch(dr0)
#        plt.gca().add_patch(dr1)
#        #plt.gca().add_patch(crv)
#        patch = patches.PathPatch(crv, facecolor='none', lw=2)
#        plt.gca().add_patch(patch)
#        plt.show()


#    return res_pt, res_norm

#============================== END OF FILE ==================================