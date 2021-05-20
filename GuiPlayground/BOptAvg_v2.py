import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from CSubd2D import *

def bspline_average_export( t0, t1, b_open, p0, p1, n0, n1 ):
    p, n = bspline_average(t0, p0, p1, n0, n1)
    return p, n, p, 0,0,0

#-----------------------------------------------------------------------------
def bspline_average_v1(t0, p0, p1, n0, n1):
    ca_pt, ca_norm, ca_center, ca_radius, ca_beta0, ca_beta1 = \
      circle_avg(t0, 1. - t0, True, p0, p1, n0, n1)

    gs = []
    n_pts = 5
    angle_inc = (np.abs(ca_beta1 - ca_beta0))/ (n_pts - 1)
    if ca_beta1 > ca_beta0:
        angle_koef = 1. 
    else:
        if ca_beta0 - ca_beta1 > 180.:
            angle_koef = 1.
            angle_inc = (360. - np.abs(ca_beta1 - ca_beta0))/ (n_pts - 1)
        else:
            angle_koef = -1.
    for i in range(n_pts):
        gs.append(ca_beta0 + angle_koef * i * angle_inc)
    ctrl_pts = []
    for g in gs:
        g = g * np.pi / 180.0
        ctrl_pts.append(ca_center + np.array([ np.cos(g), np.sin(g)])*ca_radius)

    degree = 3
    ctrl_pts = np.array(ctrl_pts)
    n_ctrl_pts = len(ctrl_pts)
    ctrl_pts_x = ctrl_pts[:,0]
    ctrl_pts_y = ctrl_pts[:,1]

    knots = range(len(ctrl_pts_x))
    bspl_x = si.make_interp_spline(knots,ctrl_pts_x,degree)
    bspl_y = si.make_interp_spline(knots,ctrl_pts_y,degree)

    ctrl_pts_x, ctrl_pts_y = bspl_x.tck[1], bspl_y.tck[1]

    der0_dir = -angle_koef * np.array([n0[1], -n0[0]])    
    der1_dir = -angle_koef * np.array([n1[1], -n1[0]])    

    der0_len = ((ctrl_pts_x[0] - ctrl_pts_x[1])**2 + \
                (ctrl_pts_y[0] - ctrl_pts_y[1])**2 ) **0.5
    ctrl_pt_2 = p0 + der0_dir * der0_len
    bspl_x.c[1] = ctrl_pt_2[0]
    bspl_y.c[1] = ctrl_pt_2[1]

    der1_len = ((ctrl_pts_x[-1] - ctrl_pts_x[-2])**2 + \
                (ctrl_pts_y[-1] - ctrl_pts_y[-2])**2 ) **0.5
    ctrl_pt_m2 = p1 - der1_dir * der1_len
    bspl_x.c[-2] = ctrl_pt_m2[0]
    bspl_y.c[-2] = ctrl_pt_m2[1]

    res_t = (knots[-1] - knots[0]) * t0
    res_pt  = np.array([bspl_x(res_t), bspl_y(res_t)])
    res_norm = np.array([bspl_y.derivative()(res_t), -bspl_x.derivative()(res_t)] )
    res_norm /= np.linalg.norm(res_norm)
    if np.linalg.norm(res_norm - ca_norm) > np.linalg.norm(res_norm + ca_norm):
        res_norm = -res_norm

    #Experiment
    res_norm = ca_norm
    

    #res_norm_0 = np.array([bspl_y.derivative()(0), -bspl_x.derivative()(0)] )
    #res_norm_0 /= np.linalg.norm(res_norm_0)
    #res_norm_K = np.array([bspl_y.derivative()(knots[-1]), -bspl_x.derivative()(knots[-1])] )
    #res_norm_K /= np.linalg.norm(res_norm_K)
    #if not vec_eeq(res_norm_0, n0):
    #    res_norm = -res_norm

    #DEBUG = True
    DEBUG = False
    if DEBUG:
        ipl_t = np.linspace(knots[0], knots[-1], 1000)
        x_i, y_i = bspl_x(ipl_t), bspl_y(ipl_t)
        fig = plt.figure()
        plt.plot(ctrl_pts_x, ctrl_pts_y, '-og')
        plt.plot(x_i, y_i, 'r')
        plt.xlim([min(ctrl_pts_x) - 0.3, max(ctrl_pts_x) + 0.3])
        plt.ylim([min(ctrl_pts_y) - 0.3, max(ctrl_pts_y) + 0.3])
        cr1 = plt.Circle( (ca_center[0], ca_center[1]), radius=ca_radius, fc='y', ec='none')
        dn0 = plt.Arrow( p0[0], p0[1], n0[0]*.1, n0[1]*.1, width=0.03, fc='b', ec='none' )
        dn1 = plt.Arrow( p1[0], p1[1], n1[0]*.1, n1[1]*.1, width=0.03, fc='b', ec='none' )
        gpR = plt.Circle( (res_pt[0], res_pt[1]), radius=0.02, fc='r', ec='none')
        gnR = plt.Arrow( res_pt[0], res_pt[1], res_norm[0]*.1, res_norm[1]*.1, width=0.03, fc='r', ec='none' )

        plt.gca().add_patch(cr1)
        plt.gca().add_patch(dn0)
        plt.gca().add_patch(dn1)
        plt.gca().add_patch(gpR)
        plt.gca().add_patch(gnR)

        plt.axis('equal')
        plt.show()
    return res_pt, res_norm
#-----------------------------------------------------------------------------
def bspline_average_v2(t0, p0, p1, n0, n1):
    ca_pt, ca_norm, ca_center, ca_radius, ca_beta0, ca_beta1, ca_betaR = \
      circle_avg(t0, 1. - t0, True, p0, p1, n0, n1)

    gs = []
    d_beta0_R = np.abs(ca_beta0 - ca_betaR)
    d_beta1_R = np.abs(ca_beta1 - ca_betaR)
    b_betaR_closer_beta0 = d_beta1_R > d_beta0_R
    b_extrapol = False
    n_extrap_pts = int(2 + (np.abs(1.0 / t0) if t0 < 0.0 else 1.0 / (t0-1.)))
    n_pts = 0
    angle_inc = 0.0
    angle_src = 0.0
    angle_dst = 0.0
    if ca_beta1 > ca_beta0:
        if ca_beta1 - ca_beta0 > 180.:
            # zero inside the interval
            b_extrapol = ca_beta0 < ca_betaR < ca_beta1 
            if b_extrapol:
                angle_src = ca_beta1 if b_betaR_closer_beta0 else ca_betaR
                angle_dst = ca_betaR if b_betaR_closer_beta0 else ca_beta0
                angle_inc = d_beta0_R if b_betaR_closer_beta0 else d_beta1_R
                n_pts = n_extrap_pts 
            else:
                angle_src = ca_beta1
                angle_dst = ca_beta0
                n_pts = 6
                angle_inc = (360. - np.abs(ca_beta1 - ca_beta0))/ (n_pts - 1)
            t0 = 1. - t0
        else:
            b_extrapol = not ca_beta0 < ca_betaR < ca_beta1 
            if b_extrapol:
                angle_src = ca_betaR if b_betaR_closer_beta0 else ca_beta0
                angle_dst = ca_betaR if not b_betaR_closer_beta0 else ca_beta1
                angle_inc = d_beta0_R if b_betaR_closer_beta0 else d_beta1_R
                n_pts = n_extrap_pts 
            else:
                angle_src = ca_beta0
                angle_dst = ca_beta1
                n_pts = 6
                angle_inc = (np.abs(ca_beta1 - ca_beta0))/ (n_pts - 1)

    else:
        if ca_beta0 - ca_beta1 > 180.:
            # zero inside the interval
            b_extrapol = ca_beta1 < ca_betaR < ca_beta0 
            if b_extrapol:
                angle_src = ca_betaR if b_betaR_closer_beta0 else ca_beta0
                angle_dst = ca_beta1 if b_betaR_closer_beta0 else ca_betaR
                angle_inc = d_beta0_R if b_betaR_closer_beta0 else d_beta1_R
                n_pts = n_extrap_pts 
            else:
                angle_src = ca_beta0
                angle_dst = ca_beta1
                n_pts = 6
                angle_inc = (360. - np.abs(ca_beta0 - ca_beta1))/ (n_pts - 1)
        else:
            b_extrapol = not ca_beta1 < ca_betaR < ca_beta0 
            if b_extrapol:
                angle_src = ca_betaR if not b_betaR_closer_beta0 else ca_beta1
                angle_dst = ca_betaR if b_betaR_closer_beta0 else ca_beta0
                angle_inc = d_beta0_R if b_betaR_closer_beta0 else d_beta1_R
                n_pts = n_extrap_pts 
            else:
                angle_src = ca_beta1
                angle_dst = ca_beta0
                n_pts = 6
                angle_inc = (np.abs(ca_beta1 - ca_beta0))/ (n_pts - 1)
            t0 = 1. - t0

    for i in range(n_pts):
        gs.append((angle_src + i * angle_inc)%360)
    ctrl_pts = []
    for g in gs:
        g = g * np.pi / 180.0
        ctrl_pts.append(ca_center + np.array([ np.cos(g), np.sin(g)])*ca_radius)

    degree = 3
    ctrl_pts = np.array(ctrl_pts)
    n_ctrl_pts = len(ctrl_pts)
    ctrl_pts_x = ctrl_pts[:,0]
    ctrl_pts_y = ctrl_pts[:,1]

    knots = range(len(ctrl_pts_x))
    bspl_x = si.make_interp_spline(knots,ctrl_pts_x,degree)
    bspl_y = si.make_interp_spline(knots,ctrl_pts_y,degree)

    ctrl_pts_x, ctrl_pts_y = bspl_x.tck[1], bspl_y.tck[1]

    der0_dir = -np.array([n0[1], -n0[0]])    
    der1_dir = -np.array([n1[1], -n1[0]])    

    anchor_idx_src = 1 if b_extrapol and b_betaR_closer_beta0 else 0
    anchor_idx_dst = -2 if b_extrapol and not b_betaR_closer_beta0 else -1

    #der0_len = ((ctrl_pts_x[0] - ctrl_pts_x[1])**2 + \
    #            (ctrl_pts_y[0] - ctrl_pts_y[1])**2 ) **0.5
    #der0_len = np.linalg.norm(ctrl_pts[0] - ctrl_pts[1])
    #if b_extrapol and b_betaR_closer_beta0:
    #    anchor_pt = np.array(bspl_x.c[1], bspl_y.c[1])
    #    ctrl_pt_0 = anchor_pt + der0_dir * der0_len
    #    ctrl_pt_2 = anchor_pt - der0_dir * der0_len
    #    bspl_x.c[2] = ctrl_pt_2[0]
    #    bspl_y.c[2] = ctrl_pt_2[1]
    #    bspl_x.c[0] = ctrl_pt_0[0]
    #    bspl_y.c[0] = ctrl_pt_0[1]
    #else:
    #    anchor_pt = np.array(bspl_x.c[0], bspl_y.c[0])
    #    ctrl_pt_1 = anchor_pt + der0_dir * der0_len
    #    bspl_x.c[1] = ctrl_pt_1[0]
    #    bspl_y.c[1] = ctrl_pt_1[1]

    #bspl_x.c[n_ctrl_pts/2] = ctrl_pts[n_ctrl_pts/2][0]
    #bspl_y.c[n_ctrl_pts/2] = ctrl_pts[n_ctrl_pts/2][1]

    ##der1_len = ((ctrl_pts_x[-1] - ctrl_pts_x[-2])**2 + \
    ##            (ctrl_pts_y[-1] - ctrl_pts_y[-2])**2 ) **0.5
    #der1_len = np.linalg.norm(ctrl_pts[-2] - ctrl_pts[-1])
    #ctrl_pt_m2 = p1 - der1_dir * der1_len
    #bspl_x.c[-2] = ctrl_pt_m2[0]
    #bspl_y.c[-2] = ctrl_pt_m2[1]

    knots_dist = knots[-1] - knots[0]
    res_t = knots_dist * t0
    res_pt  = np.array([bspl_x(res_t), bspl_y(res_t)])
    res_norm = np.array([bspl_y.derivative()(res_t), -bspl_x.derivative()(res_t)] )
    res_norm /= np.linalg.norm(res_norm)
    if np.linalg.norm(res_norm - ca_norm) > np.linalg.norm(res_norm + ca_norm):
        res_norm = -res_norm

    #Experiment
    #res_norm = ca_norm
    

    #res_norm_0 = np.array([bspl_y.derivative()(0), -bspl_x.derivative()(0)] )
    #res_norm_0 /= np.linalg.norm(res_norm_0)
    #res_norm_K = np.array([bspl_y.derivative()(knots[-1]), -bspl_x.derivative()(knots[-1])] )
    #res_norm_K /= np.linalg.norm(res_norm_K)
    #if not vec_eeq(res_norm_0, n0):
    #    res_norm = -res_norm

    #DEBUG = True
    DEBUG = 'IN_DEBUG' in globals()
    if DEBUG:
        param0 = min(knots[0], t0*knots_dist)
        param1 = max(knots[-1], t0*knots_dist)
        ipl_t = np.linspace(param0, param1, 1000)
        x_i, y_i = bspl_x(ipl_t), bspl_y(ipl_t)
        fig = plt.figure()
        plt.plot(ctrl_pts_x, ctrl_pts_y, '-og')
        plt.plot(x_i, y_i, 'r')
        plt.xlim([min(ctrl_pts_x) - 0.3, max(ctrl_pts_x) + 0.3])
        plt.ylim([min(ctrl_pts_y) - 0.3, max(ctrl_pts_y) + 0.3])
        cr1 = plt.Circle( (ca_center[0], ca_center[1]), radius=ca_radius, fc='y', ec='none')
        dn0 = plt.Arrow( p0[0], p0[1], n0[0]*.1, n0[1]*.1, width=0.03, fc='b', ec='none' )
        dn1 = plt.Arrow( p1[0], p1[1], n1[0]*.1, n1[1]*.1, width=0.03, fc='b', ec='none' )
        gpR = plt.Circle( (res_pt[0], res_pt[1]), radius=0.02, fc='r', ec='none')
        gnR = plt.Arrow( res_pt[0], res_pt[1], res_norm[0]*.1, res_norm[1]*.1, width=0.03, fc='r', ec='none' )

        plt.gca().add_patch(cr1)
        plt.gca().add_patch(dn0)
        plt.gca().add_patch(dn1)
        plt.gca().add_patch(gpR)
        plt.gca().add_patch(gnR)

        plt.axis('equal')
        plt.show()
    return res_pt, res_norm
 
#-----------------------------------------------------------------------------
def bspline_average(t0, p0, p1, n0, n1):
    orig_t0 = t0
    ca_pt, ca_norm, ca_center, ca_radius, ca_beta0, ca_beta1 = \
      circle_avg(t0, 1. - t0, True, p0, p1, n0, n1)

    theta = np.abs(ca_beta0 - ca_beta1)
    compl_theta = 360. - theta
    theta_ratio  = theta / compl_theta

    gs = []
    n_total_pts = 24
    n_int_pts = max(4, int(n_total_pts * theta_ratio))
    angle_inc = 0.0
    angle_src = 0.0
    angle_dst = 0.0
    b_passes_zero = False
    if ca_beta1 > ca_beta0:
        if ca_beta1 - ca_beta0 > 180.:
            # zero inside the interval
            b_passes_zero = True
            angle_src = ca_beta1
            angle_dst = ca_beta0
            angle_inc = (360. - np.abs(ca_beta1 - ca_beta0))/ (n_int_pts - 1)
            t0 = 1. - t0
        else:
            angle_src = ca_beta0
            angle_dst = ca_beta1
            angle_inc = (np.abs(ca_beta1 - ca_beta0))/ (n_int_pts - 1)

    else:
        if ca_beta0 - ca_beta1 > 180.:
            # zero inside the interval
            b_passes_zero = True
            angle_src = ca_beta0
            angle_dst = ca_beta1
            angle_inc = (360. - np.abs(ca_beta0 - ca_beta1))/ (n_int_pts - 1)
        else:
            angle_src = ca_beta1
            angle_dst = ca_beta0
            angle_inc = (np.abs(ca_beta1 - ca_beta0))/ (n_int_pts - 1)
            t0 = 1. - t0

    for i in range(n_int_pts):
        gs.append((angle_src + i * angle_inc)%360)
    compl_num_pts = n_total_pts - n_int_pts + 2 # * int( (360. / np.abs(angle_dst - angle_src) ))
    #compl_num_pts = n_pts
    if b_passes_zero:
        compl_angle_inc = np.abs(angle_dst - angle_src) / (compl_num_pts - 1)
    else:
        compl_angle_inc = (360. - np.abs(angle_dst - angle_src)) / (compl_num_pts - 1)
    for i in range(1, compl_num_pts-1):
        gs.append((angle_dst + i * compl_angle_inc)%360)

    ctrl_pts = []
    for g in gs:
        g = g * np.pi / 180.0
        ctrl_pts.append(ca_center + np.array([ np.cos(g), np.sin(g)])*ca_radius)


    #fig = plt.figure()
    #ctrl_pts = np.array(ctrl_pts)
    #ctrl_pts_x = ctrl_pts[:,0]
    #ctrl_pts_y = ctrl_pts[:,1]
    #plt.plot(ctrl_pts_x, ctrl_pts_y, '-og')
    #plt.xlim([min(ctrl_pts_x) - 0.3, max(ctrl_pts_x) + 0.3])
    #plt.ylim([min(ctrl_pts_y) - 0.3, max(ctrl_pts_y) + 0.3])
    #cr1 = plt.Circle( (ca_center[0], ca_center[1]), radius=ca_radius, fc='y', ec='none')
    #dn0 = plt.Arrow( p0[0], p0[1], n0[0]*.1, n0[1]*.1, width=0.03, fc='b', ec='none' )
    #dn1 = plt.Arrow( p1[0], p1[1], n1[0]*.1, n1[1]*.1, width=0.03, fc='b', ec='none' )
    #plt.gca().add_patch(cr1)
    #plt.gca().add_patch(dn0)
    #plt.gca().add_patch(dn1)
    #plt.axis('equal')
    #plt.show()

    n_points = len(ctrl_pts)
    der0_dir = np.array([n0[1], -n0[0]])    
    der1_dir = np.array([n1[1], -n1[0]])    
    z2o_ccw = vec_eeq(p0, ctrl_pts[0])
    if not z2o_ccw:
        der0_dir, der1_dir = der1_dir, der0_dir
    i0 = 0
    i0m1 = (i0-1+n_points)%n_points
    i0p1 = (i0+1)%n_points
    ctrl_pts[i0p1] = ctrl_pts[i0] - der0_dir * np.linalg.norm(ctrl_pts[i0p1] - ctrl_pts[i0])
    ctrl_pts[i0m1] = ctrl_pts[i0] + der0_dir * np.linalg.norm(ctrl_pts[i0m1] - ctrl_pts[i0])
    i1 = n_int_pts - 1
    i1m1 = (i1-1+n_points)%n_points
    i1p1 = (i1+1)%n_points
    ctrl_pts[i1p1] = ctrl_pts[i1] - der1_dir * np.linalg.norm(ctrl_pts[i1p1] - ctrl_pts[i1])
    ctrl_pts[i1m1] = ctrl_pts[i1] + der1_dir * np.linalg.norm(ctrl_pts[i1m1] - ctrl_pts[i1])

    degree = 3
    ctrl_pts = ctrl_pts + ctrl_pts[0:degree + 1]
    ctrl_pts = np.array(ctrl_pts)
    n_points = len(ctrl_pts)
    ctrl_pts_x = ctrl_pts[:,0]
    ctrl_pts_y = ctrl_pts[:,1]

    knots = range(len(ctrl_pts_x))
    ipl_t = np.linspace(1.0, len(ctrl_pts) - degree, 1000)

    x_tup = si.splrep(knots, ctrl_pts_x, k=degree, per=1)
    y_tup = si.splrep(knots, ctrl_pts_y, k=degree, per=1)
    x_list = list(x_tup)
    xl = ctrl_pts_x.tolist()
    x_list[1] = [0.0] + xl + [0.0, 0.0, 0.0, 0.0]

    y_list = list(y_tup)
    yl = ctrl_pts_y.tolist()
    y_list[1] = [0.0] + yl + [0.0, 0.0, 0.0, 0.0]

    x_i = si.splev(ipl_t, x_list)
    y_i = si.splev(ipl_t, y_list)

    # experiment
    #bspl_x = si.make_interp_spline(knots,ctrl_pts_x,degree)
    #bspl_y = si.make_interp_spline(knots,ctrl_pts_y,degree)
    #bspl_x = si.BSpline(x_list)
    #bspl_y = si.BSpline(ipl_t, ctrl_pts_y, degree)
    #x_i, y_i = bspl_x(ipl_t), bspl_y(ipl_t)
    # --- 

    knots_dist = knots[-1] - knots[0]
    res_t = knots_dist * t0

    knots_dist = knots[n_int_pts-1] - knots[0]
    ut0 = t0 if 0. <= t0 <= 1. else (1. - t0)
    res_t = knots_dist * ut0
    res_t = knots[n_total_pts] + res_t if res_t < 0. else res_t
    #res_t = knots[-1] - (knots[-1] - knots[n_pts])*ut0
    res_x = si.splev([res_t], x_list)
    der_x = si.splev([res_t], x_list, der=1)
    res_y = si.splev([res_t], y_list)
    der_y = si.splev([res_t], y_list, der=1)

    #res_pt  = np.array([bspl_x(res_t), bspl_y(res_t)])
    #res_norm = np.array([bspl_y.derivative()(res_t), -bspl_x.derivative()(res_t)] )
    res_pt  = np.array([res_x, res_y])
    res_norm = np.array([der_y, -der_x] )    
    res_norm /= np.linalg.norm(res_norm)
    #if np.linalg.norm(res_norm - ca_norm) > np.linalg.norm(res_norm + ca_norm):
    #    res_norm = -res_norm

    #fig = plt.figure()
    #plt.plot(x, y, '-og')
    #plt.plot(x_i, y_i, 'r')
    #plt.xlim([min(x) - 0.3, max(x) + 0.3])
    #plt.ylim([min(y) - 0.3, max(y) + 0.3])    
    #plt.title('Splined f(x(t), y(t))')
    #cr1 = plt.Circle( (ca_center[0], ca_center[1]), radius=ca_radius, fc='y', ec='none')
    #plt.gca().add_patch(cr1)
    #plt.axis('equal')
    #plt.show()

    #DEBUG = True
    DEBUG = 'IN_DEBUG' in globals()
    if DEBUG:
        param0 = min(knots[0], t0*knots_dist)
        param1 = max(knots[-1], t0*knots_dist)
        ipl_t = np.linspace(param0, param1, 1000)
        fig = plt.figure()
        plt.plot(ctrl_pts_x, ctrl_pts_y, '-oc')
        plt.plot(x_i, y_i, 'r')
        plt.xlim([min(ctrl_pts_x) - 0.3, max(ctrl_pts_x) + 0.3])
        plt.ylim([min(ctrl_pts_y) - 0.3, max(ctrl_pts_y) + 0.3])
        cr1 = plt.Circle( (ca_center[0], ca_center[1]), radius=ca_radius, fc='y', ec='none')
        dn0 = plt.Arrow( p0[0], p0[1], n0[0]*.15, n0[1]*.15, width=0.05, fc='b', ec='none' )
        dn1 = plt.Arrow( p1[0], p1[1], n1[0]*.1, n1[1]*.1, width=0.03, fc='g', ec='none' )
        gpR = plt.Circle( (res_pt[0], res_pt[1]), radius=0.02, fc='r', ec='none')
        gnR = plt.Arrow( res_pt[0], res_pt[1], res_norm[0]*.1, res_norm[1]*.1, width=0.03, fc='r', ec='none' )
        plt.title('t = ' + str(orig_t0) )
        plt.gca().add_patch(cr1)
        plt.gca().add_patch(dn0)
        plt.gca().add_patch(dn1)
        plt.gca().add_patch(gpR)
        plt.gca().add_patch(gnR)

        plt.axis('equal')
        plt.show()
    return res_pt, res_norm

#-----------------------------------------------------------------------------
if __name__ == "__main__":
    global IN_DEBUG
    IN_DEBUG = True
    points  = [(  1.0, 0.0),
               (  0.0, 1.0)]
    normals = [(  1.0, 0.0),
               (  0.0,  1.0)]
    normals[0] /= np.linalg.norm(normals[0])
    normals[1] /= np.linalg.norm(normals[1])
    p0 = np.array([points[0][0], points[0][1]])
    p1 = np.array([points[1][0], points[1][1]])
    n0 = np.array([normals[0][0], normals[0][1]])
    n1 = np.array([normals[1][0], normals[1][1]])
    bspline_average(-0.125, p0, p1, n0, n1)
