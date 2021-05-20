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
    if 0.0 <= t0 <= 1.0:
        return bspline_average_core(t0, p0, p1, n0, n1, True,
                                    ca_pt, ca_norm, 
                                    ca_center, ca_radius, 
                                    ca_beta0, ca_beta1)
    elif 0.0 > t0:
        return bspline_average_core(np.abs(t0), p1, p0, n1, n0, False,
                                    ca_pt, ca_norm, 
                                    ca_center, ca_radius, 
                                    ca_beta1, ca_beta0)
    elif 1.0 < t0:
        return bspline_average_core(t0 - 1., p0, p1, n0, n1, False,
                                    ca_pt, ca_norm, 
                                    ca_center, ca_radius, 
                                    ca_beta0, ca_beta1)
            
#-----------------------------------------------------------------------------
def bspline_average_core(t0, p0, p1, n0, n1, b_short_arc,
                         ca_pt, ca_norm, 
                         ca_center, ca_radius, 
                         ca_beta0, ca_beta1):
    theta = np.abs(ca_beta0 - ca_beta1)
    if theta > 180.:
        theta = 360. - theta 
    if not b_short_arc: 
        theta = 360. - theta

    n_interp_pts = 7 if b_short_arc else 9
    angle_inc = theta / (n_interp_pts - 1)
    angle_src = 0.0
    angle_dst = 0.0
    b_passes_zero = False
    if ca_beta1 > ca_beta0:
        if ca_beta1 - ca_beta0 > 180.:
            # zero inside the interval
            b_passes_zero = True
            angle_src = ca_beta1
            angle_dst = ca_beta0
            t0 = 1. - t0
        else:
            angle_src = ca_beta0
            angle_dst = ca_beta1

    else:
        if ca_beta0 - ca_beta1 > 180.:
            # zero inside the interval
            b_passes_zero = True
            angle_src = ca_beta0
            angle_dst = ca_beta1
        else:
            angle_src = ca_beta1
            angle_dst = ca_beta0
            t0 = 1. - t0
    if not b_short_arc:
        angle_src, angle_dst = angle_dst, angle_src
        t0 = 1. - t0

    knots = []
    for i in range(n_interp_pts):
        knots.append(angle_src + i * angle_inc)

    interp_pts = []
    for k in knots:
        k = k * np.pi / 180.0
        interp_pts.append(ca_center + \
                          np.array([ np.cos(k), np.sin(k)]) * \
                          ca_radius)

    degree = 3
    interp_pts = np.array(interp_pts)
    interp_pts_x = interp_pts[:,0]
    interp_pts_y = interp_pts[:,1]

    bspl_x = si.make_interp_spline(knots, interp_pts_x, degree)
    bspl_y = si.make_interp_spline(knots, interp_pts_y, degree)

    bspl_p0 = np.array([bspl_x.c[0], bspl_y.c[0]])
    if vec_eeq(bspl_p0, p0):
        der0_dir = -np.array([n0[1], -n0[0]])    
        der1_dir = np.array([n1[1], -n1[0]])   
    else:
        der1_dir = np.array([n0[1], -n0[0]])    
        der0_dir = -np.array([n1[1], -n1[0]])   
 
    der0_len = ((bspl_x.c[ 0] - bspl_x.c[ 1])**2 + (bspl_y.c[ 0] - bspl_y.c[ 1])**2)**0.5
    der1_len = ((bspl_x.c[-1] - bspl_x.c[-2])**2 + (bspl_y.c[-1] - bspl_y.c[-2])**2)**0.5

    bspl_x.c[ 1] = bspl_x.c[0]  + der0_dir[0] * der0_len
    bspl_y.c[ 1] = bspl_y.c[0]  + der0_dir[1] * der0_len
    bspl_x.c[-2] = bspl_x.c[-1] + der1_dir[0] * der1_len
    bspl_y.c[-2] = bspl_y.c[-1] + der1_dir[1] * der1_len
    ipl_t = np.linspace(knots[0], knots[-1], 1000)
    x_i, y_i = bspl_x(ipl_t), bspl_y(ipl_t)


    if not b_short_arc: 
        if t0 < 0.5:
            res_t = knots[0] + (360. - theta) * t0
        else:
            res_t = knots[-1] - (360. - theta) * (1. - t0)
    else:
        res_t = knots[0] + theta * t0
    
    res_pt  = np.array([bspl_x(res_t), bspl_y(res_t)])
    res_norm = np.array([bspl_y.derivative()(res_t), -bspl_x.derivative()(res_t)] )
    res_norm /= np.linalg.norm(res_norm)
    
    if np.linalg.norm(res_norm - ca_norm) > np.linalg.norm(res_norm + ca_norm):
        res_norm = -res_norm

    #DEBUG = True
    DEBUG = 'IN_DEBUG' in globals()
    if DEBUG and IN_DEBUG:
        #param0 = min(knots[0], t0*knots_dist)
        #param1 = max(knots[-1], t0*knots_dist)
        #ipl_t = np.linspace(param0, param1, 1000)
        fig = plt.figure()
        plt.plot(bspl_x.c, bspl_y.c, '-oc')
        plt.plot(interp_pts_x, interp_pts_y, '-g')
        #plt.plot(tck[1][0], tck[1][1], '-oc')
        plt.plot(x_i, y_i, 'r')
        plt.xlim([min(interp_pts_x) - 0.3, max(interp_pts_x) + 0.3])
        plt.ylim([min(interp_pts_y) - 0.3, max(interp_pts_y) + 0.3])
        cr1 = plt.Circle( (ca_center[0], ca_center[1]), radius=ca_radius, fc='y', ec='none')
        dn0 = plt.Arrow( p0[0], p0[1], n0[0]*.15, n0[1]*.15, width=0.05, fc='b', ec='none' )
        dn1 = plt.Arrow( p1[0], p1[1], n1[0]*.1, n1[1]*.1, width=0.03, fc='g', ec='none' )
        gpR = plt.Circle( (res_pt[0], res_pt[1]), radius=0.02, fc='r', ec='none')
        gnR = plt.Arrow( res_pt[0], res_pt[1], res_norm[0]*.1, res_norm[1]*.1, width=0.03, fc='r', ec='none' )
        plt.title('t = ' + str(t0) )
        plt.gca().add_patch(cr1)
        plt.gca().add_patch(dn0)
        plt.gca().add_patch(dn1)
        plt.gca().add_patch(gpR)
        plt.gca().add_patch(gnR)

        plt.axis('equal')
        plt.show()
    return res_pt, res_norm

#-----------------------------------------------------------------------------
def create_input_on_a_polygon5():
    c30 = 3.0**0.3/2.0
    pts = [np.array([-2, -2]),
           np.array([ 4, -3]),
           np.array([ 2,  0]),
           np.array([ 5,  3.5]),
           np.array([-1,  3])]
    nrm = [np.array([-1, -1]),
           np.array([ 1, -1]),
           np.array([ 1,  0]),
           #np.array([ -1,  0]),
           #np.array([  0,  1]),
           #np.array([  0,  -1]),

           #np.array([  c30, -0.5]), #rd
           #np.array([ c30, 0.5]), #ru
           #np.array([ -c30, 0.5]), #lu
           #np.array([ -c30, -0.5]), #ld


           np.array([ 1,  1]),
           np.array([-1,  1])]
    nnrm = []
    for n in nrm:
        nnrm.append( n/np.linalg.norm(n))
    return pts, nnrm
#-----------------------------------------------------------------------------
def create_input_on_a_square():
    pts = [np.array([ 0., 0.]),
           np.array([ 10., 0.]),
           np.array([ 10., 10.]),
           np.array([ 0., 10.])]
    nrm = [np.array([  0., -1.]),
           np.array([  1.,  0.]),
           np.array([  0.,  1.]),
           np.array([ -1.,  0.])]
    return pts, nrm

#-----------------------------------------------------------------------------
def create_input_on_a_square_diff_norms():
    pts = [np.array([ 0., 0.]),
           np.array([ 1., 0.]),
           np.array([ 1., 1.]),
           np.array([ 0., 1.])]
    nrm = [np.array([  0., -1.]),
           np.array([  1.,  1.]),
           np.array([  0.,  -1.]),
           np.array([  1.,  -1.])]
    nnrm = []
    for n in nrm:
        nnrm.append( n/np.linalg.norm(n))
    return pts, nnrm

#-----------------------------------------------------------------------------
def plot_pts_and_norms(pts, nrm, b_open, draw_norms, clr, 
                       linestyle='', linewidth=1.0, cnvs = plt ):
    n = len(pts)
    nn = n-1 if b_open else n
    for i in range(nn):
        curr_pt = pts[i]
        next_pt = pts[(i+1)%n]
        if linestyle.startswith('da'):
            cnvs.plot([curr_pt[0], next_pt[0]], 
                        [curr_pt[1], next_pt[1]], 
                        color=clr, linestyle=linestyle,
                        linewidth=linewidth, dashes=(1,15))
        else:
            cnvs.plot([curr_pt[0], next_pt[0]], 
                        [curr_pt[1], next_pt[1]], 
                        color=clr, linestyle=linestyle, 
                        linewidth=linewidth)

        if draw_norms:
            curr_norm = nrm[i]
            #wdth = 0.15 if i != 2 else 0.3
            wdth = 0.3
            #colr = clr if i != 2 else 'r'
            colr = clr
            #curr_norm /= 2.0 #if i == 2 else 1.0
            gnr = cnvs.Arrow(curr_pt[0], curr_pt[1], 
                               curr_norm[0], curr_norm[1], 
                               width=wdth, fc=colr, ec='none' )
            cnvs.gca().add_patch(gnr)
    if draw_norms and b_open:
        curr_norm = nrm[-1]
        curr_pt = pts[-1]
        gnr = cnvs.Arrow(curr_pt[0], curr_pt[1], 
                           curr_norm[0], curr_norm[1], 
                           width=0.05, fc=clr, ec='none' )
        cnvs.gca().add_patch(gnr)
#-----------------------------------------------------------------------------
def build_curves():
    n_of_iterations = 1
    b_open = False
    #subd_pts, subd_nrm = create_input_on_a_polygon5()
    subd_pts, subd_nrm = create_input_on_a_square()
    #subd_pts, subd_nrm = create_input_on_a_square_diff_norms()
    orig_pts = subd_pts[:]
    bsubd_4p_pts, bsubd_4p_nrm = subd_pts[:], subd_nrm[:]
    isubd_4p_pts, isubd_4p_nrm = subd_pts[:], subd_nrm[:]
    csubd_4p_pts, csubd_4p_nrm = subd_pts[:], subd_nrm[:]

    for k in range(n_of_iterations):
        bsubd_4p_pts, bsubd_4p_nrm = subd_4PT_one_step(bsubd_4p_pts, bsubd_4p_nrm, 
                                                        b_open, bspline_average_export)
        csubd_4p_pts, csubd_4p_nrm = subd_4PT_one_step(csubd_4p_pts, csubd_4p_nrm, 
                                                        b_open, circle_avg)
        isubd_4p_pts, isubd_4p_nrm = double_polygon(isubd_4p_pts, isubd_4p_nrm,
                                                    True, b_open,
                                                    bspline_average_export)

    fig = plt.figure()#figsize=(8,8), dpi=100, frameon = False)
    #frame1 = plt.gca()
    #frame1.axes.get_xaxis().set_visible(False)
    #frame1.axes.get_yaxis().set_visible(False)

    plot_pts_and_norms(orig_pts, subd_nrm, b_open, False, clr='k', linewidth=0.4, linestyle='dotted')
    plot_pts_and_norms(bsubd_4p_pts, bsubd_4p_nrm, b_open, True, clr='r', linewidth=1.0, linestyle='solid')
    plot_pts_and_norms(isubd_4p_pts, isubd_4p_nrm, b_open, True, clr='c', linewidth=1.0, linestyle='solid')
    plot_pts_and_norms(csubd_4p_pts, csubd_4p_nrm, b_open, True, clr='g', linewidth=1.0, linestyle='solid')

    plt.axis('equal')
    plt.xlim([-10, 10])
    plt.ylim([-10, 10])
    #plt.axis('off')
    plt.show()

#-----------------------------------------------------------------------------
if __name__ == "__main__":

    global IN_DEBUG
    IN_DEBUG = True
    #build_curves()
    points  = [(  0.0, 0.0),
               (  1.0, 0.0)]
    normals = [( -1.0, 1.0),
               (  1.0, 1.0)]
    normals[0] /= np.linalg.norm(normals[0])
    normals[1] /= np.linalg.norm(normals[1])
    p0 = np.array([points[0][0], points[0][1]])
    p1 = np.array([points[1][0], points[1][1]])
    n0 = np.array([normals[0][0], normals[0][1]])
    n1 = np.array([normals[1][0], normals[1][1]])
    bspline_average(1.125, p0, p1, n0, n1)
