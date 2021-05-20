import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
import matplotlib.path as Path
import matplotlib.patches as patches

from CSubd2D import *
from BezierAvg import eval_cubic_bezier
#from  BOptAvg import *
from SmoothnessAnalysis import *

#-----------------------------------------------------------------------------
def sliding_circ_average_export(t0, t1, b_open, p0, p1, n0, n1):
    q_pt = get_intersection_point(p0, p1, n0, n1)
    r0 = get_dist( q_pt, p0)
    r1 = get_dist( q_pt, p1)
    res_radius = r0*t0 + r1 * t1
    res_t = get_weighted_angle(t0, t1, n0, n1)
    res_norm = np.array([M.cos(res_t), M.sin(res_t)])
    res_pt = q_pt + res_norm * res_radius 
    return res_pt, res_norm, res_pt, 0,0,0



#-----------------------------------------------------------------------------
def cheb_nodes(N, a, b):
    N = N - 2
    jj = 2.*np.arange(N) + 1
    c = np.cos(np.pi * jj / (2 * N) )[::-1]
    x = 0.5*(b-a)*c + 0.5 * (a+b)
    x = np.append(np.insert(x, 0, a), b)
    return x

#-----------------------------------------------------------------------------
def bspline_average_export_v3( t0, t1, b_open, p0, p1, n0, n1 ):
    p, n = bspline_average_v3(t0, p0, p1, n0, n1)
    return p, n, p, 0,0,0

#-----------------------------------------------------------------------------
def extrapolate_bezier(a, b, c, d, t):
    v = 1. / (1. + 1./t)
    u = 1. - v
    A = a
    factor = 1. / u
    b_a_len = get_dist(a, b)
    b_a_uvec = (b - a)/b_a_len
    P = a + b_a_uvec * (b_a_len * factor)
    c_b_len = get_dist(b, c)
    c_b_uvec = (c - b)/c_b_len
    H = b + c_b_uvec * (c_b_len * factor)
    d_c_len = get_dist(d, c)
    d_c_uvec = (d - c)/d_c_len
    E = c + d_c_uvec * (d_c_len * factor)
    H_P_len = get_dist(H, P)
    H_P_uvec = (H - P)/H_P_len
    Q = P + H_P_uvec * (H_P_len * factor)
    E_H_len = get_dist(E, H)
    E_H_uvec = (E - H)/E_H_len
    F = H + E_H_uvec * (E_H_len * factor)
    F_Q_len = get_dist(F, Q)
    F_Q_uvec = (F-Q)/F_Q_len
    G = Q + F_Q_uvec * (F_Q_len * factor)
    return A, P, Q, G

#-----------------------------------------------------------------------------
def correct_derivative_v3(anchor_pt, norm_dir, der_length, b_ccw_rot):
    if b_ccw_rot:
        der_dir = np.array([ -norm_dir[1],  norm_dir[0] ])
    else:
        der_dir = np.array([  norm_dir[1], -norm_dir[0] ])

    res_ctrl_pt = anchor_pt + der_dir * der_length
    return res_ctrl_pt

#-----------------------------------------------------------------------------
def bspline_average_v3(t0, p0, p1, n0, n1):
    p0_p1_dist = get_dist( p0, p1 )
    if p0_p1_dist < 0.001:
        return p1, n1

    theta = get_angle_between(n0, n1)
    der_length = p0_p1_dist/(3. * (np.cos(theta/4.) ** 2.))
    a = p0 
    b = correct_derivative_v3(p0, n0, der_length, b_ccw_rot = False)
    c = correct_derivative_v3(p1, n1, der_length, b_ccw_rot = True)
    d = p1

    res_t = t0
    #if t0 < 0.:
    #    d, c, b, a = extrapolate_bezier(d, c, b, a, np.abs(t0))
    #    #a, b, c, d = extrapolate_bezier(d, c, b, a, np.abs(t0))
    #    res_t = 0.0
    #elif t0 > 1.:
    #    a, b, c, d = extrapolate_bezier(a, b, c, d, t0-1.)
    #    res_t = 1.

    res_pt, res_der = eval_cubic_bezier( a, b, c, d, 1. - res_t )
    res_norm = np.array([-res_der[1], res_der[0]])
    res_norm /= np.linalg.norm(res_norm)
    #if np.linalg.norm(res_norm - ca_norm) > \
    #   np.linalg.norm(res_norm + ca_norm):
    #    res_norm = -res_norm    

    DEBUG = 'IN_DEBUG' in globals()
    if DEBUG and IN_DEBUG:
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
        
        ca_pt, ca_norm, ca_center, ca_radius, ca_beta0, ca_beta1 = \
          circle_avg(0.5, 0.5, True, p0, p1, n0, n1)
        fig = plt.figure() 
        plt.xlim([(ca_center[0] - ca_radius)*1.2, (ca_center[0] + ca_radius)*1.2])
        plt.ylim([(ca_center[1] - ca_radius)*1.2, (ca_center[1] + ca_radius)*1.2])
        cr1 = plt.Circle( (ca_center[0], ca_center[1]), radius=ca_radius, fc='y', ec='none')
        norm_len_factor = 0.2
        nr000 = plt.Arrow(p0[0], p0[1], 
                          n0[0]*norm_len_factor, n0[1]*norm_len_factor, 
                          width=0.03, fc='b', ec='none' )
        nr_res = plt.Arrow(res_pt[0], res_pt[1], 
                          res_norm[0]*norm_len_factor, res_norm[1]*norm_len_factor, 
                          width=0.03, fc='r', ec='none' )
        nr_ca = plt.Arrow(ca_pt[0], ca_pt[1], 
                          ca_norm[0]*norm_len_factor, ca_norm[1]*norm_len_factor, 
                          width=0.03, fc='g', ec='none' )
        nr100 = plt.Arrow(p1[0], p1[1], 
                          n1[0]*norm_len_factor, n1[1]*norm_len_factor, 
                          width=0.03, fc='b', ec='none' )
        plt.plot([a[0], b[0], c[0], d[0]], 
                 [a[1], b[1], c[1], d[1]], '-c' )
        plt.gca().add_patch(cr1)
        plt.gca().add_patch(nr000)
        plt.gca().add_patch(nr_res)
        plt.gca().add_patch(nr_ca)
        plt.gca().add_patch(nr100)
        patch = patches.PathPatch(crv, facecolor='none', lw=1)
        plt.gca().add_patch(patch)
        plt.axis('equal')
        #plt.axis([-1.0, 2., -1.0, 2])
        plt.show()
    return res_pt, res_norm


#-----------------------------------------------------------------------------
def create_input_on_a_polygon6():
    angs = np.linspace(0.0, 2*np.pi, 6,endpoint = False)
    pts = []
    nrm = []
    radius = 2.
    for a in angs:
        curr_norm = np.array([np.cos(a), np.sin(a)])
        nrm.append(curr_norm)
        curr_pt = curr_norm.copy()
        pts.append(curr_pt*radius)
    c45 = np.cos(2.**0.5/2.0)
    nrm[1] = np.array([-c45, c45])
    #nrm[1] = np.array([0., 1.])
    #nrm[2] = np.array([c45, c45])
    nrm[4] = np.array([-1., 0.])
    return pts, nrm

#-----------------------------------------------------------------------------
def create_input_on_a_polygon5():
    c30 = 3.0**0.3/2.0
    pts = [np.array([-2., -2.] ),
           np.array([-1.,  3.] ),
           np.array([ 5.,  3.5]),
           np.array([ 2.,  0.] ),
           np.array([ 4., -3.] )
           ]
    nrm = [np.array([-1., -1.]),
           np.array([-1.,  1.]),
           np.array([ 1.,  1.]),

           np.array([ 1.,  0.]),
           #np.array([ -1,  0]),

           np.array([ 1., -1.])
           #np.array([ 0., -1.]),
           #np.array([  0,  1]),
           #np.array([  0,  -1]),
           #np.array([  c30, -0.5]), #rd
           #np.array([ c30, 0.5]), #ru
           #np.array([ -c30, 0.5]), #lu
           #np.array([ -c30, -0.5]), #ld
           ]

    nnrm = []
    for n in nrm:
        nnrm.append( n/np.linalg.norm(n))
    return pts, nnrm
#-----------------------------------------------------------------------------
def create_input_on_a_polygon5_rot():
    c30 = 3.0**0.3/2.0
    pts = [np.array([-2., -2.] ),
           np.array([-1.,  3.] ),
           np.array([ 5.,  3.5]),
           np.array([ 2.,  0.] ),
           np.array([ 4., -3.] )
           ]

    nrm = [np.array([1., -1.]),
           np.array([-1.,  -1.]),
           np.array([ -1.,  1.]),

           #np.array([ 1.,  0.]),
           np.array([ 0,  -1]),

           np.array([ 1., 1.])
           #np.array([ 0., -1.]),
           #np.array([  0,  1]),
           #np.array([  0,  -1]),
           #np.array([  c30, -0.5]), #rd
           #np.array([ c30, 0.5]), #ru
           #np.array([ -c30, 0.5]), #lu
           #np.array([ -c30, -0.5]), #ld
           ]


#-----------------------------------------------------------------------------
def create_input_on_a_square():
    pts = [np.array([ 0., 0.]),
           np.array([ 4., 0.]),
           np.array([ 4., 4.]),
           np.array([ 0., 4.])]
    nrm = [np.array([ -3., -1.]),
           np.array([  1., -3.]),
           np.array([  3.,  1.]),
           np.array([  -1.,  3.])]
    nnrm = []
    for n in nrm:
        nnrm.append( n/np.linalg.norm(n))
    return pts, nnrm

#-----------------------------------------------------------------------------
def create_input_on_a_square_clockwise():
    pts = [np.array([ 0., 0.]),
           np.array([ 0., 4.]),
           np.array([ 4., 4.]),
           np.array([ 4., 0.])]
    nrm = [np.array([ -1., 0.]),
           np.array([  0., 1.]),
           np.array([  1.,  0.]),
           np.array([  0.,  -1.])]
    nnrm = []
    for n in nrm:
        nnrm.append( n/np.linalg.norm(n))
    return pts, nnrm

#-----------------------------------------------------------------------------
def create_input_on_a_square_diff_norms():
    pts = [np.array([  0.,  0.]),
           np.array([ 10.,  0.]),
           np.array([ 10., 10.]),
           np.array([  0., 10.])]
    s30 = 0.5
    c30 = (3.0**0.5) / 2.0
    nrm = [np.array([ -1.0, -0.1]),
           np.array([ -1.0, -1.0]),
           np.array([  1.0,  0.9]),
           np.array([  0.1,  1.0])]
    #nrm = [np.array([ 0., -1.]),
    #       np.array([ 1.,  0.]),
    #       np.array([ 0.,  1.]),
    #       np.array([-1.,  0.])]
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
    n_of_iterations = 3
    bspline_average_export = bspline_average_export_v3
    b_open = False
    #subd_pts, subd_nrm = create_input_on_a_polygon5()
    #subd_pts, subd_nrm = create_input_on_a_polygon5_rot()
    #subd_pts, subd_nrm = create_input_on_a_polygon6()
    #subd_pts, subd_nrm = create_input_on_a_square()
    subd_pts, subd_nrm = create_input_on_a_square_clockwise()
    #subd_pts, subd_nrm = create_input_on_a_square_diff_norms()
    #subd_nrm = init_normals(subd_pts, b_open)
    orig_pts = subd_pts[:]

    bsubd_INS_pts,  bsubd_INS_nrm  = subd_pts[:], subd_nrm[:]
    bsubd_MLR2_pts, bsubd_MLR2_nrm = subd_pts[:], subd_nrm[:]
    bsubd_MLR3_pts, bsubd_MLR3_nrm = subd_pts[:], subd_nrm[:]
    bsubd_MLR5_pts, bsubd_MLR5_nrm = subd_pts[:], subd_nrm[:]
    bsubd_4pt_pts,  bsubd_4pt_nrm  = subd_pts[:], subd_nrm[:]

    csubd_INS_pts,  csubd_INS_nrm  = subd_pts[:], subd_nrm[:]
    csubd_MLR3_pts, csubd_MLR3_nrm = subd_pts[:], subd_nrm[:]
    csubd_MLR5_pts, csubd_MLR5_nrm = subd_pts[:], subd_nrm[:]
    csubd_4pt_pts,  csubd_4pt_nrm  = subd_pts[:], subd_nrm[:]

    lsubd_INS_pts,  lsubd_INS_nrm  = subd_pts[:], subd_nrm[:]
    lsubd_MLR3_pts, lsubd_MLR3_nrm = subd_pts[:], subd_nrm[:]
    lsubd_MLR5_pts, lsubd_MLR5_nrm = subd_pts[:], subd_nrm[:]
    lsubd_4pt_pts,  lsubd_4pt_nrm  = subd_pts[:], subd_nrm[:]

    opt_INS_pts,    opt_INS_nrm    = subd_pts[:], subd_nrm[:]
    opt_MLR3_pts,   opt_MLR3_nrm   = subd_pts[:], subd_nrm[:]
    opt_MLR5_pts,   opt_MLR5_nrm   = subd_pts[:], subd_nrm[:]
    opt_4pt_pts,    opt_4pt_nrm    = subd_pts[:], subd_nrm[:]

    corn_cut_pts,  corn_cut_nrm  = subd_pts[:], subd_nrm[:]

    for k in range(n_of_iterations):
        #--- Bezier Average
        #bsubd_INS_pts, bsubd_INS_nrm   = double_polygon(bsubd_INS_pts, bsubd_INS_nrm,
        #                                               True, b_open,
        #                                              bspline_average_export)
        #bsubd_MLR2_pts, bsubd_MLR2_nrm = subd_LR_one_step(bsubd_MLR2_pts, bsubd_MLR2_nrm, 
        #                                                  b_open, bspline_average_export, n_deg = 2)
        bsubd_MLR3_pts, bsubd_MLR3_nrm = subd_LR_one_step(bsubd_MLR3_pts, bsubd_MLR3_nrm, 
                                                          b_open, bspline_average_export)
        #bsubd_MLR5_pts, bsubd_MLR5_nrm = subd_LR_one_step(bsubd_MLR5_pts, bsubd_MLR5_nrm, 
        #                                                 b_open, bspline_average_export, n_deg = 5)
        #bsubd_4pt_pts, bsubd_4pt_nrm = subd_4PT_one_step(bsubd_4pt_pts, bsubd_4pt_nrm, 
        #                                                 b_open, bspline_average_export)
        
        #--- Circle Average
        #csubd_INS_pts, csubd_INS_nrm    = double_polygon(csubd_INS_pts, csubd_INS_nrm,
        #                                                 True, b_open,
        #                                                 circle_avg)
        csubd_MLR3_pts, csubd_MLR3_nrm  = subd_LR_one_step(csubd_MLR3_pts, csubd_MLR3_nrm, 
                                                           b_open, circle_avg)
        #csubd_MLR5_pts, csubd_MLR5_nrm  = subd_LR_one_step(csubd_MLR5_pts, csubd_MLR5_nrm, 
        #                                                   b_open, circle_avg, n_deg = 5)
        #csubd_4pt_pts, csubd_4pt_nrm = subd_4PT_one_step(csubd_4pt_pts, csubd_4pt_nrm, 
        #                                                 b_open, circle_avg)

        #--- Linear Average
        #lsubd_INS_pts, lsubd_INS_nrm    = double_polygon(lsubd_INS_pts, lsubd_INS_nrm,
        #                                                 True, b_open,
        #                                                 linear_avg)
        #lsubd_MLR3_pts, lsubd_MLR3_nrm  = subd_LR_one_step(lsubd_MLR3_pts, lsubd_MLR3_nrm, 
        #                                                   b_open, linear_avg)
        #lsubd_MLR5_pts, lsubd_MLR5_nrm  = subd_LR_one_step(lsubd_MLR5_pts, lsubd_MLR5_nrm, 
        #                                                   b_open, linear_avg, n_deg = 5)
        #lsubd_4pt_pts, lsubd_4pt_nrm = subd_4PT_one_step(lsubd_4pt_pts, lsubd_4pt_nrm, 
        #                                                 b_open, linear_avg)

        #--- Optimizer Average
        #opt_INS_pts, opt_INS_nrm    = double_polygon(opt_INS_pts, opt_INS_nrm,
        #                                              True, b_open,
        #                                              bopt_average_export)
        #opt_MLR3_pts, opt_MLR3_nrm  = subd_LR_one_step(opt_MLR3_pts, opt_MLR3_nrm, 
        #                                               b_open, bopt_average_export)
        #opt_MLR5_pts, opt_MLR5_nrm   = subd_LR_one_step(opt_MLR5_pts, opt_MLR5_nrm, 
        #                                                b_open, bopt_average_export, n_deg = 5)
        #opt_4pt_pts, opt_4pt_nrm    = subd_4PT_one_step(opt_4pt_pts, opt_4pt_nrm, 
        #                                                b_open, bopt_average_export)

        #--- Corner cutting
        #corn_cut_pts, corn_cut_nrm  = subd_CornerCutting_one_step( corn_cut_pts, corn_cut_nrm, 
        #                                                           b_open, bspline_average_export, 1./3.)

    fig = plt.figure()#figsize=(8,8), dpi=100, frameon = False)
    #frame1 = plt.gca()
    #frame1.axes.get_xaxis().set_visible(False)
    #frame1.axes.get_yaxis().set_visible(False)

    #plot_pts_and_norms(bsubd_INS_pts, bsubd_INS_nrm, b_open, False, clr='c', linewidth=1.0, linestyle='solid')
    #plot_pts_and_norms(bsubd_MLR2_pts, bsubd_MLR2_nrm, b_open, True, clr='#9d68e6', linewidth=1.0, linestyle='solid')
    plot_pts_and_norms(bsubd_MLR3_pts, bsubd_MLR3_nrm, b_open, True, clr='#ff6e61', linewidth=1.0, linestyle='solid')
    #plot_pts_and_norms(bsubd_MLR5_pts, bsubd_MLR5_nrm, b_open, False, clr='#9d9bd9', linewidth=1.0, linestyle='solid')
    #plot_pts_and_norms(bsubd_4pt_pts, bsubd_4pt_nrm, b_open, True, clr='#4441a9', linewidth=1.0, linestyle='solid')

    #plot_pts_and_norms(csubd_INS_pts, csubd_INS_nrm, b_open, True, clr='#43a941', linewidth=1.0, linestyle='solid')
    plot_pts_and_norms(csubd_MLR3_pts, csubd_MLR3_nrm, b_open, True, clr='b', linewidth=1.0, linestyle='solid') #'#90d876'
    #plot_pts_and_norms(csubd_MLR5_pts, csubd_MLR5_nrm, b_open, True, clr='#43a941', linewidth=1.0, linestyle='solid')
    #plot_pts_and_norms(csubd_4pt_pts, csubd_4pt_nrm, b_open, True, clr='b', linewidth=1.0, linestyle='solid')

    #plot_pts_and_norms(lsubd_INS_pts,  lsubd_INS_nrm, b_open, False, clr='#ff93de', linewidth=1.0, linestyle='solid')
    #plot_pts_and_norms(lsubd_MLR3_pts, lsubd_MLR3_nrm, b_open, False, clr='g', linewidth=1.0, linestyle='solid')
    #plot_pts_and_norms(lsubd_MLR5_pts, lsubd_MLR5_nrm, b_open, False, clr='#d286bb', linewidth=1.0, linestyle='solid')
    #plot_pts_and_norms(lsubd_4pt_pts,  lsubd_4pt_nrm, b_open, False, clr='#ff93de', linewidth=1.0, linestyle='solid')

    #plot_pts_and_norms(opt_ins_pts, opt_ins_nrm, b_open, True, clr='y', linewidth=1.0, linestyle='solid')
    #plot_pts_and_norms(opt_MLR3_pts, opt_MLR3_nrm, b_open, True, clr='g', linewidth=1.0, linestyle='solid')
    #plot_pts_and_norms(opt_MLR5_pts, opt_MLR5_nrm, b_open, True, clr='#cfcf30', linewidth=1.0, linestyle='solid')
    #plot_pts_and_norms(opt_4pt_pts, opt_4pt_nrm, b_open, True, clr='y', linewidth=1.0, linestyle='solid')

    #plot_pts_and_norms(corn_cut_pts, corn_cut_nrm, b_open, True, clr='#45ff02', linewidth=1.0, linestyle='solid')

    plot_pts_and_norms(orig_pts, subd_nrm, b_open, True, clr='k', linewidth=1.0, linestyle='dotted')

    plt.axis('equal')
    plt.xlim([-4, 6])
    plt.ylim([-5, 6])
    plt.axis('off')
    plt.show()


    #xs1 = compute_x_axes_values(bsubd_MLR2_pts, b_open, n_of_iterations)
    #slopes_ins = compute_slopes(bsubd_INS_pts, b_open, n_of_iterations)
    #curv_ins = compute_curvature(bsubd_INS_pts, b_open, n_of_iterations)
    #slopes_mlr2 = compute_slopes(bsubd_MLR2_pts, b_open, n_of_iterations)
    #curv_mlr2 = compute_curvature(bsubd_MLR2_pts, b_open, n_of_iterations)
    #slopes_mlr3 = compute_slopes(bsubd_MLR3_pts, b_open, n_of_iterations)
    #curv_mlr3 = compute_curvature(bsubd_MLR3_pts, b_open, n_of_iterations)
    #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

    #ax1.plot(xs1, slopes_ins, 'c')
    #ax1.set_xlabel('Slope')
    #ax1.xaxis.set_label_position('top') 
    #ax1.set_ylabel('MLR1')

    #ax2.plot(xs1, curv_ins, 'r')
    #ax2.xaxis.set_label_position('top') 
    #ax2.set_xlabel('Curvature')

    #ax3.plot(xs1, slopes_mlr2, 'c', label = 'MLR2')
    #ax3.set_ylabel('MLR2')

    #ax4.plot(xs1, curv_mlr2, 'r')
    
    ##ax5.plot(xs1, slopes_mlr3, 'c', label = 'MLR3')
    ##ax5.set_ylabel('MLR3')
    ##ax6.plot(xs1, curv_mlr3, 'r')
    #plt.tight_layout()
    #plt.show()

#-----------------------------------------------------------------------------
def bezier_test():
    a = np.array([0.0,  0.0])
    b = np.array([1.0,  1.0])
    c = np.array([2.0, -1.0])
    d = np.array([3.0,  0.0])
    print 'Orig half-length = ', estimate_bezier_length_v2(a,b,c,d)/2.
    lb = split_cubic_bezier_get_left(a,b,c,d, 0.5)
    print 'Left length      = ', estimate_bezier_length_v2(lb[0],lb[1],lb[2],lb[3])

    res000_pt, res000_norm = eval_cubic_bezier(a, b, c, d, 0.0 )
    res025_pt, res025_norm = eval_cubic_bezier(a, b, c, d, 0.25)
    res050_pt, res050_norm = eval_cubic_bezier(a, b, c, d, 0.5 )
    res075_pt, res075_norm = eval_cubic_bezier(a, b, c, d, 0.75)
    res100_pt, res100_norm = eval_cubic_bezier(a, b, c, d, 1.0 )
    res000_norm = np.array([res000_norm[1], -res000_norm[0]])
    res025_norm = np.array([res025_norm[1], -res025_norm[0]])
    res050_norm = np.array([res050_norm[1], -res050_norm[0]])
    res075_norm = np.array([res075_norm[1], -res075_norm[0]])
    res100_norm = np.array([res100_norm[1], -res100_norm[0]])
    print get_angle(res000_norm[0], res000_norm[1]) * 180. / np.pi
    print get_angle(res025_norm[0], res025_norm[1])* 180. / np.pi
    print get_angle(res050_norm[0], res050_norm[1])* 180. / np.pi
    print get_angle(res075_norm[0], res075_norm[1])* 180. / np.pi
    print get_angle(res100_norm[0], res100_norm[1])* 180. / np.pi
    #DEBUG = True
    DEBUG = 'IN_DEBUG' in globals()
    if DEBUG and IN_DEBUG:
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
        nr000 = plt.Arrow(res000_pt[0], res000_pt[1], 
                          res000_norm[0]*.1, res000_norm[1]*.1, width=0.03, fc='g', ec='none' )
        nr025 = plt.Arrow(res025_pt[0], res025_pt[1], 
                          res025_norm[0]*.1, res025_norm[1]*.1, width=0.03, fc='g', ec='none' )
        nr050 = plt.Arrow(res050_pt[0], res050_pt[1], 
                          res050_norm[0]*.1, res050_norm[1]*.1, width=0.03, fc='g', ec='none' )
        nr075 = plt.Arrow(res075_pt[0], res075_pt[1], 
                          res075_norm[0]*.1, res075_norm[1]*.1, width=0.03, fc='g', ec='none' )
        nr100 = plt.Arrow(res100_pt[0], res100_pt[1], 
                          res100_norm[0]*.1, res100_norm[1]*.1, width=0.03, fc='g', ec='none' )
        plt.gca().add_patch(nr000)
        plt.gca().add_patch(nr025)
        plt.gca().add_patch(nr050)
        plt.gca().add_patch(nr075)
        plt.gca().add_patch(nr100)
        patch = patches.PathPatch(crv, facecolor='none', lw=2)
        plt.gca().add_patch(patch)
        plt.axis('equal')
        plt.axis([-1.0, 2., -1.0, 2])
        plt.show()

#-----------------------------------------------------------------------------
def one_pair_test():
    points  = [(  0.0, 0.0),
               (  0.0, 1.0)]
    normals = [(  1.0, 1.0),
               (  1.0, 1.0)]
    normals[0] /= np.linalg.norm(normals[0])
    normals[1] /= np.linalg.norm(normals[1])
    p0 = np.array([points[0][0], points[0][1]])
    p1 = np.array([points[1][0], points[1][1]])
    n0 = np.array([normals[0][0], normals[0][1]])
    n1 = np.array([normals[1][0], normals[1][1]])
    bspline_average_v3(0.5, p0, p1, n0, n1)
    #bspline_average_v1(1.125, p0, p1, n0, n1)
    #bspline_average_v2(-0.125, p0, p1, n0, n1)

#-----------------------------------------------------------------------------
def one_line_test():
    n_of_pnps  = 5
    orig_pts   = [np.array([float(i), 0.]) for i in range(n_of_pnps)]
    orig_pts.append(np.array([float(n_of_pnps-1)/2., -3.]))
    orig_norms = [np.array([0., 1. if i%2==0 else -1.]) for i in range(n_of_pnps)]
    orig_norms.append(np.array([0., -1.]))
    n_of_iterations = 5
    b_open = False
    bsubd_INS_pts, bsubd_INS_nrm = orig_pts[:], orig_norms[:]
    bsubd_MLR2_pts, bsubd_MLR2_nrm = orig_pts[:], orig_norms[:]
    for k in range(n_of_iterations):
        bsubd_INS_pts, bsubd_INS_nrm   = double_polygon(bsubd_INS_pts, 
                                                        bsubd_INS_nrm,
                                                        True, 
                                                        b_open,
                                                  bspline_average_export_v3)
        bsubd_MLR2_pts, bsubd_MLR2_nrm = subd_LR_one_step(bsubd_MLR2_pts, 
                                                          bsubd_MLR2_nrm, 
                                                          b_open, 
                                                  bspline_average_export_v3, 
                                                     n_deg = 2)
    plot_pts_and_norms(bsubd_INS_pts, bsubd_INS_nrm, b_open, 
                       True, clr='c', linewidth=1.0, linestyle='solid')
    plot_pts_and_norms(bsubd_MLR2_pts, bsubd_MLR2_nrm, b_open, 
                       True, clr='b', linewidth=1.0, linestyle='solid')
    plot_pts_and_norms(orig_pts, orig_norms, b_open, 
                       True, clr='k', linewidth=1.0, linestyle='dotted')
    plt.axis('equal')
    plt.xlim([-5, 6])
    plt.ylim([-5, 6])
    #plt.axis('off')
    plt.show()
#-----------------------------------------------------------------------------
def max_dist_test():
    global IN_DEBUG
    n0_angs = np.arange(0.0, 2*np.pi, np.pi/180., float)
    n1_angs = np.arange(0.0, 2*np.pi, np.pi/180., float)
    p0 = np.array([0., 0.])
    p1 = np.array([1., 0.])
    p0_r_max = (0.0, -1, -1)
    p1_r_max = (0.0, -1, -1)
    b_inter = False
    j = b_two_segments_intersect(p0, np.array([0., 1.]), p1, np.array([-1.1, 0.1]))
    for n0_ang_idx in range(len(n0_angs)):
        for n1_ang_idx in range(len(n1_angs)):
            n0_ang = n0_angs[n0_ang_idx]
            n1_ang = n1_angs[n1_ang_idx]
            n0 = np.array([np.cos(n0_ang), np.sin(n0_ang)])
            n1 = np.array([np.cos(n1_ang), np.sin(n1_ang)])
            
            ca_pt, ca_norm, ca_center, ca_radius, ca_beta0, ca_beta1 = \
                        circle_avg(0.5, 0.5, True, p0, p1, n0, n1)
            theta = get_angle_between(n0, n1)
            der_length = 4./3. * np.tan(theta/4.) * ca_radius
            t0 = np.array([n0[1], -n0[0]]) * der_length
            t1 = np.array([-n1[1], n1[0]]) * der_length
            b_inter = b_two_segments_intersect(p0, t0, p1, t1)
            if b_inter:
                print 'Intersection. Der len = ', der_length, 't0 = (', t0[0], t0[1], ') t1 = (', t1[0], t1[1], ')'
    if not b_inter:
        print 'No intersectons found'
    #        r_pt, r_norm = bspline_average_v3(0.5, p0, p1, n0, n1)
    #        dist_p0_r = get_dist(p0, r_pt)
    #        if dist_p0_r > p0_r_max[0]:
    #            p0_r_max = (dist_p0_r, n0_ang_idx, n1_ang_idx)
    #            #print 'Updating P0 max, with ', dist_p0_r, 'at', n0_ang_idx, ',', n1_ang_idx
    #        dist_p1_r = get_dist(p1, r_pt)
    #        if dist_p1_r > p1_r_max[0]:
    #            p1_r_max = (dist_p1_r, n0_ang_idx, n1_ang_idx)
    #            #print 'Updating P1 max, with ', dist_p1_r, 'at', n0_ang_idx, ',', n1_ang_idx
    #print 'P0, max dist ', p0_r_max[0], ' n0 ', n0_angs[p0_r_max[1]], ' n1 ', n1_angs[p0_r_max[2]] 
    #print 'P1, max dist ', p1_r_max[0], ' n0 ', n0_angs[p1_r_max[1]], ' n1 ', n1_angs[p1_r_max[2]] 
    #IN_DEBUG = True
    #n0_ang = n0_angs[p0_r_max[1]]
    #n1_ang = n1_angs[p0_r_max[2]]
    #n0 = np.array([np.cos(n0_ang), np.sin(n0_ang)])
    #n1 = np.array([np.cos(n1_ang), np.sin(n1_ang)])
    #bspline_average_v3(0.5, p0, p1, n0, n1)
#-----------------------------------------------------------------------------
def derivative_length_graph():
    x_s = np.linspace(0.0, np.pi-0.1, 100)
    y_s = [4./3. * np.tan(x/4.) * 1. / (2*np.cos(x/2.)) for x in x_s]
    y1_s = [np.tan(x/4.) for x in x_s]
    fig = plt.figure()
    plt.plot(x_s, y_s, 'c')
    plt.plot(x_s, y1_s, 'g')
    plt.axis('equal')
    plt.show()#block=False)

#-----------------------------------------------------------------------------
def naive_normals_case():
   pts = [np.array([ -10., 0.] ),
           np.array([ 10., 0.] ),
           np.array([ 10,  1.] )]
   b_open = True
   nrms = init_normals( pts, b_open )
   a = 5
#-----------------------------------------------------------------------------
def solve_equation():
    xs = np.arange(M.atan(1./3.), M.atan(2./3.), 0.0000001)
    def equation(x):
        return np.tan(x) - 1./(3. * (np.cos( ( np.pi/2. + x ) / 4. ) )**2. )
    ys = [equation(x) for x in xs]
    def pair_comp(p):
        return np.abs(p[1])
    min_pair = min(zip(xs,ys), key=pair_comp)
    print 'Minimal x =', min_pair[0],'y =', min_pair[1]
    plt.plot(xs, ys)
    plt.show()
#-----------------------------------------------------------------------------
if __name__ == "__main__":
    global IN_DEBUG
    IN_DEBUG = True
    IN_DEBUG = False
    solve_equation()
    #naive_normals_case()
    #max_dist_test()
    #bezier_test()
    #build_curves()
    #one_pair_test()
    #one_line_test()
    #derivative_length_graph()
#============================ END OF FILE ====================================
