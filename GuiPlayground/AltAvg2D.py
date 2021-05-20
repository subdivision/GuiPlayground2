import math as M
import numpy as np
import matplotlib.pyplot as plt
from CSubd2D import get_weighted_angle, get_angle, get_angle_between, vec_eeq, eeq

#------------------------------------------------------------------------------
def line_intersect(p0, p1, m0=None, m1=None, q0=None, q1=None):
    ''' intersect 2 lines given 2 points and 
        (either associated slopes or one extra point)
        Inputs:
        p0 - first point of first line [x,y]
        p1 - fist point of second line [x,y]
        m0 - slope of first line
        m1 - slope of second line
        q0 - second point of first line [x,y]
        q1 - second point of second line [x,y]
    '''
    if m0 is  None:
        if q0 is None:
            raise ValueError('either m0 or q0 is needed')
        dy = q0[1] - p0[1]
        dx = q0[0] - p0[0]
        lhs0 = [-dy, dx]
        rhs0 = p0[1] * dx - dy * p0[0]
    else:
        lhs0 = [-m0, 1]
        rhs0 = p0[1] - m0 * p0[0]

    if m1 is  None:
        if q1 is None:
            raise ValueError('either m1 or q1 is needed')
        dy = q1[1] - p1[1]
        dx = q1[0] - p1[0]
        lhs1 = [-dy, dx]
        rhs1 = p1[1] * dx - dy * p1[0]
    else:
        lhs1 = [-m1, 1]
        rhs1 = p1[1] - m1 * p1[0]

    a = np.array([lhs0, 
                  lhs1])

    b = np.array([rhs0, 
                  rhs1])
    try:
        px = np.linalg.solve(a, b)
    except:
        px = np.array([np.nan, np.nan])

    return px

#------------------------------------------------------------------------------
def create_circle(p0, p1, n0):
    slope0 = n0[1]/n0[0]
    slope1 = (p1[1] - p0[1]) / (p1[0] - p0[0])
    slope1 = -1.0 / slope1
    mid_pt = (p0 + p1)/2.0
    cntr = line_intersect(p0, mid_pt, m0 = slope0, m1 = slope1)
    radius = np.linalg.norm(p0 - cntr)
    return cntr, radius

#------------------------------------------------------------------------------
def create_circle_v2(p0, p1, n0):
    p1p0_vec = p1-p0
    p1p0_len = np.linalg.norm(p1p0_vec)
    p1p0_vec /= p1p0_len
    alpha = get_angle_between(n0, p1p0_vec)
    alpha = M.pi - alpha
    radius = (p1p0_len/2.)/M.cos(alpha)
    cntr = p0 - n0*radius
    return cntr, radius

#------------------------------------------------------------------------------
def create_circle_v3(p, n, q):
    natural_norm = p-q
    radius = np.linalg.norm(natural_norm)
    center = p + natural_norm
    natural_norm /= radius
    center = q if vec_eeq(natural_norm,n) else center
    return center, radius 

#------------------------------------------------------------------------------
def find_circle_center(p0, p1, radius, old_center):
    mid_pt = (p0 + p1)/2.0
    mid_perp = np.array([-(p0 - p1)[1], (p0 - p1)[0]])
    segm_length = np.linalg.norm(p0-p1)
    mid_perp /= np.linalg.norm(mid_perp)
    length_to_center = (radius ** 2 - (segm_length/2.) ** 2) ** 0.5
    cntr1 = mid_pt + mid_perp * length_to_center
    cntr2 = mid_pt - mid_perp * length_to_center
    test1a = np.linalg.norm(p0 - cntr1)
    test1b = np.linalg.norm(p1 - cntr1)
    test2a = np.linalg.norm(p0 - cntr2)
    test2b = np.linalg.norm(p1 - cntr2)
    d1 = np.linalg.norm(old_center - cntr1)
    d2 = np.linalg.norm(old_center - cntr2)
    return (cntr1 if d1 < d2 else cntr2)

#------------------------------------------------------------------------------
def alt_average(t, p0, p1, n0, n1):
    c0,r0 = create_circle(p0, p1, n0)
    c1,r1 = create_circle(p1, p0, n1)
    
    #new_cntr = c0*(1. - t) + c1*t    
    #-- v1
    # new_radius = r0*(1. - t) + r1*t
    #-- v2
    #d0 = np.linalg.norm(p0 - new_cntr)
    #d1 = np.linalg.norm(p1 - new_cntr)
    #new_radius = d0*(1. - t) + d1*t
    #-- v3
    res_radius = r0*(1. - t) + r1*t
    res_cntr = find_circle_center(p0, p1, res_radius, c0)
    vec_p0 = p0 - res_cntr
    vec_p0 /= np.linalg.norm(vec_p0)
    vec_p1 = p1 - res_cntr
    vec_p1 /= np.linalg.norm(vec_p1)
    res_ang = get_weighted_angle(1. - t, t, vec_p0, vec_p1)
    #res_ang = get_weighted_angle(1. - t, t, n0, n1)
    res_norm = np.array([M.cos(res_ang), M.sin(res_ang)])
    res_pt = res_cntr + res_norm * res_radius

    #draw_debug = True
    draw_debug = False

    if draw_debug:
        ar0 = plt.Arrow( p0[0], p0[1], n0[0], n0[1], width=0.03, fc='b', ec='none' )
        ar1 = plt.Arrow( p1[0], p1[1], n1[0], n1[1], width=0.03, fc='g', ec='none' )
        pt0 = plt.Circle( (p0[0], p0[1]), radius=0.02, fc='b', ec='none')
        pt1 = plt.Circle( (p1[0], p1[1]), radius=0.02, fc='g', ec='none')
        cr0 = plt.Circle( (c0[0], c0[1]), radius=r0, fc='none', ec='b')
        cr1 = plt.Circle( (c1[0], c1[1]), radius=r1, fc='none', ec='g')
        crr = plt.Circle( (res_cntr[0], res_cntr[1]), radius=res_radius, fc='none', ec='r')
        arr = plt.Arrow( res_pt[0], res_pt[1], res_norm[0], res_norm[1], width=0.03, fc='r', ec='none' )
        ptr = plt.Circle( (res_pt[0], res_pt[1]), radius=0.02, fc='r', ec='none')
        plt.gca().add_patch(ar0)
        plt.gca().add_patch(ar1)
        plt.gca().add_patch(pt0)
        plt.gca().add_patch(pt1)
        plt.gca().add_patch(cr0)
        plt.gca().add_patch(cr1)
        plt.gca().add_patch(crr)
        plt.gca().add_patch(arr)
        plt.gca().add_patch(ptr)
        plt.axis('scaled')
        plt.show()

    return res_pt, res_norm

#-----------------------------------------------------------------------------
def alt_average_v2(t, p0, p1, n0, n1):

    np0, np1, nn0, nn1, gamma, ofx, ofy = put_to_zero(p0, p1, n0, n1)
    
    c0,r0 = create_circle_v2(np0, np1, nn0)
    c1,r1 = create_circle_v2(np1, np0, nn1)
    
    res_radius = r0*(1. - t) + r1*t
    res_cntr = find_circle_center(np0, np1, res_radius, c0)
    vec_p0 = np0 - res_cntr
    vec_p0 /= np.linalg.norm(vec_p0)
    vec_p1 = np1 - res_cntr
    vec_p1 /= np.linalg.norm(vec_p1)
    res_ang2 = get_weighted_angle(1. - t, t, vec_p0, vec_p1)
    res_ang = get_weighted_angle(1. - t, t, nn0, nn1)
    theta = get_angle_between( nn0, nn1 )
    res_ang1 = theta*(1.-t)
    res_norm = np.array([M.cos(res_ang), M.sin(res_ang)])
    res_pt = res_cntr + res_norm * res_radius

    draw_debug = True
    #draw_debug = False

    if draw_debug:
        ar0 = plt.Arrow( np0[0], np0[1], nn0[0], nn0[1], width=0.03, fc='b', ec='none' )
        ar1 = plt.Arrow( np1[0], np1[1], nn1[0], nn1[1], width=0.03, fc='g', ec='none' )
        pt0 = plt.Circle( (np0[0], np0[1]), radius=0.02, fc='b', ec='none')
        pt1 = plt.Circle( (np1[0], np1[1]), radius=0.02, fc='g', ec='none')
        cr0 = plt.Circle( (c0[0], c0[1]), radius=r0, fc='none', ec='b')
        cr1 = plt.Circle( (c1[0], c1[1]), radius=r1, fc='none', ec='g')
        crr = plt.Circle( (res_cntr[0], res_cntr[1]), radius=res_radius, fc='none', ec='r')
        arr = plt.Arrow( res_pt[0], res_pt[1], res_norm[0], res_norm[1], width=0.03, fc='r', ec='none' )
        ptr = plt.Circle( (res_pt[0], res_pt[1]), radius=0.02, fc='r', ec='none')
        plt.gca().add_patch(ar0)
        plt.gca().add_patch(ar1)
        plt.gca().add_patch(pt0)
        plt.gca().add_patch(pt1)
        plt.gca().add_patch(cr0)
        plt.gca().add_patch(cr1)
        plt.gca().add_patch(crr)
        plt.gca().add_patch(arr)
        plt.gca().add_patch(ptr)
        plt.axis('scaled')
        plt.show()

    res_pt, res_norm = put_back(res_pt, res_norm, gamma, ofx, ofy )
    p0, n0 = put_back(p0, n0, gamma, ofx, ofy )
    p1, n1 = put_back(p1, n1, gamma, ofx, ofy )

    return res_pt, res_norm

#-----------------------------------------------------------------------------
def put_to_zero(p0, p1, n0, n1):
    offset_x = p0[0]
    offset_y = p0[1]
    gamma = p1-p0
    p1p0_len = np.linalg.norm(gamma)
    gamma /= p1p0_len

    cw_rot_mtrx = np.array([[ gamma[0], gamma[1]],
                            [-gamma[1], gamma[0]]])
    nn0 = cw_rot_mtrx.dot(n0)
    nn1 = cw_rot_mtrx.dot(n1)
    np0 = np.array([0., 0.])
    np1 = np.array([p1p0_len, 0.])
    return np0, np1, nn0, nn1, gamma, offset_x, offset_y

#-----------------------------------------------------------------------------
def put_back(npt, nnorm, gamma, offset_x, offset_y):
    cw_rot_mtrx = np.array([[ gamma[0], -gamma[1]],
                            [ gamma[1],  gamma[0]]])
    norm = cw_rot_mtrx.dot(nnorm)
    pt = cw_rot_mtrx.dot(npt) + np.array([offset_x, offset_y])
    return pt, norm

#-----------------------------------------------------------------------------
def sliding_circ_average_export( t1, t2, b_open, p1, p2, n1, n2 ):
    p, n = sliding_circ_average(t1, p1,p2,n1,n2)
    return p,n, p, 0,0,0

#-----------------------------------------------------------------------------
def sliding_circ_average(t, p0, p1, n0, n1):
    m0 = n0[1]/n0[0] if not eeq(n0[0], 0) else None
    m1 = n1[1]/n1[0] if not eeq(n1[0], 0) else None
    if m0 == m1:
        res_ang = get_weighted_angle(1. - t, t, n0, n1)
        res_norm = np.array([M.cos(res_ang), M.sin(res_ang)])
        res_pt = p0*(1. - t) + p1*t
    else:  
        q0 = p0 + n0 if m0 == None else None
        q1 = p1 + n1 if m1 == None else None
        q = line_intersect(p0, p1, m0 = m0, m1 = m1, q0 = q0, q1 = q1)
        c0, r0 = create_circle_v3(p0, n0, q)
        c1, r1 = create_circle_v3(p1, n1, q)
        res_radius = r0*(1. - t) + r1*t
        res_cntr = c0*(1. - t) + c1*t
        res_ang = get_weighted_angle(1. - t, t, n0, n1)
        res_norm = np.array([M.cos(res_ang), M.sin(res_ang)])
        res_pt = res_cntr + res_norm * res_radius
    
    #draw_debug = True
    draw_debug = False

    if draw_debug:
        ar0 = plt.Arrow( p0[0], p0[1], n0[0], n0[1], width=0.03, fc='b', ec='none' )
        ar1 = plt.Arrow( p1[0], p1[1], n1[0], n1[1], width=0.03, fc='g', ec='none' )
        pt0 = plt.Circle( (p0[0], p0[1]), radius=0.02, fc='b', ec='none')
        pt1 = plt.Circle( (p1[0], p1[1]), radius=0.02, fc='g', ec='none')
        cr0 = plt.Circle( (c0[0], c0[1]), radius=r0, fc='none', ec='b')
        cr1 = plt.Circle( (c1[0], c1[1]), radius=r1, fc='none', ec='g')
        crr = plt.Circle( (res_cntr[0], res_cntr[1]), radius=res_radius, fc='none', ec='r')
        arr = plt.Arrow( res_pt[0], res_pt[1], res_norm[0], res_norm[1], width=0.03, fc='r', ec='none' )
        ptr = plt.Circle( (res_pt[0], res_pt[1]), radius=0.02, fc='r', ec='none')
        plt.gca().add_patch(ar0)
        plt.gca().add_patch(ar1)
        plt.gca().add_patch(pt0)
        plt.gca().add_patch(pt1)
        plt.gca().add_patch(cr0)
        plt.gca().add_patch(cr1)
        plt.gca().add_patch(crr)
        plt.gca().add_patch(arr)
        plt.gca().add_patch(ptr)
        plt.axis('scaled')
        plt.show()

    return res_pt, res_norm

#-----------------------------------------------------------------------------
if __name__ == "__main__":
    p0 = np.array([0.,0.])
    p1 = np.array([1.,0.])
    n0 = np.array([-(2.**0.5)/2.,(2.**0.5)/2.])
    n1 = np.array([(3.**0.5)/2., 0.5])

    pt05, norm05 = sliding_circ_average(0.5, p0, p1, n0, n1)
    pt025, norm025 = sliding_circ_average(0.25, p0, p1, n0, n1)

    pt05_05, norm05_05 = sliding_circ_average(0.5, p0, pt05, n0, norm05)


    c0,r0 = create_circle(p0, p1, n0)
    c1,r1 = create_circle(p1, p0, n1)

    pts = []
    norms = []
    n_steps = 20
    one_step = 1. / n_steps
    for i in range(n_steps+1):
        t = one_step*i
        curr_pt, curr_norm = alt_average(t, p0, p1, n0, n1)
        pts.append(curr_pt)
        norms.append(curr_norm)

    ar0 = plt.Arrow( p0[0], p0[1], n0[0], n0[1], width=0.03, fc='b', ec='none' )
    ar1 = plt.Arrow( p1[0], p1[1], n1[0], n1[1], width=0.03, fc='g', ec='none' )
    pt0 = plt.Circle( (p0[0], p0[1]), radius=0.02, fc='b', ec='none')
    pt1 = plt.Circle( (p1[0], p1[1]), radius=0.02, fc='g', ec='none')
    cr0 = plt.Circle( (c0[0], c0[1]), radius=r0, fc='none', ec='b')
    cr1 = plt.Circle( (c1[0], c1[1]), radius=r1, fc='none', ec='g')

    plt.gca().add_patch(ar0)
    plt.gca().add_patch(ar1)
    plt.gca().add_patch(pt0)
    plt.gca().add_patch(pt1)
    plt.gca().add_patch(cr0)
    plt.gca().add_patch(cr1)
    for i in range(len(pts)):
        curr_pt = pts[i]
        curr_norm = norms[i]
        #gr_pt = plt.plot([curr_pt[0],next_pt[0]], [curr_pt[1], next_pt[1]], 'r' )
        gr_nr = plt.Arrow(curr_pt[0], curr_pt[1], curr_norm[0], curr_norm[1], 
                  width=0.03, fc='r', ec='none' )
        #plt.gca().add_patch(gr_pt)
        plt.gca().add_patch(gr_nr)
    
    plt.axis('scaled')
    plt.show()
