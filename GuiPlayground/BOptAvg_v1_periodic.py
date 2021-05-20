#import numpy as np
#from scipy import interpolate

#import matplotlib.pyplot as plt

##x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
##y = np.sin(x)

##ctr =np.array( [(3 , 1), (2.5, 4), (0, 1), (-2.5, 4),
##                (-3, 0), (-2.5, -4), (0, -1), (2.5, -4), (3, -1)])
#ctr =np.array( [(-1, 0), (0, 1), (1, 0), (0, -1), (-1, 0)] )
#x=ctr[:,0]
#y=ctr[:,1]

##x=np.append(x,x[0])
##y=np.append(y,y[0])

#tck,u = interpolate.splprep([x,y],k=3, s = 0, per=1)
#u=np.linspace(0,1,num=50,endpoint=True)
#out = interpolate.splev(u,tck)

#plt.figure()
#plt.plot(x, y, 'ro', out[0], out[1], 'b')
##plt.legend(['Points', 'Interpolated B-spline', 'True'],loc='best')
#plt.axis([min(x)-1, max(x)+1, min(y)-1, max(y)+1])
#plt.title('B-Spline interpolation')
#cr1 = plt.Circle( (0., 0.), radius=1., fc='none', ec='y')
#plt.gca().add_patch(cr1)

#plt.axis('equal')
#plt.show()



#=============================================================================

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from CSubd2D import *

points  = [(  0.0, 0.0),
           (  1.0, 0.0)]
normals = [(  -1.0,  1.0),
           (   1.0,  1.0)]
normals[0] /= np.linalg.norm(normals[0])
normals[1] /= np.linalg.norm(normals[1])
p0 = np.array([points[0][0], points[0][1]])
p1 = np.array([points[1][0], points[1][1]])
n0 = np.array([normals[0][0], normals[0][1]])
n1 = np.array([normals[1][0], normals[1][1]])

ca_pt, ca_norm, ca_center, ca_radius, ca_beta0, ca_beta1 = \
  circle_avg(0.5, 0.5, True, p0, p1, n0, n1)

sq2 = (2. ** 0.5)/2.
#ctrl_pts = [ca_center + np.array([ -1.,  0.])*ca_radius, \
#            ca_center + np.array([-sq2,  sq2])*ca_radius, \
#            ca_center + np.array([  0.,  1.])*ca_radius, \
#            ca_center + np.array([ sq2,  sq2])*ca_radius, \
#            ca_center + np.array([  1.,  0.])*ca_radius, \
#            ca_center + np.array([ sq2, -sq2])*ca_radius, \
#            ca_center + np.array([  0.,  -1.])*ca_radius, \
#            ca_center + np.array([-sq2, -sq2])*ca_radius ]

ca_beta0, ca_beta1 = ca_beta1, ca_beta0
g0 = ca_beta0
part1_inc = (ca_beta1 - ca_beta0)/ 3.
g1 = ca_beta0 + part1_inc
g2 = ca_beta0 + 2. * part1_inc
g3 = ca_beta1
part2_inc = (360 - (ca_beta1 - ca_beta0))/ 5.
g4 = ca_beta1 + 1. * part2_inc
g5 = ca_beta1 + 2. * part2_inc
g6 = ca_beta1 + 3. * part2_inc
g7 = ca_beta1 + 4. * part2_inc

gs = [g0, g1, g2, g3, g4, g5, g6, g7]
ctrl_pts = []
for g in gs:
    g = g * np.pi / 180.0
    ctrl_pts.append(ca_center + np.array([ np.cos(g), np.sin(g)])*ca_radius)

degree = 3
ctrl_pts = ctrl_pts + ctrl_pts[0:degree+1]
ctrl_pts = np.array(ctrl_pts)
n_ctrl_pts = len(ctrl_pts)
ctrl_pts_x = ctrl_pts[:,0]
ctrl_pts_y = ctrl_pts[:,1]

t = range(len(ctrl_pts_x))
t = [(g-g0)/360. for g in gs]
for i in range(0, degree+1):
    t.append(1.0 + (gs[i] - gs[0])/360.)
bspl_x = si.make_interp_spline(t,ctrl_pts_x,degree)
bspl_y = si.make_interp_spline(t,ctrl_pts_y,degree)

ctrl_pts_x, ctrl_pts_y = bspl_x.tck[1], bspl_y.tck[1]
#ipl_t = np.linspace(1.0, len(ctrl_pts) - degree, 1000)
ipl_t = np.linspace(t[0], t[-1], 1000)

x_i, y_i = bspl_x(ipl_t), bspl_y(ipl_t)


#t = range(len(ctrl_pts_x))
#ipl_t = np.linspace(1.0, len(ctrl_pts) - degree, 1000)

#x_tup = si.splrep(t, ctrl_pts_x, k=degree, per=1)
#y_tup = si.splrep(t, ctrl_pts_y, k=degree, per=1)
#x_list = list(x_tup)
#xl = ctrl_pts_x.tolist()
#x_list[1] = xl + [0.0, 0.0, 0.0, 0.0]

#y_list = list(y_tup)
#yl = ctrl_pts_y.tolist()
#y_list[1] = yl + [0.0, 0.0, 0.0, 0.0]

#x_i = si.splev(ipl_t, x_list)
#y_i = si.splev(ipl_t, y_list)


fig = plt.figure()
plt.plot(ctrl_pts_x, ctrl_pts_y, '-og')
plt.plot(x_i, y_i, 'r')
plt.xlim([min(ctrl_pts_x) - 0.3, max(ctrl_pts_x) + 0.3])
plt.ylim([min(ctrl_pts_y) - 0.3, max(ctrl_pts_y) + 0.3])
plt.title('Splined f(x(t), y(t))')
cr1 = plt.Circle( (ca_center[0], ca_center[1]), radius=ca_radius, fc='y', ec='none')
plt.gca().add_patch(cr1)
plt.axis('equal')
plt.show()


#=============================================================================
#import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#import matplotlib.path as Path
#import matplotlib.patches as patches

#import scipy.optimize
#import functools
#from CSubd2D import *

##------------------------------------------------------------------------------
#def eval_quadr_bezier(a, b, c, t):
#  t2 = t * t
#  mt = 1 - t
#  mt2 = mt * mt
#  return a*mt2 + b*2*mt*t + c*t2, 2*(a-b)*mt + 2*(b-c)*t

##------------------------------------------------------------------------------
#def eval_cubic_bezier(a, b, c, d, t):
#    t2 = t * t
#    t3 = t2 * t
#    mt = 1-t
#    mt2 = mt * mt
#    mt3 = mt2 * mt
#    der1, _ = eval_quadr_bezier(3.*(b-a), 3.*(c-b), 3.*(d-c), t)
#    return a*mt3 + b*3*mt2*t + c*3*mt*t2 + d*t3, der1

#def curr_bezier_curve(p0, p1, n0, n1, params):
#    dlen0 = np.abs(params[0])
#    dlen1 = np.abs(params[1])
#    der0 = np.array([n0[1], -n0[0]])    
#    der1 = np.array([n1[1], -n1[0]])    
#    b = p0 + der0 * dlen0
#    c = p1 - der1 * dlen0
#    return b, c

#def get_dist_to_circle(p, c, r):
#    dist_to_center = np.linalg.norm(p-c)
#    dist_to_circle = np.abs(dist_to_center - r)
#    return dist_to_circle

#def error(params, points, normals, ca_center, ca_radius):
#    p0 = np.array([points[0][0], points[0][1]])
#    p1 = np.array([points[1][0], points[1][1]])
#    n0 = np.array([normals[0][0], normals[0][1]])
#    n1 = np.array([normals[1][0], normals[1][1]])
#    curr_b, curr_c = curr_bezier_curve(p0, p1, n0, n1, params)
#    penalty = 0.0
#    for i in range(0, 6):
#        t = 0.2 * i
#        curr_pt, _ = eval_cubic_bezier(p0, curr_b, curr_c, p1, t)
#        penalty += get_dist_to_circle(curr_pt, ca_center, ca_radius)
#    #curr_pt = eval_cubic_bezier(p0, curr_b, curr_c, p1, 0.5)
#    #penalty += get_dist_to_circle(curr_pt, ca_center, ca_radius)
#    return penalty


#points  = [(  0.0, 0.0),
#           (  1.0, 0.0)]
#normals = [(   0.0,  1.0),
#           (  -1.0, -1.0)]
#normals[0] /= np.linalg.norm(normals[0])
#normals[1] /= np.linalg.norm(normals[1])
#p0 = np.array([points[0][0], points[0][1]])
#p1 = np.array([points[1][0], points[1][1]])
#n0 = np.array([normals[0][0], normals[0][1]])
#n1 = np.array([normals[1][0], normals[1][1]])

#_, _, ca_center, ca_radius, _, _ = \
#  circle_avg(0.5, 0.5, True, p0, p1, n0, n1)
#fun = functools.partial(error, \
#                        points=points, normals=normals, \
#                        ca_center=ca_center, ca_radius=ca_radius)
#params0 = [1., 1.]
#res = scipy.optimize.minimize(fun, params0, bounds=((0.0, None),(0.0, None)))

#if not res.success:
#    print res
    

#dlen0 = np.abs(res.x[0])
#dlen1 = np.abs(res.x[1])
#der0 = np.array([n0[1], -n0[0]])    
#der1 = np.array([n1[1], -n1[0]])   
#a = p0 
#b = p0 + der0 * dlen0
#c = p1 - der1 * dlen0
#d = p1
#res_pt, res_norm = eval_cubic_bezier(a,b,c,d, 0.5)
#res_norm = np.array([res_norm[1], -res_norm[0]])

#DEBUG = True
##DEBUG = False
#if DEBUG:
#    verts = [
#        (a[0], a[1]), # P0
#        (b[0], b[1]), # P1
#        (c[0], c[1]), # P2
#        (d[0], d[1]), # P3
#        ]
        
#    codes = [Path.Path.MOVETO,
#                Path.Path.CURVE4,
#                Path.Path.CURVE4,
#                Path.Path.CURVE4,
#                ]
        
#    crv = Path.Path(verts, codes)
        
#    fig = plt.figure()
#    dn0 = plt.Arrow( p0[0], p0[1], n0[0]*.1, n0[1]*.1, width=0.03, fc='b', ec='none' )
#    dn1 = plt.Arrow( p1[0], p1[1], n1[0]*.1, n1[1]*.1, width=0.03, fc='b', ec='none' )
#    dn2 = plt.Arrow( res_pt[0], res_pt[1], res_norm[0]*.1, res_norm[1]*.1, 
#                    width=0.05, fc='r', ec='none' )
#    plt.gca().add_patch(dn0)
#    plt.gca().add_patch(dn1)
#    plt.gca().add_patch(dn2)
#    patch = patches.PathPatch(crv, facecolor='none', lw=2)
#    plt.gca().add_patch(patch)
#    plt.axis('equal')
#    plt.axis([-1.0, 2., -1.0, 2])
#    plt.show()

